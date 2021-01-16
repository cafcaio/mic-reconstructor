#include "ClusterGLP.h"

//imprime vector<int> em linha
void printVector(vector<int> vector) {
    for (int i = 0; i < vector.size(); i++) {
        cout << vector[i] << "\t";
    }
    cout << "\n";
}

ClusterGLP::ClusterGLP(Microstructure& mic_, int numAxis_) : mic(mic_) { //assume como referência
    isReference = true;

    numAxis = numAxis_;

    rmax = mic.max_diagonal_distance();

    firstCalc();

    target.assign(rmax, 0);

    for (int r = 0; r < rmax; r++) {
        target[r] = pointCurr[r] / mic.lattice[r];
    }


    currEnergy = getStartEnergy();

}

ClusterGLP::ClusterGLP(Microstructure& mic_, vector<double> target_, int numAxis_) : mic(mic_) { //lê do vetor
    numAxis = numAxis_;

    rmax = mic.max_diagonal_distance();

    firstCalc();
    target = target_;
    currEnergy = getStartEnergy();

}

ClusterGLP::ClusterGLP(Microstructure& mic_, string targetPath, int numAxis_) : mic(mic_) {
    numAxis = numAxis_;

    rmax = mic.max_diagonal_distance();

    firstCalc();
    readFile(targetPath);
    currEnergy = getStartEnergy();

}

void ClusterGLP::firstCalc() {

    mic.allocateLattice();

    labels.assign(mic.n, 0);
    pointCurr.assign(rmax, 0);

    clusters.reserve(mic.n2);
    clusters.push_back({}); //0-index means unlabeled black pixel

    pointsByCluster.reserve(mic.n2);
    pointsByCluster.push_back({}); //0-index means unlabeled black pixel


    //cluster labelling using depth-first search

    int nextClusterLabel = 0;

    for (int i = 0; i < mic.n2; i++) {
        int ind = mic.phase1[i];
        if (labels[ind] == 0) {
            nextClusterLabel++;

            clusters.push_back({});
            clusters[nextClusterLabel].reserve(mic.n2);

            pointsByCluster.push_back({});
            pointsByCluster[nextClusterLabel].assign(rmax, 0);

            visit(ind, nextClusterLabel);
        }
    }



    // cluster functions computation
    int dist, a, b;

    for (int i = 1; i < clusters.size(); i++) {

        for (int j = 0; j < clusters[i].size(); j++) {
            a = clusters[i][j];
            for (int k = 0; k < j; k++) {
                b = clusters[i][k];
                dist = mic.dist(a, b);
                pointCurr[dist] += 2;
                pointsByCluster[i][dist] += 2;

            }

        }

        pointCurr[0] += clusters[i].size();
        pointsByCluster[i][0] += clusters[i].size();
    }
}


double ClusterGLP::getStartEnergy() {
    double sum = 0;
    double val;
    for (int r = 0; r < rmax; r++) {
        val = (target[r] - pointCurr[r] / mic.lattice[r]);
        sum += val * val;
    }
    return sum;
}

double ClusterGLP::getAuxEnergy() {
    double sum = 0;
    double val;
    for (int r = 0; r < rmax; r++) {
        val = (target[r] - pointAux[r] / mic.lattice[r]);
        sum += val * val;
    }
    return sum;
}


double ClusterGLP::energyDiff() {
    return auxEnergy - currEnergy;
}

void ClusterGLP::swap() {

    // current <- aux

    pointCurr = pointAux;
    currEnergy = auxEnergy;
    clusters = clustersAux;
    labels = labelsAux;
    deletedClusterLabels = deletedClusterLabelsAux;
    pointsByCluster = pointsByClusterAux;

}


void ClusterGLP::readFile(string path) {

    ifstream file(path);
    string value;

    string line;

    target.assign(rmax, 0);

    for (int r = 0; r < rmax; r++) {
        file >> value >> value;
        target[r] = stod(value);
        getline(file, line);
    }

    file.close();
}

int ClusterGLP::getNextLabel() {
    int label;
    if (deletedClusterLabelsAux.empty()) {

        label = clustersAux.size();

        clustersAux.push_back({});
        clustersAux[label].reserve(mic.n2);

        pointsByClusterAux.push_back({});
        pointsByClusterAux[label].assign(rmax, 0);

    }

    else {
        label = deletedClusterLabelsAux.back();
        deletedClusterLabelsAux.pop_back();
    }

    return label;
}


vector<int> ClusterGLP::neighborsLabels(int n) {
    vector<int> neighbors = mic.neighborsIndexes(n);
    set<int> ans;

    for (int i = 0; i < neighbors.size(); i++) {
        int ind = neighbors[i];
        ans.insert(labelsAux[ind]);
    }

    vector<int> ansv(ans.begin(), ans.end());

    return ansv;
}

void ClusterGLP::visit(int n, int currentLabel) {
    if (labels[n] != 0) return;
    labels[n] = currentLabel;
    clusters[currentLabel].push_back(n);
    vector<int> neighbors = mic.neighborsIndexes(n);
    for (int i = 0; i < neighbors.size(); i++) {
        if (labels[neighbors[i]] == 0) {
            visit(neighbors[i], currentLabel);
        }
    }
}

void ClusterGLP::visitAux(int n, int currentLabel) {
    if (labelsAux[n] != 0) return;
    labelsAux[n] = currentLabel;
    clustersAux[currentLabel].push_back(n);
    vector<int> neighbors = mic.neighborsIndexes(n); //pode ser vazio
    for (int i = 0; i < neighbors.size(); i++) {
        if (labelsAux[neighbors[i]] == 0) { //ocupado e não visitado
            visitAux(neighbors[i], currentLabel);
        }
    }
}

void ClusterGLP::vacateClusterAux(int label) {
    if (!clustersAux[label].empty()) {
        clustersAux[label].clear();
        pointsByClusterAux[label].assign(rmax, 0);
        deletedClusterLabelsAux.push_back(label);
    }
}

void ClusterGLP::update(int n0_, int n1_) {
    //n0_ é sempre 0 (pixel branco - índice do mic.arr)
    //n1_ é sempre 1 (pixel preto)

    pointAux = pointCurr;
    clustersAux = clusters;
    labelsAux = labels;
    deletedClusterLabelsAux = deletedClusterLabels;
    pointsByClusterAux = pointsByCluster;

    //início tratamento n0_ (pixel branco que vai virar preto)
    mic.arr[n0_] = 1;


    vector<int> nl = neighborsLabels(n0_); //nl pode ser vazio!

    int newLabel;

    if (nl.empty()) {
        newLabel = getNextLabel();
    }
    else {

        newLabel = *min_element(nl.begin(), nl.end());
    }

    //cross-check pair from now interconnected clusters
    for (int ii = 0; ii < nl.size(); ii++) {
        int la = nl[ii];

        for (int jj = 0; jj < ii; jj++) {
            int lb = nl[jj];
            int dist, a, b;

            for (int j = 0; j < clustersAux[la].size(); j++) {
                a = clustersAux[la][j];

                for (int k = 0; k < clustersAux[lb].size(); k++) {
                    b = clustersAux[lb][k];
                    dist = mic.dist(a, b);
                    pointAux[dist] += 2;
                    pointsByClusterAux[newLabel][dist] += 2;
                }
            }
        }
    }

    for (int i = 0; i < nl.size(); i++) {
        int l = nl[i];
        if (l != newLabel) {
            //reassign labels of points of other clusters to main one
            for (int j = 0; j < clustersAux[l].size(); j++) {
                labelsAux[clustersAux[l][j]] = newLabel;
            }

            //copy pairs of other clusters to main one
            for (int k = 0; k < pointsByClusterAux[l].size(); k++) {
                pointsByClusterAux[newLabel][k] += pointsByClusterAux[l][k];
            }

            //update cluster table
            clustersAux[newLabel].insert(clustersAux[newLabel].end(), clustersAux[l].begin(), clustersAux[l].end());

            //clear old clusters
            vacateClusterAux(l);
        }
    }

    //add pairs with new black pixel

    for (int i = 0; i < clustersAux[newLabel].size(); i++) { //note que nunca ocorre mic.dist(n0_, n0_)
        int b = clustersAux[newLabel][i];
        int dist = mic.dist(n0_, b);
        pointAux[dist] += 2;
        pointsByClusterAux[newLabel][dist] += 2;

    }

    pointAux[0]++; //adicionando contribuição do novo ponto
    pointsByClusterAux[newLabel][0]++;

    labelsAux[n0_] = newLabel;
    clustersAux[newLabel].push_back(n0_);



    // ----------------- fim tratamento pixel branco que virou preto -----------------



    // ---------------- inicio tratamento pixel preto que virará branco ----------------------

    mic.arr[n1_] = 0; //muito importante para visitAux não revisitar n1_

    int lab = labelsAux[n1_]; //lab é a label que contém n1_ e todos os pontos do cluster em que n1_ estava

    vector<int> neighbors = mic.neighborsIndexes(n1_);

    //removing contribution of soon-to-be dismantled cluster of the general count
    for (int i = 0; i < pointsByClusterAux[lab].size(); i++) {
        pointAux[i] -= pointsByClusterAux[lab][i];
    }

    //making points of soon-to-be dismantled cluster label-less
    for (int i = 0; i < clustersAux[lab].size(); i++) {
        int a = clustersAux[lab][i];
        labelsAux[a] = 0;
    }

    vector<int> saved = pointsByClusterAux[lab];

    //remove cluster
    vacateClusterAux(lab);


    vector<int> clusterLabelsToRecompute = {};

    //visit neighbors again using DFS (possibly increasing the number of connected components)
    for (int i = 0; i < neighbors.size(); i++) {
        int b = neighbors[i];
        if (labelsAux[b] == 0) {
            int l = getNextLabel();
            clusterLabelsToRecompute.push_back(l);
            visitAux(b, l);
        }
    }

    if (clusterLabelsToRecompute.size() == 1) {

        for (int i = 0; i < pointsByClusterAux[lab].size(); i++) {
            pointsByClusterAux[lab][i] = saved[i];
            pointAux[i] += pointsByClusterAux[lab][i];
        }

        for (int i = 0; i < clustersAux[lab].size(); i++) { //note que nunca ocorre mic.dist(n1_, n1_)
            int b = clustersAux[lab][i];
            int dist = mic.dist(n1_, b);
            pointAux[dist] -= 2;
            pointsByClusterAux[lab][dist] -= 2;

        }

        pointsByClusterAux[lab][0]--;
        pointAux[0]--;

    }
    else {

        for (int ii = 0; ii < clusterLabelsToRecompute.size(); ii++) {
            int l = clusterLabelsToRecompute[ii];

            for (int i = 0; i < clustersAux[l].size(); i++) {
                int a = clustersAux[l][i];
                for (int j = 0; j < i; j++) {
                    int b = clustersAux[l][j];
                    int dist = mic.dist(a, b);
                    pointAux[dist] += 2;
                    pointsByClusterAux[l][dist] += 2;
                }
            }
            pointAux[0] += clustersAux[l].size();
            pointsByClusterAux[l][0] += clustersAux[l].size();

        }
    }



    mic.arr[n0_] = 0;
    mic.arr[n1_] = 1;

    auxEnergy = getAuxEnergy();

}

double ClusterGLP::getCurrEnergy() {
    return currEnergy;
}

void ClusterGLP::writeToFile(string path_) {
    ofstream file;
    file.open(path_);
    int pos = 0;

    for (int r = 0; r < rmax; r++)
    {
        file << r << "\t" << pointCurr[r] / mic.lattice[r] << "\n";
        pos++;
    }
    file.close();
}

void ClusterGLP::writeLabelsToFile(string path_) {
    ofstream file;
    file.open(path_);

    //Escreve o cabeçalho para o Tecplot
    file << "TITLE='Grid'\n"
        << "VARIABLES= \"X\", \"Y\", \"Z\", \"PHASE\" \n" << "ZONE"
        << " T=\"\", I=" << mic.nx << " J= " << mic.ny << " K= " << mic.nz << " F=POINT\n";

    for (int i = 0; i < mic.n; i++) {
        file << mic.xx[i] << "\t" << mic.yy[i] << "\t" << mic.zz[i] << "\t" << labels[i] << "\n";
    }

    file.close();
}