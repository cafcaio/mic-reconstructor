#include "ClusterGLP.h"

ClusterGLP::ClusterGLP(Microstructure& mic_) : mic(mic_) { //assume como referência
    isReference = true;
    firstIteration = true; //marks if it is first iteration (for vector copying purposes)


    rmax = mic.max_diagonal_distance();

    firstCalc();

    target.assign(rmax, 0);

    for (int r = 0; r < rmax; r++) {
        target[r] = pointCurr[r] / mic.lattice[r];
    }


    currEnergy = getStartEnergy();

}

ClusterGLP::ClusterGLP(Microstructure& mic_, vector<double> target_) : mic(mic_) { //lê do vetor
    isReference = false; //marks correlation function as reference (won't be updated and such)
    firstIteration = true; //marks if it is first iteration (for vector copying purposes)


    rmax = mic.max_diagonal_distance();
    //rmax = target_.size();

    firstCalc();
    target = target_;

    fillTarget();

    currEnergy = getStartEnergy();

}

ClusterGLP::ClusterGLP(Microstructure& mic_, string targetPath) : mic(mic_) {
    isReference = false; //marks correlation function as reference (won't be updated and such)
    firstIteration = true; //marks if it is first iteration (for vector copying purposes)



    rmax = mic.max_diagonal_distance();

    firstCalc();
    readFile(targetPath);

    fillTarget();

    currEnergy = getStartEnergy();

}

void ClusterGLP::fillTarget() {
    if (target.size() < mic.max_diagonal_distance()) {
        cout << "A função de correlação de cluster é pequena demais para este tamanho de microestrutura." << endl;
        cout << "Tamanho ideal para 2D: " << int(target.size() / sqrt(2)) << endl;
        cout << "Tamanho ideal para 3D: " << int(target.size() / sqrt(3)) << endl;
    }

    while (target.size() < mic.max_diagonal_distance()) {
        target.push_back(0);
    }
}

void ClusterGLP::firstCalc() {

    mic.allocateSqrtTable();
    mic.allocateLattice();


    labels.assign(mic.n, 0);
    pointCurr.assign(rmax, 0);

    clusters.reserve(mic.n2);
    clusters.push_back(vector<int>()); //0-index means unlabeled black pixel

    pointsByCluster.reserve(mic.n2);
    pointsByCluster.push_back(vector<int>()); //0-index means unlabeled black pixel



    //cluster labelling using depth-first search

    int nextClusterLabel = 0;

    for (int i = 0; i < mic.n2; i++) {
        int ind = mic.phase1[i];
        if (labels[ind] == 0) {
            nextClusterLabel++;

            clusters.push_back(vector<int>());
            clusters[nextClusterLabel].reserve(mic.n2);

            pointsByCluster.push_back(vector<int>());
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
        //cout << "val " << val << endl;
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

        clustersAux.push_back(vector<int>());
        clustersAux[label].reserve(mic.n2);

        pointsByClusterAux.push_back(vector<int>());
        pointsByClusterAux[label].assign(rmax, 0);

    }

    else {
        label = deletedClusterLabelsAux.back();
        deletedClusterLabelsAux.pop_back();

        pointsByClusterAux[label].assign(rmax, 0);



    }

    changedLabels.insert(label);
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

        //clustersAux[label].clear();
        //pointsByClusterAux[label].assign(rmax, 0);

        vector<int>().swap(clustersAux[label]);
        vector<int>().swap(pointsByClusterAux[label]);

        deletedClusterLabelsAux.push_back(label);

        changedLabels.insert(label);
    }
}

void joinClustersAux(int label) {

}

void ClusterGLP::update(int n0_, int n1_) {
    //n0_ é sempre 0 (pixel branco - índice do mic.arr)
    //n1_ é sempre 1 (pixel preto)


    restoreAux();
    
    //início tratamento n0_ (pixel branco que vai virar preto)
    mic.arr[n0_] = 1;


    vector<int> nl = neighborsLabels(n0_); //nl pode ser vazio!
    //cout << "numero de vizinhos " << nl.size() << endl;

    int newLabel;

    if (nl.empty()) {
        //então o novo ponto preto forma um novo cluster isolado
        newLabel = getNextLabel();

    }
    else {

        newLabel = nl[0];

        //cout << "entrei aqui " << endl;
        //cout << "tamanho do PBclusters com label newLabel: " << pointsByCluster[newLabel].size() << endl;
        //cout << "tamanho do PBclustersAUX com label newLabel: " << pointsByClusterAux[newLabel].size() << endl;


        //pointsByClusterAux[newLabel].assign(rmax, 0); //é necessário zerar tudo?
        changedLabels.insert(newLabel);

        //cross-check pair from now interconnected clusters

        ////programa serial
        //for (int ii = 0; ii < nl.size(); ii++) {
        //    int la = nl[ii];

        //    for (int jj = 0; jj < ii; jj++) {
        //        int lb = nl[jj];
        //        int dist, a, b;

        //        for (int j = 0; j < clustersAux[la].size(); j++) {
        //            for (int k = 0; k < clustersAux[lb].size(); k++) {
        //                a = clustersAux[la][j];
        //                b = clustersAux[lb][k];
        //                dist = mic.dist(a, b);
        //                pointAux[dist] += 2;
        //                pointsByClusterAux[newLabel][dist] += 2;
        //            }
        //        }
        //    }
        //}



        //programa paralelo



            for (int ii = 0; ii < nl.size(); ii++) {
                int la = nl[ii];

                for (int jj = 0; jj < ii; jj++) {
                    int lb = nl[jj];

#pragma omp parallel
                    {
                        vector<int> pointAuxP(rmax, 0);
                        vector<int> pointsByClusterAuxP(rmax, 0);

#pragma omp for nowait
                    for (int n = 0; n < clustersAux[la].size()*clustersAux[lb].size(); n++) {

                            int j = n / clustersAux[lb].size();
                            int k = n % clustersAux[lb].size();

                            int a = clustersAux[la][j];
                            int b = clustersAux[lb][k];
                            int dist = mic.dist(a, b);
                            pointAuxP[dist] += 2;
                            pointsByClusterAuxP[dist] += 2;
                    }

#pragma omp critical
                    {
                        for (int i = 0; i < rmax; i++) {
                            pointAux[i] += pointAuxP[i];
                            pointsByClusterAux[newLabel][i] += pointsByClusterAuxP[i];
                        }
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

                //update cluster table, appending old clusters to the new master one
                clustersAux[newLabel].insert(clustersAux[newLabel].end(), clustersAux[l].begin(), clustersAux[l].end());

                //clear old clusters
                vacateClusterAux(l);
            }
        }

        //add pairs with new black pixel
        //parallelize
        for (int i = 0; i < clustersAux[newLabel].size(); i++) { //note que nunca ocorre mic.dist(n0_, n0_)
            int b = clustersAux[newLabel][i];
            int dist = mic.dist(n0_, b);
            pointAux[dist] += 2;
            pointsByClusterAux[newLabel][dist] += 2;

        }
    }

    
    //cout << "newLabel " << newLabel << endl;
    //cout << pointsByClusterAux[newLabel].size() << endl;


    pointAux[0]++; //adicionando contribuição do novo ponto
    pointsByClusterAux[newLabel][0]++;

    labelsAux[n0_] = newLabel;
    clustersAux[newLabel].push_back(n0_);



    // ----------------- fim tratamento pixel branco que virou preto -----------------



    // ---------------- inicio tratamento pixel preto que virará branco ----------------------

    mic.arr[n1_] = 0; //muito importante para visitAux não revisitar n1_
    int lab = labelsAux[n1_]; //lab é a label que contém n1_ e todos os pontos do cluster em que n1_ estava


    if (mic.freeEnergy8(n1_)==0 ) {
        //pixel inside cluster

        changedLabels.insert(lab);


        //delete n1_ from clustersAux[lab] with back-swap trick
        auto it = find(clustersAux[lab].begin(), clustersAux[lab].end(), n1_);
        *it = move(clustersAux[lab].back());
        clustersAux[lab].pop_back();


        for (int i = 0; i < clustersAux[lab].size(); i++) { //note que nunca ocorre mic.dist(n1_, n1_)
            int b = clustersAux[lab][i];
            int dist = mic.dist(n1_, b);
            pointAux[dist] -= 2;
            pointsByClusterAux[lab][dist] -= 2;
        }

        labelsAux[n1_] = 0;
        pointAux[0]--;

        pointsByClusterAux[lab][0]--;


    }
    else {


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


        vector<int> clusterLabelsToRecompute;

        //visit neighbors again using DFS (possibly increasing the number of connected components)
        for (int i = 0; i < neighbors.size(); i++) {
            int b = neighbors[i];
            if (labelsAux[b] == 0) {
                int l = getNextLabel();
                clusterLabelsToRecompute.push_back(l);
                changedLabels.insert(l);

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

                pointAux[0] += clustersAux[l].size();
                pointsByClusterAux[l][0] += clustersAux[l].size();

#pragma omp parallel
                {

                    vector<int> pointAuxP(rmax, 0);
                    vector<int> pointsByClusterAuxP(rmax, 0);

#pragma omp for nowait
                    for (int k = 0; k < clustersAux[l].size()*(clustersAux[l].size()-1)/2; k++) {


                            int i = k / clustersAux[l].size();
                            int j = k % clustersAux[l].size();
                            if (j <= i) i = clustersAux[l].size() - i - 2, j = clustersAux[l].size() - j - 1;

                            int a = clustersAux[l][i];
                            int b = clustersAux[l][j];
                            int dist = mic.dist(a, b);
                            pointAuxP[dist] += 2;
                            pointsByClusterAuxP[dist] += 2;
                        }


#pragma omp critical
                    {
                        for (int i = 0; i < rmax; i++) {
                            pointAux[i] += pointAuxP[i];
                            pointsByClusterAux[l][i] += pointsByClusterAuxP[i];
                        }
                    }

                }


            }

        }
    }



    mic.arr[n0_] = 0;
    mic.arr[n1_] = 1;

    auxEnergy = getAuxEnergy();

}

double ClusterGLP::getCurrEnergy() {
    return currEnergy;
}

void ClusterGLP::restoreAux() {

    //pointAux = pointCurr;
    //labelsAux = labels;
    //deletedClusterLabelsAux = deletedClusterLabels;
    //pointsByClusterAux = pointsByCluster;
    //clustersAux = clusters;

    if (firstIteration) {

        pointAux = pointCurr;
        labelsAux = labels;
        deletedClusterLabelsAux = deletedClusterLabels;
        pointsByClusterAux = pointsByCluster;
        clustersAux = clusters;

        firstIteration = false;
    }

    else {

        pointAux = pointCurr;
        labelsAux = labels;
        deletedClusterLabelsAux = deletedClusterLabels;

        int s = clusters.size();
        clustersAux.resize(s);
        pointsByClusterAux.resize(s);


        for (auto f : changedLabels) {
            if (f < s) {
                clustersAux[f] = clusters[f];
                pointsByClusterAux[f] = pointsByCluster[f];
            }
        }
    }

    changedLabels.clear();

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

