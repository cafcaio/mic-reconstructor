#include "TwoPointGLP.h"

TwoPointGLP::TwoPointGLP(Microstructure& mic_, int numAxis_) : mic(mic_) { //assume como referência
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

TwoPointGLP::TwoPointGLP(Microstructure& mic_, vector<double> target_, int numAxis_) : mic(mic_) { //lê do vetor
    numAxis = numAxis_;


    rmax = mic.max_diagonal_distance();

    firstCalc();
    target = target_;
    currEnergy = getStartEnergy();

}





TwoPointGLP::TwoPointGLP(Microstructure& mic_, string targetPath, int numAxis_) : mic(mic_) {
    numAxis = numAxis_;

    rmax = mic.max_diagonal_distance();

    firstCalc();
    readFile(targetPath);
    currEnergy = getStartEnergy();

}



void TwoPointGLP::firstCalc() {

    pointCurr.assign(rmax, 0);
    mic.allocateLattice();

    int dist;


    for (int i = 0; i < mic.n2; i++) {
        //primeiro ponto
        for (int j = 0; j < mic.n2; j++) {
            //todos os outros segundos pontos
            dist = mic.dist(mic.phase1[i], mic.phase1[j]); //phase1 = só pontos pretos
            pointCurr[dist]++;
        }
    }
}


double TwoPointGLP::getStartEnergy() {
    double sum = 0;
    double val;
    for (int r = 0; r < rmax; r++) {
        val = (target[r] - pointCurr[r] / mic.lattice[r]);
        sum += val * val;
    }
    return sum;
}

double TwoPointGLP::getAuxEnergy() {
    double sum = 0;
    double val;
    for (int r = 0; r < rmax; r++) {
        val = (target[r] - pointAux[r] / mic.lattice[r]);
        sum += val * val;
    }
    return sum;
}


double TwoPointGLP::energyDiff() {
    return auxEnergy - currEnergy;
}

void TwoPointGLP::swap() {
    pointCurr = pointAux;
    currEnergy = auxEnergy;
}


void TwoPointGLP::readFile(string path) {

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

void TwoPointGLP::update(int n0_, int n1_) {
    //n0_ é sempre da fase 1 (índice do array principal de Microstructure)
    //n1_ é sempre da fase 2

    pointAux = pointCurr;

    int dist;

    for (int i = 0; i < mic.n2; i++) {
        dist = mic.dist(n1_, mic.phase1[i]);
        if (dist != 0) pointAux[dist] = pointAux[dist] - 2;
    }


    for (int i = 0; i < mic.n2; i++) {
        dist = mic.dist(n0_, mic.phase1[i]);
        if (mic.phase1[i] != n1_)  pointAux[dist] = pointAux[dist] + 2;
    }


    auxEnergy = getAuxEnergy();

}

double TwoPointGLP::getCurrEnergy() {
    return currEnergy;
}

void TwoPointGLP::writeToFile(string path_) {
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


