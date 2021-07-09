#include "CCFOrtho.h"


void CCFOrtho::initializeVectors() {
    aux.assign(4, vector<int>());
    curr.assign(4, vector<int>());
    target.assign(4, vector<double>());
}


CCFOrtho::CCFOrtho(Microstructure& mic_, int numAxis_) : mic(mic_) { //lê mic e assume que é referência (sem target)

    initializeVectors();

    isReference = true;

    numAxis = numAxis_;
    if (numAxis == 1) {
        rmax = mic.nx / 2;
    }
    else if (numAxis == 2) {
        rmax = min(mic.nx / 2, mic.ny / 2);
    }
    else if (numAxis == 3) {
        rmax = min(mic.nx / 2, min(mic.ny / 2, mic.nz / 2));
    }

    firstCalc();

    target[0].assign(rmax, 0);
    target[1].assign(rmax, 0);
    target[2].assign(rmax, 0);
    target[3].assign(rmax, 0);

    for (int r = 0; r < rmax; r++) {
        for (int i = 0; i < 4; i++) {
            target[i][r] = double(curr[i][r]) / (numAxis * mic.n);
        }
    }


}


CCFOrtho::CCFOrtho(Microstructure& mic_, vector<vector<double>> target_, int numAxis_) : mic(mic_) { //aceita vetor como entrada

    initializeVectors();

    isReference = false;

    numAxis = numAxis_;
    if (numAxis == 1) {
        rmax = mic.nx / 2;
    }
    else if (numAxis == 2) {
        rmax = min(mic.nx / 2, mic.ny / 2);
    }
    else if (numAxis == 3) {
        rmax = min(mic.nx / 2, min(mic.ny / 2, mic.nz / 2));
    }

    firstCalc();
    target = target_;
    currEnergy = getStartEnergy();

}

CCFOrtho::CCFOrtho(Microstructure& mic_, string targetPath, int numAxis_) : mic(mic_) { // lê função do arquivo

    initializeVectors();

    isReference = false;

    numAxis = numAxis_;
    if (numAxis == 1) {
        rmax = mic.nx / 2;
    }
    else if (numAxis == 2) {
        rmax = min(mic.nx / 2, mic.ny / 2);
    }
    else if (numAxis == 3) {
        rmax = min(mic.nx / 2, min(mic.ny / 2, mic.nz / 2));
    }

    firstCalc();
    readFile(targetPath);
    currEnergy = getStartEnergy();

}

void CCFOrtho::firstCalc() {
    curr[0].assign(rmax, 0);
    curr[1].assign(rmax, 0);
    curr[2].assign(rmax, 0);
    curr[3].assign(rmax, 0);

    for (int axis = 0; axis < numAxis; axis++) {
        for (int i = 0; i < mic.n; i++) {
            mic.setPoint0(i);
            mic.setPoint1(i);
            for (int r = 0; r < rmax; r++) {
                if (mic.phasePoint0() == 0 && mic.phasePoint1() == 0) curr[0][r]++;
                if (mic.phasePoint0() == 0 && mic.phasePoint1() == 1) curr[1][r]++;
                if (mic.phasePoint0() == 1 && mic.phasePoint1() == 0) curr[2][r]++;
                if (mic.phasePoint0() == 1 && mic.phasePoint1() == 1) curr[3][r]++;
                mic.advancePoint1(axis);
            }
        }
    }
}

double CCFOrtho::getStartEnergy() {
    double sum = 0;
    for (int r = 0; r < rmax; r++) {
        for (int i = 0; i < 4; i++) {
            sum += pow(target[i][r] - double(curr[i][r]) / (numAxis * mic.n), 2);
        }
    }
    return sum;
}

double CCFOrtho::getAuxEnergy() {
    double sum = 0;
    for (int r = 0; r < rmax; r++) {
        for (int i = 0; i < 4; i++) {
            sum += pow(target[i][r] - double(aux[i][r]) / (numAxis * mic.n), 2);
        }
    }
    return sum;
}


double CCFOrtho::energyDiff() {
    return auxEnergy - currEnergy;
}

void CCFOrtho::swap() {
    curr = aux;
    currEnergy = auxEnergy;
}


void CCFOrtho::readFile(string path) {

    ifstream file(path);
    string trash;
    string value00;
    string value01;
    string value10;
    string value11;

    string line;

    target[0].assign(rmax, 0);
    target[1].assign(rmax, 0);
    target[2].assign(rmax, 0);
    target[3].assign(rmax, 0);

    for (int i = 0; i < rmax; i++) {
        file >> trash >> value00 >> value01 >> value10 >> value11;
        target[0][i] = stof(value00);
        target[1][i] = stof(value01);
        target[2][i] = stof(value10);
        target[3][i] = stof(value11);
        getline(file, line);
    }

    file.close();
}

void CCFOrtho::update(int n0_, int n1_) {
    //p0_ é sempre da fase 1
    //p1_ é sempre da fase 2

    aux = curr;

    //troca
    mic.arr[n0_] = 1;
    mic.arr[n1_] = 0;

    for (int axis = 0; axis < numAxis; axis++) {

        mic.setPoint0(n0_);
        mic.setPoint1(n0_);

        for (int r = 1; r < rmax; r++) {
            mic.backPoint0(axis);
            mic.advancePoint1(axis);
            if (mic.phasePoint0() == 0) { //00 -> 01
                aux[0][r]--;
                aux[1][r]++;
            }
            else { //10 -> 11
                aux[2][r]--;
                aux[3][r]++;
            }

            if (mic.phasePoint1() == 0) { //00 -> 10
                aux[0][r]--;
                aux[2][r]++;
            }
            else { //01 -> 11
                aux[1][r]--;
                aux[3][r]++;
            }
        }

        mic.setPoint0(n1_);
        mic.setPoint1(n1_);

        for (int r = 1; r < rmax; r++) {
            mic.backPoint0(axis);
            mic.advancePoint1(axis);
            if (mic.phasePoint0() == 0) { //01 -> 00
                aux[1][r]--;
                aux[0][r]++;
            }
            else { //11 -> 10
                aux[3][r]--;
                aux[2][r]++;
            }

            if (mic.phasePoint1() == 0) { //10 -> 00
                aux[2][r]--;
                aux[0][r]++;
            }
            else { //11 -> 01
                aux[3][r]--;
                aux[1][r]++;
            }
        }

    }

    //destroca
    mic.arr[n0_] = 0;
    mic.arr[n1_] = 1;

    auxEnergy = getAuxEnergy();
}

double CCFOrtho::getCurrEnergy() {
    return currEnergy;
}

void CCFOrtho::writeToFile(string path_) {
    ofstream file;
    file.open(path_);
    double den = mic.n * numAxis;
    int pos = 0;

    for (int r = 0; r < rmax; r++)
    {
        file << r << "\t" << curr[0][r] / den << "\t" << curr[1][r] / den << "\t" << curr[2][r] / den << "\t" << curr[3][r] / den << "\n";
        pos++;
    }
    file.close();
}

void CCFOrtho::writeToFileNonPeriodic(string path_) {
    ofstream file;
    vector<vector<double>> curr = evalNonPeriodic();
    file.open(path_);
    int pos = 0;

    for (int r = 0; r < rmax; r++)
    {
        file << r << "\t" << curr[0][r] << "\t" << curr[1][r] << "\t" << curr[2][r] << "\t" << curr[3][r] << "\n";
        pos++;
    }
    file.close();
}

vector<vector<double>> CCFOrtho::evalNonPeriodic() {

    vector<vector<double>> nonPeriodic;
    nonPeriodic.assign(4, vector<double>());
    vector<double> den;
    for (int i = 0; i < 4; i++) nonPeriodic[i].assign(rmax, 0);
    den.assign(rmax, 0);

    int diff = 0;

    vector<int> sizes;
    sizes.push_back(mic.nx);
    sizes.push_back(mic.ny);
    sizes.push_back(mic.nz);

    for (int axis = 0; axis < numAxis; axis++) {

        for (int i = 0; i < mic.n; i++) {
            mic.setPoint0(i);
            mic.setPoint1(i);
            diff = 0;
            while (diff < rmax) {

                if (mic.phasePoint0() == 0 && mic.phasePoint1() == 0) nonPeriodic[0][diff]++;
                else if (mic.phasePoint0() == 0 && mic.phasePoint1() == 1) nonPeriodic[1][diff]++;
                else if (mic.phasePoint0() == 1 && mic.phasePoint1() == 0) nonPeriodic[2][diff]++;
                else if (mic.phasePoint0() == 1 && mic.phasePoint1() == 1) nonPeriodic[3][diff]++;

                mic.advancePoint1(axis);
                den[diff]++;
                diff++;

                if (mic.coordPoint1(axis) == mic.axisSize(axis) - 1) break;
            }
        }
    }

    for (int r = 0; r < rmax; r++) {
        for (int i = 0; i < 4; i++) {
            nonPeriodic[i][r] = nonPeriodic[i][r] / den[r];
        }
    }

    return nonPeriodic;

}
