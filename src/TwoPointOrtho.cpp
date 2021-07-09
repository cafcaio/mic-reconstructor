#include "TwoPointOrtho.h"

TwoPointOrtho::TwoPointOrtho(Microstructure& mic_, int numAxis_) : mic(mic_) { //assume como referência
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

    target.assign(rmax, 0);

    for (int r = 0; r < rmax; r++) {
        target[r] = double(curr[r]) / (numAxis * mic.n);
    }


    currEnergy = getStartEnergy();

}

TwoPointOrtho::TwoPointOrtho(Microstructure& mic_, vector<double> target_, int numAxis_) : mic(mic_) { //lê do vetor

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





TwoPointOrtho::TwoPointOrtho(Microstructure& mic_, string targetPath, int numAxis_) : mic(mic_) {

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

void TwoPointOrtho::firstCalc() {

    curr.assign(rmax, 0);

    for (int axis = 0; axis < numAxis; axis++) {
        for (int i = 0; i < mic.n; i++) {
            mic.setPoint0(i);
            mic.setPoint1(i);
            for (int r = 0; r < rmax; r++) {
                if (mic.phasePoint0() == 1 && mic.phasePoint1() == 1) {
                    curr[r]++;
                }
                mic.advancePoint1(axis);
            }
        }
    }
}

double TwoPointOrtho::getStartEnergy() {
    double sum = 0;
    for (int r = 0; r < rmax; r++) {
        sum += pow(target[r] - double(curr[r]) / (numAxis * mic.n), 2);
    }
    return sum;
}

double TwoPointOrtho::getAuxEnergy() {
    double sum = 0;
    for (int r = 0; r < rmax; r++) {
        sum += pow(target[r] - double(aux[r]) / (numAxis * mic.n), 2);
    }
    return sum;
}


double TwoPointOrtho::energyDiff() {
    return auxEnergy - currEnergy;
}

void TwoPointOrtho::swap() {
    curr = aux;
    currEnergy = auxEnergy;
}


void TwoPointOrtho::readFile(string path) {

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

void TwoPointOrtho::update(int n0_, int n1_) {
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

            if (mic.phasePoint0() == 1) aux[r]++;
            if (mic.phasePoint1() == 1) aux[r]++;
        }

        mic.setPoint0(n1_);
        mic.setPoint1(n1_);

        for (int r = 1; r < rmax; r++) {
            mic.backPoint0(axis);
            mic.advancePoint1(axis);

            if (mic.phasePoint0() == 1) aux[r]--;
            if (mic.phasePoint1() == 1) aux[r]--;
        }

        //destroca       
        mic.arr[n0_] = 0;
        mic.arr[n1_] = 1;

        auxEnergy = getAuxEnergy();


    }
}

double TwoPointOrtho::getCurrEnergy() {
    return currEnergy;
}

void TwoPointOrtho::writeToFile(string path_) {
    ofstream file;
    double den = numAxis * mic.n;
    file.open(path_);
    int pos = 0;

    for (int r = 0; r < rmax; r++)
    {
        file << r << "\t" << curr[r] / den << "\n";
        pos++;
    }
    file.close();
}

void TwoPointOrtho::writeToFileNonPeriodic(string path_) {
    ofstream file;
    vector<double> curr = evalNonPeriodic();
    file.open(path_);
    int pos = 0;

    for (int r = 0; r < rmax; r++)
    {
        file << r << "\t" << curr[r] << "\n";
        pos++;
    }
    file.close();
}

vector<double> TwoPointOrtho::evalNonPeriodic() {

    vector<double> nonPeriodic;
    vector<double> den;
    nonPeriodic.assign(rmax, 0);
    den.assign(rmax, 0);

    int diff = 0;

    for (int axis = 0; axis < numAxis; axis++) {

        for (int i = 0; i < mic.n; i++) {
            mic.setPoint0(i);
            mic.setPoint1(i);
            diff = 0;
            while (diff < rmax) {
                if (mic.phasePoint0() == 1 && mic.phasePoint1() == 1) nonPeriodic[diff]++;

                mic.advancePoint1(axis);
                den[diff]++;
                diff++;

                if (mic.coordPoint1(axis) == mic.axisSize(axis) - 1) break;

            }
        }
    }

    for (int r = 0; r < rmax; r++) {
        nonPeriodic[r] = nonPeriodic[r] / den[r];
    }

    return nonPeriodic;
}



