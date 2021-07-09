#include "LinealOrtho.h"

LinealOrtho::LinealOrtho(Microstructure& mic_, int numAxis_) : mic(mic_) { //assume como referência
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

LinealOrtho::LinealOrtho(Microstructure& mic_, vector<double> target_, int numAxis_) : mic(mic_) { //lê do vetor
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


LinealOrtho::LinealOrtho(Microstructure& mic_, string targetPath, int numAxis_) : mic(mic_) { //lê do arquivo
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

void LinealOrtho::firstCalc() {

    curr.assign(rmax, 0);
    int lmax = 0;

    for (int i = 0; i < mic.n; i++) {
        mic.setPoint0(i);
        if (mic.phasePoint0() == 1) {
            for (int axis = 0; axis < numAxis; axis++) {
                mic.setPoint0(i);
                lmax = 0;
                mic.advancePoint0(axis);

                while (mic.phasePoint0() == 1) {
                    lmax++;
                    mic.advancePoint0(axis);
                }

                for (int r = 0; r <= lmax; r++) curr[r]++;

            }
        }
    }
}

double LinealOrtho::getStartEnergy() {
    double sum = 0;
    for (int r = 0; r < rmax; r++) {
        sum += pow(target[r] - double(curr[r]) / (numAxis * mic.n), 2);
    }
    return sum;
}

double LinealOrtho::getAuxEnergy(){
    double sum = 0;
    for (int r = 0; r < rmax; r++) {
        sum += pow(target[r] - double(aux[r]) / (numAxis * mic.n), 2);
    }
    return sum;
}


double LinealOrtho::energyDiff(){
    return auxEnergy - currEnergy;
}

void LinealOrtho::swap() {
    curr = aux;
    currEnergy = auxEnergy;
}


void LinealOrtho::readFile(string path) {

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

void LinealOrtho::update(int n0_, int n1_){
    //p0_ é sempre da fase 0
    //p1_ é sempre da fase 1 (fase de interesse)

    aux = curr;


    //atualização do tipo 111011 -> 1111111

    int left = 0;
    int right = 0;
    int total = 0;

    for (int axis = 0; axis < numAxis; axis++) {

        mic.setPoint0(n0_);
        mic.setPoint1(n0_);
        left = 0;
        right = 0;

        mic.advancePoint1(axis);
        while (mic.phasePoint1() == 1) {
            right++;
            mic.advancePoint1(axis);
        }

        mic.backPoint0(axis);
        while (mic.phasePoint0() == 1) {
            left++;
            mic.backPoint0(axis);
        }

        total = left + right + 1; //remover corda de tamanho total

        for (int i = 0; i < left; i++) aux[i] -= (left - i);
        for (int i = 0; i < right; i++) aux[i] -= (right - i);
        for (int i = 0; i < total; i++) aux[i] += (total - i);

        //atualização do tipo 111111 -> 1110111

        mic.setPoint0(n1_);
        mic.setPoint1(n1_);
        left = 0;
        right = 0;

        mic.advancePoint1(axis);
        while (mic.phasePoint1() == 1) {
            right++;
            mic.advancePoint1(axis);
        }

        mic.backPoint0(axis);
        while (mic.phasePoint0() == 1) {
            left++;
            mic.backPoint0(axis);
        }

        total = left + right + 1; //remover corda de tamanho total

        for (int i = 0; i < total; i++) aux[i] -= (total - i);
        for (int i = 0; i < left; i++) aux[i] += (left - i);
        for (int i = 0; i < right; i++) aux[i] += (right - i);

    }

    auxEnergy = getAuxEnergy();

}

double LinealOrtho::getCurrEnergy() {
    return currEnergy;
}

void LinealOrtho::writeToFile(string path_) {
    ofstream file;
    file.open(path_);
    int pos = 0;

    for (int r = 0; r < rmax; r++)
    {
        file << r << "\t" << double(curr[r]) / (numAxis * mic.n) << "\n";
        pos++;
    }
    file.close();
}

void LinealOrtho::writeToFileNonPeriodic(string path_) {
    ofstream file;
    vector<double> vec = evalNonPeriodic();
    file.open(path_);
    int pos = 0;

    for (int r = 0; r < rmax; r++)
    {
        file << r << "\t" << vec[r] << "\n";
        pos++;
    }
    file.close();
}

vector<double> LinealOrtho::evalNonPeriodic() {

    vector<double> nonPeriodic;
    vector<int> den;
    nonPeriodic.assign(rmax, 0);
    den.assign(rmax, 0);

    for (int r = 0; r < rmax; r++) {
        if (numAxis == 1) den[r] = mic.n - mic.nx * r;
        else if (numAxis == 2) den[r] = 2 * mic.n - (mic.nx + mic.ny) * r;
        else if (numAxis == 3) den[r] = 3 * mic.n - r * (mic.nx * mic.ny + mic.nx * mic.nz + mic.ny * mic.nz);
    }


    int lmax = 0;

    for (int i = 0; i < mic.n; i++) {
        mic.setPoint0(i);
        if (mic.phasePoint0() == 1) {
            for (int axis = 0; axis < numAxis; axis++) {
                mic.setPoint0(i);
                lmax = 0;
                mic.advancePoint0(axis);

                while (mic.phasePoint0() == 1 && mic.coordPoint0(axis) != mic.axisSize(axis) - 1) {
                    lmax++;
                    mic.advancePoint0(axis);
                }

                for (int r = 0; r <= lmax; r++) nonPeriodic[r]++;

            }
        }
    }

    for (int r = 0; r < rmax; r++) {
        nonPeriodic[r] = nonPeriodic[r] / den[r];
    }

    return nonPeriodic;
}
