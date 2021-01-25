#include "Reconstructor.h"

struct Reconstructor::LogData {
    vector<double> totalEnergyLog;
    vector<double> timeLog;
    vector<int> numIterationsLog;
    vector<int> negEnergySwapsLog;
    vector<int> posEnergySwapsLog;
    vector<double> temperatureLog;
    int size = 0;

    void reserve(int s) {
        totalEnergyLog.reserve(s);
        timeLog.reserve(s);
        numIterationsLog.reserve(s);
        negEnergySwapsLog.reserve(s);
        posEnergySwapsLog.reserve(s);
        temperatureLog.reserve(s);
        size = s;
    }
};

//mutators
void Reconstructor::setConsolePrintFreq(int freq_) {
    consolePrintFreq = freq_;
}

void Reconstructor::setSurfaceOptStart(double s_) {
    surfaceOptStart = s_;
}

void Reconstructor::setWeights(vector<double> w) {
    double sum = 0;

    //normaliza pesos
    for (int i = 0; i < nfuns; i++) sum += w[i];
    for (int i = 0; i < nfuns; i++) weights[i] = w[i] / sum;

}

double Reconstructor::energyDiffAfterUpdate(int n0_, int n1_) {
    double dE = 0;
    for (int i = 0; i < nfuns; i++) {
        funcs[i]->update(n0_, n1_);
        dE += weights[i] * funcs[i]->energyDiff();
    }

    return dE;
}

void Reconstructor::swapAll(int i0_, int i1_) {
    for (int i = 0; i < nfuns; i++) {
        funcs[i]->swap();
    }
    mic.swap(i0_, i1_);
}

double Reconstructor::totalEnergy() {
    double sum = 0;
    for (int i = 0; i < nfuns; i++) {
        sum += weights[i] * funcs[i]->getCurrEnergy();
    }
    return sum;
}

void Reconstructor::setMicWriteFreq(int freq_, string fileLabel_ ) {
    micWriteFreq = freq_;
    micWriteLabel = fileLabel_;
}

Reconstructor::Reconstructor(Microstructure& mic_, vector<CorrFunction*>& funcs_, double tol_, int tmax_) : mic(mic_), funcs(funcs_) {

    tmax = tmax_;
    tol = tol_;
    nfuns = funcs.size();
    weights.assign(nfuns, 1);

    //normaliza pesos (unitários por padrão)
    setWeights(weights);

}


void Reconstructor::reconstruct() {
    //RNG
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);

    double timeElapsed = clock();

    int i0, i1, n0, n1;
    int trocasPos = 0;
    int trocasZero = 0;
    bool initialDone = false;

    int pos_dE_counter = 0;
    double dE_history = 0;
    double dE_thr = 0;
    double dE = 10;
    double temp = 1000;

    double probTeste, probThr;

    //Log variables
    LogData logFile;
    if(logFreq!=0) logFile.reserve(tmax / logFreq + 1);
    int micWriteCounter = 0;


    while (t<tmax and totalEnergy()>tol) {


        if (t < tmax * surfaceOptStart) { //faz trocas totalmente aleatórias

            i0 = mic.n1 * dist(gen);
            i1 = mic.n2 * dist(gen);
            n0 = mic.phase0[i0];
            n1 = mic.phase1[i1];
        }
        else { //prefere trocas nas superfícies pretas e brancas

            i0 = mic.n1 * dist(gen); //seleciona pixel branco
            n0 = mic.phase0[i0];

            i1 = mic.n2 * dist(gen); //seleciona pixel preto
            n1 = mic.phase1[i1];

            while (mic.freeEnergy(n1) == 0) { //força pixel superficial preto
                i1 = mic.n2 * dist(gen);
                n1 = mic.phase1[i1];
            }

        }


        dE = energyDiffAfterUpdate(n0, n1);

        if (t % chainSize == 0) {
            temp = factor * temp;
        }

        if (initialDone == false) {
            if (dE > 0.0 && pos_dE_counter < initial) {
                dE_history += dE;
                pos_dE_counter++;
                if (pos_dE_counter == initial) {
                    dE_history = (double)dE_history / initial;
                    temp = -dE_history / log(prob);
                    initialDone = true;
                }
            }
        }

        if (dE <= 0) {
            trocasZero++;
            swapAll(i0, i1);
        }
        else {
            probTeste = dist(gen);
            probThr = exp(-dE / temp);
            if (probTeste < probThr) {
                trocasPos++;
                swapAll(i0, i1);
            }
        }

        //outputs and logs

        if (t % consolePrintFreq == 0) {
            cout << "======= Resumo, t = " << t << " ==============\n";
            cout << "Energia total = " << totalEnergy() << "\n";
            cout << "Num de trocas desfavoraveis = " << trocasPos << "\n";
            cout << "Num de trocas favoraveis = " << trocasZero << " (dif: " << trocasZero - trocasPos << ")" << "\n";
            cout << "Temperatura = " << temp << "\n";
            cout << "Tempo decorrido = " << (clock() - timeElapsed) / CLOCKS_PER_SEC << " s\n\n";
        }

        if (logFreq!=0 && t % logFreq == 0) {
            logFile.totalEnergyLog.push_back(totalEnergy());
            logFile.timeLog.push_back((clock() - timeElapsed) / CLOCKS_PER_SEC);
            logFile.negEnergySwapsLog.push_back(trocasZero);
            logFile.posEnergySwapsLog.push_back(trocasPos);
            logFile.temperatureLog.push_back(temp);
            logFile.numIterationsLog.push_back(t);
        }

        if (micWriteFreq != 0 && t % micWriteFreq == 0) {
            string name = micWriteLabel + "_" + to_string(micWriteCounter)+ ".vtk";
            mic.printToFileParaview(name);
            micWriteCounter++;
        }

        t = t + 1;

    }


    //finishing outputs and logs
    cout << "======= Programa finalizado, t = " << t << " ==============\n";
    cout << "Energia total = " << totalEnergy() << "\n";
    cout << "Num de trocas desfavoraveis = " << trocasPos << "\n";
    cout << "Num de trocas favoraveis = " << trocasZero << "\n";
    cout << "Temperatura = " << temp << "\n";
    cout << "Tempo total = " << (clock() - timeElapsed) / CLOCKS_PER_SEC << " s\n\n";
    if (t == tmax) {
        cout << "Encerrado pelo critério de tempo máximo" << "\n";
    }
    else {
        cout << "Encerrado pelo critério de tolerância de energia" << "\n";
    }

    if (logFreq != 0) {
        logFile.totalEnergyLog.push_back(totalEnergy());
        logFile.timeLog.push_back((clock() - timeElapsed) / CLOCKS_PER_SEC);
        logFile.negEnergySwapsLog.push_back(trocasZero);
        logFile.posEnergySwapsLog.push_back(trocasPos);
        logFile.temperatureLog.push_back(temp);
        logFile.numIterationsLog.push_back(t);
        writeLogFile(logFile);
    }

    if (micWriteFreq != 0) {
        string name = micWriteLabel + "_" + to_string(micWriteCounter) + ".vtk";
        mic.printToFileParaview(name);
        micWriteCounter++;
    }




}



void Reconstructor::setCoolingSchedule(int initial_, int chainSize_, double factor_, double prob_) {
    initial = initial_;
    chainSize = chainSize_;
    factor = factor_;
    prob = prob_;
}

void Reconstructor::writeLogFile(LogData logFile) {
    ofstream file;
    file.open(logPath);

    file << "iteration,time,totalEnergy,temperatura,negEnergySwaps,posEnergySwaps\n";

    for (int r = 0; r < logFile.size; r++)
    {
        file << logFile.numIterationsLog[r] << ","
            << logFile.timeLog[r] << ","
            << logFile.totalEnergyLog[r] << ","
            << logFile.temperatureLog[r] << ","
            << logFile.negEnergySwapsLog[r] << ","
            << logFile.posEnergySwapsLog[r] << "\n";

    }

    file.close();
}

void Reconstructor::setLogFreq(int freq_, string path_) {
    logFreq = freq_;
    logPath = path_;
}


