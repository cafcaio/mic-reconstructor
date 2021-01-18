#include "Reconstructor.h"

void Reconstructor::setPrintFreq(int freq_) {
    printEvery = freq_;
}

void Reconstructor::setSurfaceOptStart(double s_) {
    surfaceOptStart = s_;
}

void Reconstructor::updateAll(int n0_, int n1_) {
    dE = 0;
    for (int i = 0; i < nfuns; i++) {
        funcs[i]->update(n0_, n1_);
        dE += weights[i] * funcs[i]->energyDiff();
    }
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



Reconstructor::Reconstructor(Microstructure& mic_, vector<CorrFunction*>& funcs_, vector<double> weights_, double tol_, int tmax_) : mic(mic_), funcs(funcs_) {

    weights = weights_;

    double sum = 0;

    //normaliza pesos
    for (int i = 0; i < weights.size(); i++) sum += weights[i];
    for (int i = 0; i < weights.size(); i++) weights[i] = weights[i] / sum;


    tmax = tmax_;
    tol = tol_;
    nfuns = funcs.size();
    dE = 10;
}



void Reconstructor::reconstruct() {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);

    double timeElapsed = clock();

    int i0, i1, n0, n1;
    int trocasPos = 0;
    int trocasZero = 0;
    bool initialDone = false;

    double probTeste, probThr;

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

            //while (mic.freeEnergy(n0) == 0) { //força pixel superficial branco
            //    i0 = mic.n1 * dist(gen);
            //    n0 = mic.phase0[i0];
            //}

            i1 = mic.n2 * dist(gen); //seleciona pixel preto
            n1 = mic.phase1[i1];

            while (mic.freeEnergy8(n1) == 0) { //força pixel superficial preto
                i1 = mic.n2 * dist(gen);
                n1 = mic.phase1[i1];
            }

        }


        updateAll(n0, n1);

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


        if (t % printEvery == 0) {
            cout << "======= Resumo, t = " << t << " ==============\n";
            cout << "Energia total = " << totalEnergy() << "\n";
            cout << "Num de trocas desfavoraveis = " << trocasPos << "\n";
            cout << "Num de trocas favoraveis = " << trocasZero << " (dif: " << trocasZero - trocasPos << ")" << "\n";
            cout << "Temperatura = " << temp << "\n";
            cout << "Tempo decorrido = " << (clock() - timeElapsed) / CLOCKS_PER_SEC << " s\n\n";
        }


        t = t + 1;



    }
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



}

void Reconstructor::setCoolingSchedule(int initial_, int chainSize_, double factor_, double prob_) {
    initial = initial_;
    chainSize = chainSize_;
    factor = factor_;
    prob = prob_;
}
