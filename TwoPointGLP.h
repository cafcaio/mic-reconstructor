#pragma once

#include "Microstructure.h"
#include "CorrFunction.h"

using namespace std;


class TwoPointGLP : public CorrFunction {

public:
    vector<double> pointCurr;
    vector<double> pointAux;
    vector<double> target;

    Microstructure& mic;
    int rmax, numAxis;
    double currEnergy, auxEnergy;
    bool isReference = false;


    //Construtores
    TwoPointGLP(Microstructure& mic_, int numAxis_);
    TwoPointGLP(Microstructure& mic_, vector<double> target_, int numAxis_);
    TwoPointGLP(Microstructure& mic_, string targetPath, int numAxis_);

    //Cálculos iniciais
    void firstCalc();
    double getStartEnergy();

    //Funções comuns para todas as funções de correlação
    double getAuxEnergy() override;
    double energyDiff() override;
    void update(int n0_, int n1_) override;
    double getCurrEnergy() override;
    void swap() override;

    //Leitura e escrita de arquivo
    void readFile(string path);
    void writeToFile(string path_);

};

