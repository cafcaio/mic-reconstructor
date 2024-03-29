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
    int rmax;
    double currEnergy, auxEnergy;
    bool isReference;


    //Construtores
    TwoPointGLP(Microstructure& mic_);
    TwoPointGLP(Microstructure& mic_, vector<double> target_);
    TwoPointGLP(Microstructure& mic_, string targetPath);

    //C�lculos iniciais
    void firstCalc();
    double getStartEnergy();

    //Fun��es comuns para todas as fun��es de correla��o
    double getAuxEnergy() override;
    double energyDiff() override;
    void update(int n0_, int n1_) override;
    double getCurrEnergy() override;
    void swap() override;

    //Leitura e escrita de arquivo
    void readFile(string path);
    void writeToFile(string path_);

};

