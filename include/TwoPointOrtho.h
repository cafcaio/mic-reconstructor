#pragma once

#include "Microstructure.h"
#include "CorrFunction.h"

using namespace std;

class TwoPointOrtho : public CorrFunction {

public:
    vector<int> aux;
    vector<int> curr;
    vector<double> target;

    Microstructure& mic;
    int rmax, numAxis;
    double currEnergy, auxEnergy;
    bool isReference;


    //Construtores
    TwoPointOrtho(Microstructure& mic_, int numAxis_);
    TwoPointOrtho(Microstructure& mic_, vector<double> target_, int numAxis_);
    TwoPointOrtho(Microstructure& mic_, string targetPath, int numAxis_);

    //Cálculos iniciais
    void firstCalc();
    double getStartEnergy();

    //Funções comuns para todas as funções de correlação
    double getAuxEnergy() override;
    double energyDiff() override;
    void update(int n0_, int n1_) override;
    void swap() override;
    double getCurrEnergy() override;



    //Leitura e escrita noi arquivo
    void readFile(string path);
    void writeToFile(string path_);
    
    //calcula e escreve S2 considerando estrutura não-periódica
    void writeToFileNonPeriodic(string path_); 
    vector<double> evalNonPeriodic();

};

