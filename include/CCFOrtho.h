#pragma once

#include "Microstructure.h"
#include "CorrFunction.h"

using namespace std;

class CCFOrtho : public CorrFunction {

    //convencao 
    //target[0][i] = par 0 0
    //target[1][i] = par 0 1
    //target[2][i] = par 1 0
    //target[3][i] = par 1 1

public:
    vector<vector<int>> aux;
    vector<vector<int>> curr;
    vector<vector<double>> target;
    Microstructure& mic;
    int rmax, numAxis;
    double currEnergy, auxEnergy;
    bool isReference;


    //initialize vectors
    void initializeVectors();

    //Construtures
    CCFOrtho(Microstructure& mic_, int numAxis_);
    CCFOrtho(Microstructure& mic_, vector<vector<double>> target_, int numAxis_);
    CCFOrtho(Microstructure& mic_, string targetPath, int numAxis_);

    //Cálculos iniciais
    void firstCalc();
    double getStartEnergy();

    //Funções comuns a todas as funções de correlação
    void swap() override;
    void update(int n0_, int n1_) override;
    double energyDiff() override;
    double getCurrEnergy() override;
    double getAuxEnergy() override;


    //Leitura e escrita no arquivo
    void readFile(string path);
    void writeToFile(string path_);
    void writeToFileNonPeriodic(string path_);
    vector<vector<double>> evalNonPeriodic();

};

