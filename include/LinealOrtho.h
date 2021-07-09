#pragma once

#include "Microstructure.h"
#include "CorrFunction.h"

using namespace std;

class LinealOrtho : public CorrFunction {

public:
    vector<int> aux;
    vector<int> curr;
    vector<double> target;

    Microstructure& mic;
    int rmax, numAxis;
    double currEnergy, auxEnergy;
    bool isReference;


    //Construtores
    LinealOrtho(Microstructure& mic_, int numAxis_);
    LinealOrtho(Microstructure& mic_, vector<double> target_, int numAxis_);
    LinealOrtho(Microstructure& mic_, string targetPath, int numAxis_);

    //Cálculos iniciais
    void firstCalc();
    double getStartEnergy();

    //Funções comuns a todas as funções de correlação
    double getAuxEnergy() override;
    double energyDiff() override;
    void update(int n0_, int n1_) override;
    double getCurrEnergy() override;
    void swap();

    //Leitura e escrita no arquivo
    void readFile(string path);
    void writeToFile(string path_);

    //Calcula e escreve L2 considerando que a estrutura não é periódica
    void writeToFileNonPeriodic(string path_);
    vector<double> evalNonPeriodic();

};