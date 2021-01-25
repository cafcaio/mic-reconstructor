#pragma once

//Classe (interface) que gerencia microestruturas
//incluindo referências e reconstruções

class CorrFunction {

public:
    virtual void update(int n0_, int n1_) = 0;
    virtual void swap() = 0;
    virtual double energyDiff() = 0;
    virtual double getCurrEnergy() = 0;
    virtual double getAuxEnergy() = 0;


};


