#pragma once

#include "CorrFunction.h"
#include "Microstructure.h"

using namespace std;

class Reconstructor {
public:
    double tol;
    int tmax;
    int t = 0;
    int nfuns = 0;
    vector<CorrFunction*> funcs;
    vector<double> weights;
    Microstructure& mic;
    double E_thr = 0;
    int printEvery = 10000;

    int initial = 2000;
    int chainSize = 2000;
    double factor = 0.92;
    int pos_dE_counter = 0;
    double dE_history = 0;
    double dE_thr = 0;
    double dE = 0;
    double prob = 0;
    double temp = 1000;
    double surfaceOptStart = 0.5;

    void setPrintFreq(int freq_);

    void setSurfaceOptStart(double s_);

    void updateAll(int n0_, int n1_);

    void swapAll(int i0_, int i1_);

    double totalEnergy();



    Reconstructor(Microstructure& mic_, vector<CorrFunction*>& funcs_, vector<double> weights_, double tol_, int tmax_);


    void reconstruct();

    void setCoolingSchedule(int initial_, int chainSize_, double factor_, double prob_);


};

