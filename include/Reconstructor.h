#pragma once

#include "CorrFunction.h"
#include "Microstructure.h"

using namespace std;

class Reconstructor {
public:

    struct LogData; //struct for keeping log data for posterior file writing
    double tol; //energy tolerance
    int tmax; //maximum number of iterations
    int t = 0; //current iteration (or "time")
    int nfuns = 0; //number of correlation functions
    vector<CorrFunction*> funcs; //vector with references of correlation functions
    vector<double> weights; //weights to be applied to each correlation function (normalized in firstCalc() method)
    Microstructure& mic; //microstructure to be operated during reconstruction
    
    string micWriteLabel = "mics\\mic"; //default file name for intermediate microstructure file write
    string logPath = "RecStats.dat"; //default file name of reconstruction parameters log write

    //Cooling schedule parameters
    int initial = 2000; //initial iterations that counts toward first temperature computation
    int chainSize = 2000; //number of iterations between each cooling step
    double factor = 0.92; //factor by which temperature is multiplied in each cooling step
    double prob = 0.5; //probability of acceptance of positive energy difference swap at start of reconstruction
    double surfaceOptStart = 0.5; //starting point (in % of total iterations) of surface optimization

    //Logging frequencies
    int consolePrintFreq = 10000; //frequency to write on console
    int logFreq = 10000; //frequency to write reconstruction statistics
    int micWriteFreq = 10000; //frequency to write intermediate microstructure between start and finish of reconstruction


    //Constructor of Reconstructor object ----------------------------------------------------------
    Reconstructor(Microstructure& mic_, vector<CorrFunction*>& funcs_, double tol_, int tmax_);

    
    //Mutators ------------------------------------------------------------------------------------

    //Set and normalize weights
    void setWeights(vector<double> w);

    //Set console printing frequency
    void setConsolePrintFreq(int freq_);

    //Set reconstruction statistics printing frequency
    void setLogFreq(int freq_, string path_ = "RecStats.dat");

    //Set intermediate microstructure write-to-file frequency
    void setMicWriteFreq(int freq_, string fileLabel_ = "mics\\mic");

    //Set starting point of surface optimization (in % of total iterations)
    void setSurfaceOptStart(double s_);

    //Set cooling schedule parameters
    void setCoolingSchedule(int initial_, int chainSize_, double factor_, double prob_);

    
    //Main methods ----------------------------------------------------------------------
    
    //Update all correlation functions after swapping points n0_ and n1_ and returns energy
    double energyDiffAfterUpdate(int n0_, int n1_);

    //Swaps points i0_ and i1_ in microstructure and correlation functions
    void swapAll(int i0_, int i1_);

    //Returns total energy of all correlation functions, considering the appropriate weights
    double totalEnergy();

    //Performs reconstruction until termination
    void reconstruct();

    //Writes reconstruction statistics to file "logFile"
    void writeLogFile(LogData logFile);



};

