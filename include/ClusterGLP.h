#pragma once

#include "Microstructure.h"
#include "CorrFunction.h"

using namespace std;
class ClusterGLP : public CorrFunction {

public:
    vector<double> pointCurr, pointAux; //current and auxiliary histogram of separation of black pixels
    vector<double> target; //predefined target cluster function
    vector<int> labels, labelsAux; //current and auxiliary cluster labels
    vector<vector<int>> clusters, clustersAux; //current and auxiliary lists of points in a given cluster
    vector<vector<int>> pointsByCluster, pointsByClusterAux; //histogram of separation of black pixels by cluster
    vector<int> deletedClusterLabels, deletedClusterLabelsAux; //stock of deleted labels for posterior use
    set<int> changedLabels; //set (i.e. non-repeating) labels for faster restore of current arrays
    Microstructure& mic; //microstructure to be operated
    int rmax; //maximum cluster size to be evaluated
    int numAxis; //number of axis considered
    double currEnergy; //current quadratic error ("energy")
    double auxEnergy; //auxiliary quadratic error ("energy")
    bool isReference = false; //marks correlation function as reference (won't be updated and such)
    bool firstIteration = true; //marks if it is first iteration (for vector copying purposes)


    //Constructors ------------------------------------------------------------------------------

    //Reference constructor
    ClusterGLP(Microstructure& mic_, int numAxis_);

    //Working microstructure constructors
    ClusterGLP(Microstructure& mic_, vector<double> target_, int numAxis_);
    ClusterGLP(Microstructure& mic_, string targetPath, int numAxis_);





    //Main functions ------------------------------------------------------------------------------
   
    //Calculates initial cluster function, initial energy and cluster labels
    void firstCalc();

    //Calculates initial energy
    double getStartEnergy();

    //Basic correlation function methods (CorrFunction class methods overrides) -------------------------
    double getAuxEnergy() override;
    double energyDiff() override;
    void swap() override;
    void update(int n0_, int n1_) override;
    double getCurrEnergy() override;

    //Returns vector with labels of neighbors (non-repeating)
    vector<int> neighborsLabels(int n);

    //Get available lable from stock
    int getNextLabel();

    //Recursively visits neighbors of point n and marks them with "currentLabel" using depth-first search (firstCalc() only)
    void visit(int n, int currentLabel);

    //Recursively visits neighbors of point n and marks them with "currentLabel" using depth-first search (for aux vectors only)
    void visitAux(int n, int currentLabel);

    //Clears the presence of a cluster in all apropriate auxiliary arrays
    void vacateClusterAux(int label);

    //Returns aux arrays to prior state (only affected labels in the last iterations), which cuts spent time
    void restoreAux();




    //Utils ------------------------------------------------------------------------------------------
    
    //Read file from "path" into target
    void readFile(string path);

    //Write cluster function to file in "path"
    void writeToFile(string path_);

    //Write cluster labels to file in "path"
    void writeLabelsToFile(string path_);

};


