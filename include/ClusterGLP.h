#pragma once

#include "Microstructure.h"
#include "CorrFunction.h"

using namespace std;
class ClusterGLP : public CorrFunction {

public:
    vector<double> pointCurr, pointAux;
    vector<double> target;
    vector<int> labels, labelsAux;
    vector<vector<int>> clusters, clustersAux;
    vector<vector<int>> pointsByCluster, pointsByClusterAux;
    vector<int> deletedClusterLabels, deletedClusterLabelsAux;

    Microstructure& mic;
    int rmax, numAxis;
    double currEnergy, auxEnergy;
    bool isReference = false;


    ClusterGLP(Microstructure& mic_, int numAxis_);

    ClusterGLP(Microstructure& mic_, vector<double> target_, int numAxis_);

    ClusterGLP(Microstructure& mic_, string targetPath, int numAxis_);






    void firstCalc();


    double getStartEnergy();

    double getAuxEnergy() override;


    double energyDiff() override;

    void swap() override;


    void readFile(string path);

    int getNextLabel();


    vector<int> neighborsLabels(int n);

    void visit(int n, int currentLabel);

    void visitAux(int n, int currentLabel);

    void vacateClusterAux(int label);

    void update(int n0_, int n1_) override;

    double getCurrEnergy() override;

    void writeToFile(string path_);

    void writeLabelsToFile(string path_);

};


