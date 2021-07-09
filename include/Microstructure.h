#pragma once

#include <omp.h>
#include <string>
#include <vector>
#include <random>
#include <fstream>
#include <set>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;


class Microstructure {

public:
    
    int nx, ny, nz; //microstructure dimensions
    int n; //total number of pixels
    int n1; //number of white pixels
    int n2;  //number of black pixels
    double f; //volume fraction of black pixels
    vector<int> arr; //array of phases
    vector<int> phase0; //indices of all white pixels
    vector<int> phase1; //indexes of all black pixels
    vector<int> xx, yy, zz; //coordinates of all n pixels

    vector<double> lattice; //lattice sites separation histogram
    vector<int> sqrtTable; //cached square root table (for distance computation)
    bool allocatedLattice;
    bool allocatedSqrtTable;

    //vagant points coords and indexes
    //removal intended

    int np0, np1;
    int xp0, xp1;
    int yp0, yp1;
    int zp0, zp1;

    //initialize int and bool variables
    void initialize();


    //Constructors ------------------------------------------------------------------------------
    
    //Create matrix based on ParaView file (reference microstructure)
    Microstructure(string path_);
    
    //Creates a ***random*** matrix with specified dimensions and volume fraction
    Microstructure(int nx_, int ny_, int nz_, double f_);







    //Main methods ------------------------------------------------------------------------------

    //Returns indices of all OCCUPIED 4-neighbors around pixel n
    vector<int> neighborsIndexes(int n);

    //Returns indices of all OCCUPIED 8-neighbors around pixel n
    vector<int> neighborsIndexes8(int n);

    //Returns number of UNOCCUPIED 4-neighbors around pixel n (8-neighbors in 3D)
    int freeEnergy(int n);

    //Returns number of UNOCCUPIED 8-neighbors around pixel n (26-neighbors in 3D)
    int freeEnergy8(int n);

    //Returns distance between pixels i0 and i1
    int dist(int i0, int i1);

    //Swaps pixels ind0 and ind1, where ind0 is an index of vector phase0 and ind1 is an index of vector phase1
    void swap(int ind0, int ind1);




    //Utils ---------------------------------------------------------------------------------------

    //Calculates maximum distance of a pixel pair in the matrix
    int max_diagonal_distance();

    //Allocates lattice site pair distance histogram
    void allocateLattice();

    //Allocates square root table for fast distance calculations
    void allocateSqrtTable();

    //returns vectorized index of point (x,y,z)
    //returns -1 if point exceeds matrix limits
    int index(int x, int y, int z);

    //Returns axis axis_ size (axis are 0-indexed starting with x-axis)
    int axisSize(int axis_);




    //Outputs --------------------------------

    //Write Tecplot-formatted microstructure current phase state in path directory  (.dat file recommended)
    void printToFileTecplot(string path);

    //Write Paraview-formatted microstructure current phase state in path directory  (.dat file recommended)
    void printToFileParaview(string path);



    //Points accessors and mutators ----------------------------
    //removal intended

    void setPoint0(int n0_);
    void setPoint1(int n1_);
    void advancePoint0(int axis);
    void backPoint0(int axis);
    void advancePoint1(int axis);
    void backPoint1(int axis);
    int phasePoint0();
    int phasePoint1();
    int coordPoint0(int axis);
    int coordPoint1(int axis);

};


