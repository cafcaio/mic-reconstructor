#pragma once

#include <opencv2/opencv.hpp>
#include <string>
#include <vector>
#include <random>
#include <fstream>
#include <set>
#include <algorithm>

using namespace std;
using namespace cv;

class Microstructure {

public:
    int nx, ny, nz, n, n1, n2;
    double f = 0;
    vector<int> arr, phase0, phase1;
    vector<int> xx, yy, zz;

    vector<double> lattice;

    bool hasLattice = false;

    int np0 = 0;
    int np1 = 0;
    int xp0 = 0;
    int xp1 = 0;
    int yp0 = 0;
    int yp1 = 0;
    int zp0 = 0;
    int zp1 = 0;


    Microstructure(Mat& image);
    Microstructure(int nx_, int ny_, int nz_, double f_);

    //distância diagonal máxima na microestrutura
    int max_diagonal_distance();

    //cria malha de sítios para cálculo do histograma de distâncias entre pares
    void allocateLattice();


    //retorna índices 4-vizinhos OCUPADOS do ponto n
    vector<int> neighborsIndexes(int n);

    //retorna índices vizinhos OCUPADOS do ponto n
    vector<int> neighborsIndexes8(int n);

    //retorna true se houver pelo menos dois pixels vizinhos entre si EM TORNO do ponto i
    bool hasConnectedContour(int p);

    //calcula número de pixels de fase diferente considerando uma 4-vizinhança ou 6-vizinhança (3D)
    int freeEnergy(int n);

    //calcula número de pixels de outra fase considerando uma 8-vizinhança
    int freeEnergy8(int n);

    int axisSize(int axis_);

    int dist(int i0, int i1);

    void swap(int ind0, int ind1);

    void printToFileTecplot(string path);

    //métodos que gerenciam pontos vagantes p0 e p1 para fins de CorrFunction
    void setPoint0(int n0_);

    void setPoint1(int n1_);

    void advancePoint0(int axis);

    void backPoint0(int axis);

    void advancePoint1(int axis);

    void backPoint1(int axis);

    int phasePoint0();

    int phasePoint1();

    //retorna índice vetorizado do ponto (x,y,z)
    //retorna -1 se o índice exceder os limites da matriz
    int index(int x, int y, int z);


    int coordPoint0(int axis);

    int coordPoint1(int axis);

};


