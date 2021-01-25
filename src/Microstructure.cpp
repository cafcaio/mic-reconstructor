#include "Microstructure.h"

using namespace std;
using namespace cv;


Microstructure::Microstructure(Mat& image) { //lê da imagem diretamente

    nx = image.cols;
    ny = image.rows;
    nz = 1;

    n = nx * ny * nz;

    arr.assign(n, 0);
    xx.assign(n, 0);
    yy.assign(n, 0);
    zz.assign(n, 0);
    phase0.assign(n, 0);
    phase1.assign(n, 0);

    int ind = 0;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                xx[ind] = i;
                yy[ind] = j;
                zz[ind] = k;
                ind++;
            }
        }
    }

    n1 = 0;
    n2 = 0;
    int counter0 = 0;
    int counter1 = 0;

    //binariza e insere no arr
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            if (image.at<uint8_t>(j, i) < 127) { //0 = preto, 255=branco
                arr[j * nx + i] = 1;
                n2 = n2 + 1;
                phase1[counter1] = j * nx + i;
                counter1 += 1;
            }
            else {
                n1 = n1 + 1;
                phase0[counter0] = j * nx + i;
                counter0 += 1;
            }
        }
    }

    f = double(n2) / n;

    phase0.resize(n1);
    phase1.resize(n2);
}


Microstructure::Microstructure(int nx_, int ny_, int nz_, double f_) { //cria array vazio p/ reconstrucao
    nx = nx_;
    ny = ny_;
    nz = nz_;
    f = f_;
    n = nx * ny * nz;
    n2 = n * f;
    n1 = n - n2;
    arr.assign(n, 0);
    xx.assign(n, 0);
    yy.assign(n, 0);
    zz.assign(n, 0);
    phase0.assign(n1, 0);
    phase1.assign(n2, 0);

    int ind = 0;
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                xx[ind] = i;
                yy[ind] = j;
                zz[ind] = k;
                ind++;
            }
        }
    }

    int counter = 0;
    int curr_n = 0;

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);



    while (counter < n2) {
        curr_n = n * dist(gen);
        if (arr[curr_n] == 0) {
            phase1[counter] = curr_n;
            arr[curr_n] = 1;
            counter = counter + 1;
        }
    }

    counter = 0;
    for (int i = 0; i < n; i++) {
        if (arr[i] == 0) {
            phase0[counter] = i;
            counter = counter + 1;
        }
    }
}

//distância diagonal máxima na microestrutura
int Microstructure::max_diagonal_distance() {
    return sqrt((nx - 1) * (nx - 1) + (ny - 1) * (ny - 1) + (nz - 1) * (nz - 1)) + 1;
}



//cria malha de sítios para cálculo do histograma de distâncias entre pares
void Microstructure::allocateSqrtTable() {

    if (allocatedSqrtTable == false) {
        int a = max_diagonal_distance() * max_diagonal_distance();
        sqrtTable.assign(a, 0);

        for (int i = 0; i < a; i++) {
            sqrtTable[i] = sqrt(i);
        }

    }

    allocatedSqrtTable = true;
}

void Microstructure::allocateLattice() {

    if (allocatedLattice == false) {
        lattice.assign(max_diagonal_distance(), 0);


#pragma omp parallel
        {
            vector<int> p(max_diagonal_distance(), 0);
            int distance;

#pragma omp for schedule(guided)
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < i; j++) {
                    distance = dist(i, j);
                    p[distance] += 2;
                }
            }

#pragma omp critical
            {
                for (int j = 0; j < max_diagonal_distance(); j++) {
                    lattice[j] += p[j];
                }
            }
        }

        lattice[0] += n;

    }

    allocatedLattice = true;



}

//retorna índices dos 4-vizinhos OCUPADOS do ponto n
vector<int> Microstructure::neighborsIndexes(int n) {
    vector<int> ans = {};
    ans.reserve(6);

    int ind;

    ind = index(xx[n] + 1, yy[n], zz[n]);
    if (ind != -1 && arr[ind] == 1) ans.push_back(ind);

    ind = index(xx[n] - 1, yy[n], zz[n]);
    if (ind != -1 && arr[ind] == 1) ans.push_back(ind);

    ind = index(xx[n], yy[n] + 1, zz[n]);
    if (ind != -1 && arr[ind] == 1) ans.push_back(ind);

    ind = index(xx[n], yy[n] - 1, zz[n]);
    if (ind != -1 && arr[ind] == 1) ans.push_back(ind);

    ind = index(xx[n], yy[n], zz[n] + 1);
    if (ind != -1 && arr[ind] == 1) ans.push_back(ind);

    ind = index(xx[n], yy[n], zz[n] - 1);
    if (ind != -1 && arr[ind] == 1) ans.push_back(ind);


    return ans;
}



//retorna índices dos 8-vizinhos OCUPADOS do ponto n
vector<int> Microstructure::neighborsIndexes8(int n) {
    vector<int> ans = {};

    int ind;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                ind = index(xx[n] - 1 + i, yy[n] - 1 + j, zz[n] - 1 + k);
                if (ind != n && ind != -1 && arr[ind] == 1) { //verifica se não excede a malha e se está ocupado
                    ans.push_back(ind);
                }
            }
        }
    }

    return ans;

}


//calcula número de pixels de fase diferente considerando uma 4-vizinhança ou 6-vizinhança (3D)
int Microstructure::freeEnergy(int n) {

    int sum = 0;
    int ref = 1;

    int ind;

    ind = index(xx[n] + 1, yy[n], zz[n]);
    if (ind != -1 && arr[ind] != ref) sum++;

    ind = index(xx[n] - 1, yy[n], zz[n]);
    if (ind != -1 && arr[ind] != ref) sum++;

    ind = index(xx[n], yy[n] + 1, zz[n]);
    if (ind != -1 && arr[ind] != ref) sum++;

    ind = index(xx[n], yy[n] - 1, zz[n]);
    if (ind != -1 && arr[ind] != ref) sum++;

    ind = index(xx[n], yy[n], zz[n] + 1);
    if (ind != -1 && arr[ind] != ref) sum++;

    ind = index(xx[n], yy[n], zz[n] - 1);
    if (ind != -1 && arr[ind] != ref) sum++;


    return sum;

 

}

//calcula número de pixels de outra fase considerando uma 8-vizinhança
int Microstructure::freeEnergy8(int n) {

    int sum = 0;
    int ref = 1;

    int ind;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                ind = index(xx[n] - 1 + i, yy[n] - 1 + j, zz[n] - 1 + k);
                if (ind != n && ind != -1 && arr[ind] != ref) { //verifica se não excede a malha e se está ocupado
                    sum++;
                }
            }
        }
    }

    return sum;

}

int Microstructure::axisSize(int axis_) {
    if (axis_ == 0) {
        return nx;
    }
    else if (axis_ == 1) {
        return ny;
    }
    else {
        return nz;
    }
}

int Microstructure::dist(int i0, int i1) {

    int x = xx[i0] - xx[i1];
    int y = yy[i0] - yy[i1];
    int z = zz[i0] - zz[i1];

    int a = x * x + y * y + z * z;

    if (allocatedSqrtTable) {
        return sqrtTable[a];
    }
    else {
        return sqrt(a);
    }

}

void Microstructure::swap(int ind0, int ind1) {
    //ind0 é do vetor phase0
    //ind1 é do vetor phase1

    int n0 = phase0[ind0];
    int n1 = phase1[ind1];

    arr[n0] = 1;
    arr[n1] = 0;

    phase0[ind0] = n1;
    phase1[ind1] = n0;

}

void Microstructure::printToFileTecplot(string path) {
    ofstream file;
    file.open(path);

    //Escreve o cabeçalho para o Tecplot
    file << "TITLE='Grid'\n"
        << "VARIABLES= \"X\", \"Y\", \"Z\", \"PHASE\" \n" << "ZONE"
        << " T=\"\", I=" << nx << " J= " << ny << " K= " << nz << " F=POINT\n";

    for (int i = 0; i < n; i++) {
        file << xx[i] << "\t" << yy[i] << "\t" << zz[i] << "\t" << arr[i] << "\n";
    }

    file.close();
}

void Microstructure::printToFileParaview(string path) {
    ofstream file;
    file.open(path);

    //Escreve o cabeçalho para o ParaView
    file << "# vtk DataFile Version 2.0\n"
        << "Microestrutura\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_POINTS\n"
        << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n"
        << "ORIGIN 0 0 0\n"
        << "SPACING 1 1 1\n\n"
        << "POINT_DATA " << n << "\n"
        << "SCALARS PHASE float\n"
        << "LOOKUP_TABLE default\n";

    for (int i = 0; i < n; i++) {
        file << arr[i] << "\n";
    }

    file.close();
}

//métodos que gerenciam pontos vagantes p0 e p1 para fins de CorrFunction
void Microstructure::setPoint0(int n0_) {
    np0 = n0_;
    xp0 = xx[np0];
    yp0 = yy[np0];
    zp0 = zz[np0];
}

void Microstructure::setPoint1(int n1_) {
    np1 = n1_;
    xp1 = xx[np1];
    yp1 = yy[np1];
    zp1 = zz[np1];
}

void Microstructure::advancePoint0(int axis) {
    if (axis == 0) {
        xp0++;
        if (xp0 >= nx) xp0 -= nx;
    }
    else if (axis == 1) {
        yp0++;
        if (yp0 >= ny) yp0 -= ny;
    }
    else if (axis == 2) {
        zp0++;
        if (zp0 >= nz) zp0 -= nz;
    }
    np0 = index(xp0, yp0, zp0);
}

void Microstructure::backPoint0(int axis) {
    if (axis == 0) {
        xp0--;
        if (xp0 < 0) xp0 += nx;
    }
    else if (axis == 1) {
        yp0--;
        if (yp0 < 0) yp0 += ny;
    }
    else if (axis == 2) {
        zp0--;
        if (zp0 < 0) zp0 += nz;
    }
    np0 = index(xp0, yp0, zp0);
}

void Microstructure::advancePoint1(int axis) {
    if (axis == 0) {
        xp1++;
        if (xp1 >= nx) xp1 -= nx;
    }
    else if (axis == 1) {
        yp1++;
        if (yp1 >= ny) yp1 -= ny;
    }
    else if (axis == 2) {
        zp1++;
        if (zp1 >= nz) zp1 -= nz;
    }
    np1 = index(xp1, yp1, zp1);
}

void Microstructure::backPoint1(int axis) {
    if (axis == 0) {
        xp1--;
        if (xp1 < 0) xp1 += nx;
    }
    else if (axis == 1) {
        yp1--;
        if (yp1 < 0) yp1 += ny;
    }
    else if (axis == 2) {
        zp1--;
        if (zp1 < 0) zp1 += nz;
    }
    np1 = index(xp1, yp1, zp1);
}

int Microstructure::phasePoint0() {
    return arr[np0];
}

int Microstructure::phasePoint1() {
    return arr[np1];
}

//retorna índice vetorizado do ponto (x,y,z)
//retorna -1 se o índice exceder os limites da matriz
int Microstructure::index(int x, int y, int z) {
    if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz) {
        return -1;
    }
    else {
        return nx * ny * z + nx * y + x;
    }

}


int Microstructure::coordPoint0(int axis) {
    if (axis == 0) {
        return xp0;
    }
    else if (axis == 1) {
        return yp0;
    }
    else {
        return zp0;
    }
}

int Microstructure::coordPoint1(int axis) {
    if (axis == 0) {
        return xp1;
    }
    else if (axis == 1) {
        return yp1;
    }
    else {
        return zp1;
    }
}
