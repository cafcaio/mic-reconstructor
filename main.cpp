//"Premature optimization is the root of all evil." - Donald Knuth

#include "Microstructure.h"
#include "CorrFunction.h"
#include "Reconstructor.h"
#include "ClusterGLP.h"
#include "TwoPointGLP.h"

using namespace cv;
using namespace std;


//TO DO:

//Imediatamente:
//Imprimir logs com estatísticas específicas de cada função de correlação

//Para o futuro:
//Documentar código
//Avaliar paralelização na ClusterGLP (visitAux, CCL)
//Paralelizar firstCalcs() para ambas as funções

//Testar mais microestruturas (arenitos, concretos, aços)
//Testar num computador melhor!


int main()
{
    omp_set_num_threads(omp_get_max_threads());

    //matrix dimensions
    int nx = 100;
    int ny = 100;
    int nz = 1;

    //termination criteria
    double error = 1e-7;
    int maxIterations = 2e06;

    //reconstruction parameters
    int initialIterations = 3000;
    int coolingIterations = 3000;
    double coolingRate = 0.92;
    double initialProb = 0.4;
    double surfOptPerc = 0.5;

    //Useful values:
    //(200x200, 2kk iterations): 3000,3000,0.92,0.4,0.5 (fofonod.jpg)

    
    //reference image
    Mat image = imread("imgs\\cement.png", IMREAD_GRAYSCALE);
    Mat cropped = image(Rect(0,0, 300, 300)).clone();

    threshold(cropped, cropped, 128, 255, THRESH_BINARY_INV);

    resize(cropped, cropped, Size(nx, ny), 0, 0, INTER_AREA);

    imshow("Original", image);
    imshow("Cropped", cropped);
    waitKey();
    destroyAllWindows();



    double start, end;

    Microstructure ref(cropped);
    ref.printToFileTecplot("outputs\\referencia.dat");
    //ref.printToFileParaview("outputs\\refparaview.vtk");


    //reference two point function
    TwoPointGLP s2ref(ref, 2);
    s2ref.writeToFile("outputs\\s2ref.dat");


    //reference two point cluster function
    ClusterGLP c2ref(ref, 2);
    c2ref.writeToFile("outputs\\c2ref.dat");

    //to be reconstructed microstructure and correlation functions
    Microstructure michibrido(nx, ny, 1, ref.f);

    TwoPointGLP s2rechibrido(michibrido, s2ref.target, 2);
    ClusterGLP c2rechibrido(michibrido, c2ref.target, 2);

    //use operator& for each correlation function
    vector<CorrFunction*> funcs = { &s2rechibrido, &c2rechibrido };

    //setting Reconstructor object and parameters
    Reconstructor rechibrido(michibrido, funcs, error, maxIterations);
    rechibrido.setCoolingSchedule(initialIterations, coolingIterations, coolingRate, initialProb);
    rechibrido.setWeights({ 1,1 });
    rechibrido.setConsolePrintFreq(10000);
    rechibrido.setLogFreq(10000, "outputs\\RecStats.dat");
    rechibrido.setMicWriteFreq(10000, "mics\\cement");
    rechibrido.setSurfaceOptStart(surfOptPerc);


    start = clock();
    rechibrido.reconstruct();
    end = clock();

    //write reconstructed microstructure and final correlation functions
    michibrido.printToFileTecplot("outputs\\reconstruidohibridos2c2.dat");
    s2rechibrido.writeToFile("outputs\\s2rechibrido.dat");
    c2rechibrido.writeToFile("outputs\\c2rechibrido.dat");


    cout << "Tempo total: " << (end - start) / CLOCKS_PER_SEC << " s" << "\n";


    //to check if correlation function updates are working as expected
    TwoPointGLP s2final(michibrido, 2);
    s2final.writeToFile("outputs\\s2final.dat");

    ClusterGLP c2final(michibrido, 2);
    c2final.writeToFile("outputs\\c2final.dat");

    return 0;
}