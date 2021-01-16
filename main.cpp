//"Premature optimization is the root of all evil." - Donald Knuth

#include "Microstructure.h"
#include "CorrFunction.h"
#include "Reconstructor.h"
#include "ClusterGLP.h"
#include "TwoPointGLP.h"

using namespace cv;
using namespace std;





int main()
{

    int nx = 100;
    int ny = 100;
    int nz = 1;

    double error = 1e-7;
    int maxIterations = 1.5e06;

    int initialIterations = 3000;
    int coolingIterations = 3000;
    double coolingRate = 0.92;
    double initialProb = 0.3;
    double surfOptPerc = 0.5;

    
    Mat image = imread("fofonod80.png", IMREAD_GRAYSCALE);
    Mat cropped = image(Rect(80,80, 220, 220)).clone();

    threshold(cropped, cropped, 170, 255, THRESH_BINARY);

    resize(cropped, cropped, Size(nx, ny), 0, 0, INTER_AREA);

    //imshow("Original", image);
    //imshow("Cortada", cropped);
    //waitKey();


    clock_t start, end;

    Microstructure ref(cropped);
    ref.printToFileTecplot("referencia.dat");


    //two point -------------------------------------------------------------------------
    TwoPointGLP s2ref(ref, 2);
    s2ref.writeToFile("s2ref.dat");


    ////cluster function -------------------------------------------------------------------------
    ClusterGLP c2ref(ref, 2);
    c2ref.writeToFile("c2ref.dat");

   
    Microstructure michibrido(nx, ny, 1, ref.f); //a ser reconstruída

    TwoPointGLP s2rechibrido(michibrido, s2ref.target, 2);
    ClusterGLP c2rechibrido(michibrido, c2ref.target, 2);

    //inserir referências neste vector, usando operator&
    vector<CorrFunction*> funcs = { &s2rechibrido, &c2rechibrido };

    Reconstructor rechibrido(michibrido, funcs, { 1,1 }, error, maxIterations);
    rechibrido.setCoolingSchedule(initialIterations, coolingIterations, coolingRate, initialProb);
    rechibrido.setSurfaceOptStart(surfOptPerc);

    start = clock();
    rechibrido.reconstruct();
    end = clock();

    double elapsedTimehibrido = double(end - start) / CLOCKS_PER_SEC;

    michibrido.printToFileTecplot("reconstruidohibridos2c2.dat");
    s2rechibrido.writeToFile("s2rechibrido.dat");
    c2rechibrido.writeToFile("c2rechibrido.dat");



    cout << "Tempo de reconstrucao híbrido S2comC2: " << elapsedTimehibrido << " s" << "\n";


    TwoPointGLP s2final(michibrido, 2);
    s2final.writeToFile("s2final.dat");

    ClusterGLP c2final(michibrido, 2);
    c2final.writeToFile("c2final.dat");




    return 0;
}
