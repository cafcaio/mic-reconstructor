//"Premature optimization is the root of all evil." - Donald Knuth

#include "Microstructure.h"
#include "CorrFunction.h"
#include "Reconstructor.h"
#include "ClusterGLP.h"
#include "TwoPointGLP.h"
#include <ctime>
#include <map>

using namespace std;
using namespace std::chrono;

//TO DO:

//Imediatamente:
//Imprimir logs com estatísticas específicas de cada função de correlação

//Para o futuro:
//Documentar código
//Avaliar paralelização na ClusterGLP (visitAux, CCL)
//Paralelizar firstCalcs() para ambas as funções

//Testar mais microestruturas (arenitos, concretos, aços)
//Testar num computador melhor!



// default = "YYYYMMDD-HHMMSS"
inline std::string time_stamp(const std::string& fmt = "%Y%m%d-%H%M%S")
{
    time_t timer = time(0);
    tm bt;
    localtime_s(&bt, &timer);
    char buf[64];
    std::strftime(buf, sizeof(buf), fmt.c_str(), &bt);
    string ans(buf);
    return ans;
}

//template <class T>
//class ReversibleVector {
//    vector<T> valuesOld;
//    vector<T> valuesNew;
//    set<int> changedItems;
//
//    ReversibleVector(vector<T> v) {
//        valuesOld = v;
//        valuesNew = v;
//    }
//
//    void set_value(int label, T new_value) {
//        valuesNew[label] = new_value;
//        changedItems.insert(label);
//    }
//
//    T get_value_old(int label) {
//        return valuesOld[label];
//    }
//
//    T get_value_new(int label) {
//        return valuesNew[label];
//    }
//
//
//
//
//
//};



int main()
{

    auto start = high_resolution_clock::now();

    int numThreads = omp_get_max_threads();
    omp_set_num_threads(numThreads);

    //entrada da microestrutura de referência
    Microstructure ref("testefofonodserver130.vtk");
    ref.printToFileParaview("outputs\\referencia.vtk");

    //reference two point function
    TwoPointGLP s2_ref(ref);
    ClusterGLP c2_ref(ref);

    s2_ref.writeToFile("outputs\\s2ref.dat");
    c2_ref.writeToFile("outputs\\c2ref.dat");

    
    //matrix dimensions
    int nx = 50;
    int ny = 50;
    int nz = 50;

    //termination criteria
    double error = 1e-7;
    int maxIterations = 0.5e06;

    //reconstruction parameters
    int initialIterations = 3000;
    int coolingIterations = 3000;
    double coolingRate = 0.90;
    double initialProb = 0.2;
    double surfOptPerc = 0.5;

    string timeDate = time_stamp();

    string recLabel = "fofonodserver";
    string recName = timeDate + "-" + recLabel + "-";

 

    

    //Useful values:
    //(200x200, 2kk iterations): 3000,3000,0.92,0.4,0.5 (fofonod.jpg)

    
    //to be reconstructed microstructure and correlation functions
    Microstructure mic(nx, ny, nz, ref.f);
    
    TwoPointGLP s2(mic, s2_ref.target);
    ClusterGLP c2(mic, c2_ref.target);

    //use operator& for each correlation function
    vector<CorrFunction*> funcs;
    funcs.push_back(&s2);
    funcs.push_back(&c2);

    //setting Reconstructor object and parameters
    Reconstructor rec(mic, funcs, error, maxIterations);
    rec.setCoolingSchedule(initialIterations, coolingIterations, coolingRate, initialProb);
    vector<double> w;
    w.push_back(1.0); //peso da funcao 1
    w.push_back(1.0); //peso da funcao 2

    rec.setWeights(w);
    rec.setConsolePrintFreq(10000);
    rec.setLogFreq(10000, "outputs\\RecStats.dat");

    string recLabel1 = "mics\\saida" + timeDate;
    rec.setMicWriteFreq(10000, recLabel1);
    rec.setSurfaceOptStart(surfOptPerc);


    rec.reconstruct();


    //write reconstructed microstructure and final correlation functions
    mic.printToFileParaview("outputs\\reconstructed.vtk");
    s2.writeToFile("outputs\\s2final.dat");
    c2.writeToFile("outputs\\c2final.dat");


    auto stop = high_resolution_clock::now();
    double duration = duration_cast<milliseconds>(stop - start).count() / 1000.0;
    cout << "Tempo total: " << duration << " s" << "\n";

    ClusterGLP debugc2(mic);
    debugc2.writeToFile("outputs\\debugc2.dat");


    return 0;
}