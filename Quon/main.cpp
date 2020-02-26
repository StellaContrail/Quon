#include "lib/Basis/Basis1D.h"
#include "lib/Operator/Hamiltonian1D.h"
#include "lib/Stream/QuFile.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <complex>

using namespace std;

double PotentialFunction(double x, const complex<double>* Phi, int i) {
    return 0.5*x*x + 10*norm(Phi[i]);
}

int main(void)
{
    clock_t start = clock();

    //------------------------------------------------------------------------
    quon::Basis1D basis(100, 0.2);
    basis <<= quon::BOUNDARY_CONDITION(quon::BC_MODE::PERIODIC);
    // TODO 任意の精度を持ったVectorが定義出来るようにする (REAL/COMPLEX)
    // TODO evolve()する際のDensityを混ぜ合わせる工程を選択可能にする
    // TODO 運動量を適応できるようにする
    // TODO WF_FORMで指定される波動関数の中心点を変更可能にする
    // TODO 更に抽象化する
    // TODO IOのWrapperを作り、file <<= basisのように書けるようにする
    basis <<= quon::WF_FORM::GAUSSIAN;
    quon::Hamiltonian1D H(1.0, 1.0, PotentialFunction);
    
    quon::QuFile file("test.txt");
    file <<= basis;
    for (int i = 0; i < 5000; i++) {
        basis = basis.evolve(0.01, H, true);
        
        if (i % 500 == 0) {
            cout << i << endl;
            // Output at each 5 iterations into a file
            //file <<= basis;
        }
        
    }
    file <<= basis;
    //------------------------------------------------------------------------

    clock_t end = clock();
    const double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    cout << "time : " << time << " s" << endl;

    cout << "Finished all procedures. Press any key to exit." << endl;
    //getchar();
    return 0;
}