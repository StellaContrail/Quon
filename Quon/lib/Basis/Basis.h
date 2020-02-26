#ifndef _QUON_Basis_HEADER_
#define _QUON_Basis_HEADER_

#include "BoundaryCondition.h"
#include <complex>

namespace quon {
    enum class WF_FORM {
        NONE,
        GAUSSIAN
    };

    class Basis {
    protected:
        const int N;
        const double dh;
        BOUNDARY_CONDITION boundary_condition;
        std::complex<double> *Phi = nullptr;
        const int memsize;

        const std::complex<double> iu = std::complex<double>(0.0, 1.0);
        Basis(int const _N, double const _dh, BOUNDARY_CONDITION const _bc) : N(_N), dh(_dh), boundary_condition(_bc), memsize(sizeof(std::complex<double>) * _N) {}
        Basis(int const _N, double const _dh) : N(_N), dh(_dh), memsize(sizeof(std::complex<double>) * _N) {}
        
        virtual ~Basis() {
            if (Phi != nullptr) {
                delete[] Phi;
                Phi = nullptr;
            }
        }

        virtual void operator <<= (std::complex<double>*) = 0;
        virtual void operator = (std::complex<double>*) = 0;
        virtual void operator <<= (WF_FORM) = 0;
        virtual void operator <<= (BOUNDARY_CONDITION) = 0;
        virtual double getXPos(int const) const = 0;

        std::complex<double> inline *getWaveFunction() const { return Phi; }

        BOUNDARY_CONDITION inline getBoundaryCondition() const { return boundary_condition; }
    };
}

#endif