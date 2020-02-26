#ifndef _QUON_Operator1D_HEADER_
#define _QUON_Operator1D_HEADER_

#include "Hamiltonian.h"
#include "../Basis/BoundaryCondition.h"
#include <complex>

namespace quon {
    class Basis1D;
    class Hamiltonian1D : public Hamiltonian {
    protected:
        double(*potential_function_1d)(double, const std::complex<double>*, int) = nullptr;

    public:
        Hamiltonian1D(double const _mass, double const _hbar, double(*_potential_function_1d)(double, const std::complex<double>*, int))
            : Hamiltonian(_mass, _hbar), potential_function_1d(_potential_function_1d) {}
        ~Hamiltonian1D() {
            potential_function_1d = nullptr;
        }

        inline double getKineticCoe() const { return -0.5*(this->hbar*this->hbar / this->mass); }
        inline double getHbar() const { return this->hbar; }

        typedef double(*func_ptr_1d)(double, const std::complex<double>*, int);
        inline func_ptr_1d potential_function() const { return this->potential_function_1d; }

        std::complex<double>* operator()(std::complex<double>* Phi, int size, double dh, BOUNDARY_CONDITION bc, double* xpos);
        std::complex<double>* operator()(Basis1D* basis);
    };
}

#endif