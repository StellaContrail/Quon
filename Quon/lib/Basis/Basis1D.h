#ifndef _QUON_Basis1D_HEADER_
#define _QUON_Basis1D_HEADER_

#include "Basis.h"
#include "BoundaryCondition.h"
#include "../Operator/Hamiltonian1D.h"

#include <complex>
#include <fstream>

namespace quon {
    class Basis1D : private Basis {
    private:
        double* XPos = nullptr;
    public:
        Basis1D(int _N, double _dh, BOUNDARY_CONDITION _bc) : Basis(_N, _dh, _bc) {
            this->Phi = new std::complex<double>[N];
            this->XPos = new double[N];
            for (int i = 0; i < _N; i++) {
                this->XPos[i] = this->getXPos(i);
            }
        }
        Basis1D(int _N, double _dh) : Basis(_N, _dh) {
            this->Phi = new std::complex<double>[N];
            this->XPos = new double[N];
            for (int i = 0; i < _N; i++) {
                this->XPos[i] = this->getXPos(i);
            }
        }
        Basis1D(const Basis1D& basis1d) : Basis(basis1d.N, basis1d.dh, basis1d.boundary_condition) {
            this->Phi = new std::complex<double>[this->N];
            memcpy(this->Phi, basis1d.Phi, this->memsize);
            this->XPos = new double[this->N];
            for (int i = 0; i < this->N; i++) {
                this->XPos[i] = this->getXPos(i);
            }
        }
        ~Basis1D() {
            if (XPos != nullptr) {
                delete[] XPos;
                XPos = nullptr;
            }
        }

        inline int getSize() { return N; }
        inline double getStep() { return dh; }
        inline std::complex<double>* getPhi() { return Phi; }
        inline BOUNDARY_CONDITION getBoundaryCondition() { return this->boundary_condition; }

        double getXPos(int const) const;

        void operator <<= (std::complex<double>* Phi);

        void operator = (std::complex<double>* _Phi);

        void operator <<= (WF_FORM mode);

        void operator <<= (BOUNDARY_CONDITION boundary_condition);
        
        bool isComparable(Basis1D basis);

        inline void normalize();
        inline void normalize(std::complex<double>* Phi);

        std::complex<double>* evolve(double _dt, Hamiltonian1D H, bool isImag);
    };
}
#endif