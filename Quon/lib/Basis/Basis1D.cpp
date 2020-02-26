#include "Basis1D.h"
#include "../Operator/Hamiltonian1D.h"

namespace quon {
    double Basis1D::getXPos(int const i) const {
        double xmin = -(this->N / 2)*this->dh;

        if (this->N % 2 == 0) {
            return xmin + this->dh * (i + 0.5);
        }
        else
        {
            return xmin + this->dh * i;
        }
    }

    void Basis1D::operator <<= (std::complex<double>* _Phi) {
        for (int i = 0; i < this->N; i++) {
            this->Phi[i] += _Phi[i];
        }
        if (this->Phi != _Phi) {
            delete[] _Phi;
            _Phi = nullptr;
        }
    }

    void Basis1D::operator = (std::complex<double>* _Phi) {
        if (this->Phi != _Phi) {
            memcpy(this->Phi, _Phi, this->N * sizeof(_Phi[0]));
            delete[] _Phi;
            _Phi = nullptr;
        }
    }

    bool Basis1D::isComparable(Basis1D basis) {
        return (basis.N == this->N && basis.dh == this->dh);
    }

    void Basis1D::operator <<= (WF_FORM mode) {
        int size = this->N;

        double x;
        switch (mode) {
        case WF_FORM::GAUSSIAN:
            for (int i = 0; i < size; i++) {
                x = this->XPos[i];
                this->Phi[i] = exp(-0.5*x*x);
            }
            normalize();
            break;
        case WF_FORM::NONE:
            for (int i = 0; i < size; i++) {
                this->Phi[i] = 0.0;
            }
            break;
        }
        
    }

    void Basis1D::operator <<= (BOUNDARY_CONDITION boundary_condition) {
        this->boundary_condition = boundary_condition;
    }

    inline void Basis1D::normalize() {
        normalize(this->Phi);
    }

    inline void Basis1D::normalize(std::complex<double>* Phi) {
        double sum = 0.0;

        sum += 0.5*norm(this->Phi[0])*dh;
        sum += 0.5*norm(this->Phi[this->N - 1])*dh;
        for (int i = 1; i < this->N - 1; i++) {
            sum += norm(this->Phi[i])*dh;
        }
        sum = sqrt(sum);
        for (int i = 0; i < this->N; i++) {
            Phi[i] /= sum;
        }
    }

     std::complex<double>* Basis1D::evolve(double _dt, Hamiltonian1D H, bool isImag = false) {
        const std::complex<double> dt = _dt * (isImag ? -iu : 1.0);

        const int size = this->N;
        const double dh = this->dh;
        const double hbar = H.getHbar();

        const std::complex<double>* Phi = this->Phi;
        std::complex<double>* Phi_next = new std::complex<double>[size];
        std::complex<double>* LastTerm = new std::complex<double>[size];
        std::complex<double>* HPhi;
        std::complex<double> NewCoeDiff;

        memcpy(Phi_next, Phi, memsize);
        memcpy(LastTerm, Phi, memsize);
        Phi = nullptr;

        for (int k = 1; k < 25; k++) {
            NewCoeDiff = (-iu * dt) / (k*hbar);
            HPhi = H(LastTerm, size, dh, this->boundary_condition, XPos);
            for (int i = 0; i < size; i++) {
                LastTerm[i] = NewCoeDiff * HPhi[i];
                Phi_next[i] = Phi_next[i] + LastTerm[i];
            }
            delete[] HPhi;
            HPhi = nullptr;
        }

        delete[] LastTerm;
        LastTerm = nullptr;

        if (isImag) {
            normalize(Phi_next);
        }

        return Phi_next;
    }
}