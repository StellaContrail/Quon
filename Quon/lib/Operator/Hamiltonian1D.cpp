#include "Hamiltonian1D.h"
#include "../Basis/Basis1D.h"

namespace quon {
    std::complex<double>* Hamiltonian1D::operator()(Basis1D* basis) {
        int size = basis->getSize();
        double dh = basis->getStep();

        BOUNDARY_CONDITION bc = basis->getBoundaryCondition();
        std::complex<double>* Phi = basis->getPhi();
        
        double* xpos = new double[size];
        for (int i = 0; i < size; i++) {
            xpos[i] = basis->getXPos(i);
        }

        return (*this)(Phi, size, dh, bc, xpos);
    }

    std::complex<double>* Hamiltonian1D::operator()(std::complex<double>* Phi, int size, double dh, BOUNDARY_CONDITION bc, double* xpos) {
        std::complex<double>* HPhi = new std::complex<double>[size];
        switch (bc.getMode()) {
        case BC_MODE::DIRICHLET:
            Phi[0] = bc.yleft();
            Phi[size - 1] = bc.yright();

            HPhi[0] = 0.0;
            HPhi[1] = (Phi[2] - 2.0*Phi[1] + Phi[0]) / (dh*dh);
            HPhi[size - 2] = (Phi[size - 1] - 2.0*Phi[size - 2] + Phi[size - 3]) / (dh*dh);
            HPhi[size - 1] = 0.0;
            break;
        case BC_MODE::NEUMANN:
            Phi[0] = Phi[1] - bc.yleft() * dh;
            Phi[size - 1] = Phi[size - 2] + bc.yright() * dh;

            HPhi[0] = 0.0;
            HPhi[1] = (Phi[2] - 2.0*Phi[1] + Phi[0]) / (dh*dh);
            HPhi[size - 2] = (Phi[size - 1] - 2.0*Phi[size - 2] + Phi[size - 3]) / (dh*dh);
            HPhi[size - 1] = 0.0;
            break;
        case BC_MODE::PERIODIC:
            HPhi[0] = (-Phi[2] + 16.0*Phi[1] - 30.0*Phi[0] + 16.0*Phi[size - 1] - Phi[size - 2]) / (12.0*dh*dh);
            HPhi[1] = (-Phi[3] + 16.0*Phi[2] - 30.0*Phi[1] + 16.0*Phi[0] - Phi[size - 1]) / (12.0*dh*dh);
            HPhi[size - 2] = (-Phi[0] + 16.0*Phi[size - 1] - 30.0*Phi[size - 2] + 16.0*Phi[size - 3] - Phi[size - 4]) / (12.0*dh*dh);
            HPhi[size - 1] = (-Phi[1] + 16.0*Phi[0] - 30.0*Phi[size - 1] + 16.0*Phi[size - 2] - Phi[size - 3]) / (12.0*dh*dh);
            break;
        }

        for (int i = 2; i < size - 2; i++) {
            HPhi[i] = (-Phi[i + 2] + 16.0*Phi[i + 1] - 30.0*Phi[i] + 16.0*Phi[i - 1] - Phi[i - 2]) / (12.0*dh*dh);
        }

        double(*V)(double, const std::complex<double>*, int) = this->potential_function();
        double KineticCoefficient = this->getKineticCoe();
        for (int i = 0; i < size; i++) {
            HPhi[i] *= KineticCoefficient;
            HPhi[i] += V(xpos[i], Phi, i)*Phi[i];
        }

        V = nullptr;

        return HPhi;
    }
}