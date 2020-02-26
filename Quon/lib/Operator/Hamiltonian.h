#ifndef _QUON_Hamiltonian_HEADER_
#define _QUON_Hamiltonian_HEADER_

namespace quon {
    class Hamiltonian {
    protected:
        const double mass;
        const double hbar;

        Hamiltonian(double const _mass, double const _hbar) : mass(_mass), hbar(_hbar) {};
        virtual ~Hamiltonian() {};
    };
}

#endif