<header style="text-align: center;">
    <p style="font-size: 35px;"> Quon</p>
    <span style="font-size: 15px;"> Quantum Mechanics Library with Fast Calculation/Easy Implementation for C++ Users</span>
</header>
<hr>
<div text-align: center;>
<table align="center"><tr>
<th>
<div style="text-align: center;">
    <span style="font-size: 35px; font-weight: 500;">1D Harmonic Oscillator</span>
    <br><br>
    <img src="Images/1d_ho_time.png">
    <br>
    <span style="font-weight: 700;"> Example 1: Simulation of 1D Harmonic Oscillator's time development </span>
</div>
</th></tr>
<tr><td>

```c++
#include "lib/Basis/Basis1D.h"
#include "lib/Operator/Hamiltonian1D.h"
#include "lib/Stream/QuFile.h"

double V(double x, const complex<double>* Phi, int i) {
    return 0.5*(x - 0.5)*(x - 0.5);
}

int main() {
    quon::Basis1D basis(100, 0.2);
    basis <<= quon::BOUNDARY_CONDITION(quon::BC_MODE::PERIODIC);
    basis <<= quon::WF_FORM::GAUSSIAN;
    quon::Hamiltonian1D H(1.0, 1.0, PotentialFunction);
    
    quon::QuFile file("data.txt");
    file <<= basis;
    for (int i = 0; i < 5000; i++) {
        basis = basis.evolve(0.01, H, false);
    }
    file <<= basis;
}
```
    
</td></tr></table>

<br>

<table align="center"><tr>
<th>
<div style="text-align: center;">
    <span style="font-size: 35px; font-weight: 500;">Superfluid</span>
    <br><br>
    <img src="Images/1d_gpe_static.png">
    <br>
    <span style="font-weight: 700;"> Example 2: Ground state of 1D Gross-Pitaevskii Equation </span>
</div>
</th></tr>
<tr><td>

```c++
#include "lib/Basis/Basis1D.h"
#include "lib/Operator/Hamiltonian1D.h"
#include "lib/Stream/QuFile.h"

double V(double x, const complex<double>* Phi, int i) {
    return 0.5*x*x + 10*norm(Phi[i]);
}

int main() {
    quon::Basis1D basis(100, 0.2);
    basis <<= quon::BOUNDARY_CONDITION(quon::BC_MODE::PERIODIC);
    basis <<= quon::WF_FORM::GAUSSIAN;
    quon::Hamiltonian1D H(1.0, 1.0, PotentialFunction);
    
    quon::QuFile file("data.txt");
    file <<= basis;
    for (int i = 0; i < 5000; i++) {
        basis = basis.evolve(0.01, H, true);
    }
    file <<= basis;
}
```
    
</td></tr></table>
</div>