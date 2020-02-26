#ifndef _QUON_QuFile_HEADER_
#define _QUON_QuFile_HEADER_

#include "../Basis/Basis1D.h"
#include <fstream>
#include <string>

namespace quon {
    class QuFile {
    protected:
        std::string filename;
        std::ofstream file;

    public:
        QuFile(std::string _filename) : filename(_filename) {
            file.open(filename);
            if (!this->file) {
                abort();
            }
        }

        ~QuFile() {
            this->close();
        }

        void close() {
            if (this->file) {
                this->file.close();
            }
        }

        void operator <<= (Basis1D basis1D) {
            std::complex<double>* Phi;
            Phi = basis1D.getPhi();
            for (int i = 0; i < basis1D.getSize(); i++) {
                file << basis1D.getXPos(i) << " " << norm(Phi[i]) << " " << real(Phi[i]) << " " << imag(Phi[i]) << std::endl;
            }
            file << std::endl;
            Phi = nullptr;
        }
    };
}

#endif