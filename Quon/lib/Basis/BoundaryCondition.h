#ifndef _QUON_BoundaryCondition_HEADER_
#define _QUON_BoundaryCondition_HEADER_

namespace quon {
    enum class BC_MODE {
        DIRICHLET,
        NEUMANN,
        PERIODIC,
        DEBUG
    };
    struct BC_VALUES {
        double alpha;
        double beta;

        BC_VALUES(double const _alpha = 0.0, double const _beta = 0.0) : alpha(_alpha), beta(_beta) {}
        ~BC_VALUES()
        {

        }

        bool operator == (BC_VALUES const values) const {
            return (this->alpha == values.alpha && this->beta == values.beta);
        }
    };

    class BOUNDARY_CONDITION {
    private:
        BC_MODE mode;
        BC_VALUES values;

    public:
        BOUNDARY_CONDITION(BC_MODE const _mode, BC_VALUES const _values) : mode(_mode), values(_values) { }
        BOUNDARY_CONDITION(BC_MODE const _mode, double const alpha = 0.0, double const beta = 0.0)
            : mode(_mode), values(alpha, beta) {}
        BOUNDARY_CONDITION() : mode(BC_MODE::DIRICHLET), values(0.0, 0.0) {};

        ~BOUNDARY_CONDITION() {
            
        }

        double inline yleft() const { return this->values.alpha; }

        double inline yright() const { return this->values.beta; }

        BC_MODE inline getMode() const { return this->mode; }

        bool operator == (BOUNDARY_CONDITION const condition) const {
            return (this->mode == condition.mode && this->values == condition.values);
        }
    };
}

#endif