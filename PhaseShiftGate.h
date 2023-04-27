// PhaseShiftGate.h

#ifndef PHASE_SHIFT_GATE_H
#define PHASE_SHIFT_GATE_H
#include "./QuantumGate.h"

// Phase Shift gate subclass
class PhaseShiftGate : public QuantumGate 
{
public:
    PhaseShiftGate(int targetQubit, double angle, int controlQubit = -1)
        : QuantumGate(targetQubit + 1), targetQubit_(targetQubit), angle_(angle), controlQubit_(controlQubit) {}

    std::vector<std::vector<std::complex<double>>> get_matrix() const override;

    void apply(QubitRegister& qreg) const override;

private:
    int targetQubit_;
    double angle_;
    int controlQubit_;
};

std::vector<std::vector<std::complex<double>>> PhaseShiftGate::get_matrix() const 
{
    std::vector<std::vector<std::complex<double>>> matrix = {
        {1, 0},
        {0, std::polar(1.0, angle_)}
    };
    return matrix;
}

void PhaseShiftGate::apply(QubitRegister& qreg) const 
{
    // Get the matrix for the phase shift gate
    std::vector<std::vector<std::complex<double>>> matrix = get_matrix();

    // Create a copy of the state vector
    std::vector<std::complex<double>> new_state = qreg.get_state();

    // Compute the new state vector by multiplying with the phase shift matrix
    unsigned int numStates = new_state.size();
    unsigned int targetMask = (1 << targetQubit_);
    unsigned int controlMask = controlQubit_ >= 0 ? (1 << controlQubit_) : 0;
    for (unsigned int i = 0; i < numStates; i++) {
        if ((i & targetMask) == targetMask && (controlQubit_ < 0 || (i & controlMask) == controlMask)) {
            new_state[i] *= matrix[1][1];
        }
    }

    // Set the state vector to the new state
    qreg.set_state(new_state);
}

#endif