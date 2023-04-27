// PauliZGate.h

#ifndef PAULI_Z_GATE_H
#define PAULI_Z_GATE_H
#include "./QuantumGate.h"

// Z gate subclass
class PauliZGate : public QuantumGate 
{
public:
    explicit PauliZGate(int qubit) : QuantumGate(qubit + 1), qubit_(qubit) {}

    std::vector<std::vector<std::complex<double>>> get_matrix() const;

    void apply(QubitRegister& qreg) const;

private:
    int qubit_;
};

std::vector<std::vector<std::complex<double>>> PauliZGate::get_matrix() const 
{
    std::vector<std::vector<std::complex<double>>> matrix = {
        {1, 0},
        {0, -1}
    };
    return matrix;
}

void PauliZGate::apply(QubitRegister& qreg) const 
{
    // Get the matrix for the Pauli Z gate
    std::vector<std::vector<std::complex<double>>> matrix = get_matrix();

    // Create a copy of the state vector
    std::vector<std::complex<double>> new_state = qreg.get_state();

    // Compute the new state vector by multiplying with the Pauli Z matrix
    unsigned int numStates = new_state.size();
    unsigned int mask = (1 << qubit_);
    for (unsigned int i = 0; i < numStates; i++) {
        if ((i & mask) == mask) {
            new_state[i] *= matrix[1][1];
        }
    }
    // Set the state vector to the new state
    qreg.set_state(new_state);
}

#endif
