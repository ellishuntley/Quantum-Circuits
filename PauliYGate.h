// PauliYGate.h

#ifndef PAULI_Y_GATE_H
#define PAULI_Y_GATE_H
#include "./QuantumGate.h"

// Y gate subclass
class PauliYGate : public QuantumGate 
{
public:
    explicit PauliYGate(int qubit) : QuantumGate(qubit + 1), qubit_(qubit) {}

    std::vector<std::vector<std::complex<double>>> get_matrix() const;

    void apply(QubitRegister& qreg) const;

private:
    int qubit_;
};

std::vector<std::vector<std::complex<double>>> PauliYGate::get_matrix() const 
{
    std::vector<std::vector<std::complex<double>>> matrix = {
        {0, -std::complex<double>(0, 1)},
        {std::complex<double>(0, 1), 0}
    };
    return matrix;
}

void PauliYGate::apply(QubitRegister& qreg) const 
{
    // Get the matrix for the Pauli Y gate
    std::vector<std::vector<std::complex<double>>> matrix = get_matrix();

    // Create a copy of the state vector
    std::vector<std::complex<double>> new_state = qreg.get_state();

    // Compute the new state vector by multiplying with the Pauli Y matrix
    unsigned int numStates = new_state.size();
    unsigned int mask = (1 << qubit_);
    for (unsigned int i = 0; i < numStates; i++) {
        if ((i & mask) == mask) {
            unsigned int i1 = i ^ mask;
            std::complex<double> v0 = new_state[i];
            std::complex<double> v1 = new_state[i1];
            new_state[i] = matrix[0][0] * v0 + matrix[0][1] * v1;
            new_state[i1] = matrix[1][0] * v0 + matrix[1][1] * v1;
        }
    }
    // Set the state vector to the new state
    qreg.set_state(new_state);
}

#endif