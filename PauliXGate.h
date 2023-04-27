// PauliXGate.h

#ifndef PAULI_X_GATE_H
#define PAULI_X_GATE_H
#include "./QuantumGate.h"

// X gate subclass
class PauliXGate : public QuantumGate 
{
public:
    explicit PauliXGate(int target_qubit)
        : QuantumGate(target_qubit + 1),
          target_qubit_(target_qubit) {}

    std::vector<std::vector<std::complex<double>>> get_matrix() const override;

    void apply(QubitRegister& qreg) const override;

private:
    int target_qubit_;
};

std::vector<std::vector<std::complex<double>>> PauliXGate::get_matrix() const 
{
    std::vector<std::vector<std::complex<double>>> matrix = {
        {0, 1},
        {1, 0}
    };
    return matrix;
}

void PauliXGate::apply(QubitRegister& qreg) const 
{
    // Get the matrix for the Pauli X gate
    std::vector<std::vector<std::complex<double>>> matrix = get_matrix();

    // Create a copy of the state vector
    std::vector<std::complex<double>> new_state = qreg.get_state();

    // Compute the new state vector by multiplying with the Pauli X matrix
    unsigned int numStates = new_state.size();
    unsigned int target_mask = (1 << target_qubit_);
    for (unsigned int i = 0; i < numStates; i++) {
        if ((i & target_mask) == target_mask) {
            unsigned int i1 = i ^ target_mask;
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