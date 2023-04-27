// HadamardGate.h

#ifndef HADAMARD_GATE_H
#define HADAMARD_GATE_H
#include "./QuantumGate.h"


// Hadamard gate subclass
class HadamardGate : public QuantumGate 
{
public:
    
    HadamardGate(int qubit) : QuantumGate(qubit + 1), qubit_(qubit) {}

    std::vector<std::vector<std::complex<double>>> get_matrix() const override;
    
    void apply(QubitRegister& qreg) const override;

private:
    int qubit_;
};

std::vector<std::vector<std::complex<double>>> HadamardGate::get_matrix() const
    {
        const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
        std::vector<std::vector<std::complex<double>>> matrix = {
            {inv_sqrt2, inv_sqrt2},
            {inv_sqrt2, -inv_sqrt2}
        };
        return matrix;
    }

void HadamardGate::apply(QubitRegister& qreg) const
    {
        // Get the matrix for the Hadamard gate
        std::vector<std::vector<std::complex<double>>> matrix = get_matrix();

        // Create a copy of the state vector
        std::vector<std::complex<double>> new_state = qreg.get_state();

        // Compute the new state vector by multiplying with the Hadamard matrix
        unsigned int numStates = new_state.size();
        unsigned int mask = (1 << qubit_);
        for (unsigned int i = 0; i < numStates; i++) {
            if ((i & mask) == 0) {
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