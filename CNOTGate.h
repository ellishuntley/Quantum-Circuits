// CNOTGate.h

#ifndef CNOT_GATE_H
#define CNOT_GATE_H
#include "./QuantumGate.h"

// Controlled Not gate subclass
class CNOTGate : public QuantumGate 
{
public:
    CNOTGate(int target_qubit, int control_qubit) : QuantumGate(std::max(control_qubit, target_qubit) + 1), control_qubit_(control_qubit), target_qubit_(target_qubit) {}

    std::vector<std::vector<std::complex<double>>> get_matrix() const override;

    void apply(QubitRegister& qreg) const override;

private:
    int control_qubit_;
    int target_qubit_;
};

std::vector<std::vector<std::complex<double>>> CNOTGate::get_matrix() const
{ // matrix is included for completeness, but a more efficient method is used to implement CNOT in CNOTGate::apply
    return {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 1, 0}
    };
}

void CNOTGate::apply(QubitRegister& qreg) const
{
    // Create a copy of the state vector
    std::vector<std::complex<double>> oldState = qreg.get_state();
    std::vector<std::complex<double>> new_state(oldState.size());

    // Apply the CNOT gate to the state vector
    unsigned int numStates = oldState.size();
    unsigned int control_mask = (1 << control_qubit_);
    unsigned int target_mask = (1 << target_qubit_);
    for (unsigned int i = 0; i < numStates; i++) {
        if (i & control_mask) {
            unsigned int i1 = i ^ target_mask;
            new_state[i] = oldState[i1];
            new_state[i1] = oldState[i];
        } else {
            new_state[i] = oldState[i];
        }
    }
    // Set the state vector to the new state
    qreg.set_state(new_state);
}

#endif