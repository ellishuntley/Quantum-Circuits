// ControleledGate.h

#ifndef CONTROLLED_GATE_H
#define CONTROLLED_GATE_H
#include "QubitRegister.h"
#include "QuantumGate.h"

// template class for making controlled single qubit gates without making new classes
template <class GateType>
class ControlledGate : public QuantumGate 
{
public:
    ControlledGate(int target_qubit, int control_qubit)
        : QuantumGate(std::max(target_qubit, control_qubit) + 1),
          target_qubit_(target_qubit),
          control_qubit_(control_qubit),
          base_gate_(target_qubit) {}

          virtual ~ControlledGate() {};

    std::vector<std::vector<std::complex<double>>> get_matrix() const override;
    void apply(QubitRegister& qreg) const override;

private:
    int target_qubit_;
    int control_qubit_;
    GateType base_gate_;
};

template <class GateType>
void ControlledGate<GateType>::apply(QubitRegister& qreg) const
{
    unsigned int numStates = qreg.get_state().size();
    unsigned int target_mask = (1 << target_qubit_);
    unsigned int control_mask = (1 << control_qubit_);

    std::vector<std::complex<double>> new_state = qreg.get_state();

    for (unsigned int i = 0; i < numStates; i++) {
        if ((i & control_mask) == control_mask) {
            // Apply the gate to the target qubit only when the control qubit is set
            unsigned int i1 = i ^ target_mask;
            std::complex<double> v0 = new_state[i & ~target_mask];
            std::complex<double> v1 = new_state[i1];

            std::vector<std::vector<std::complex<double>>> matrix = get_matrix();
            std::complex<double> new_v0 = matrix[0][0] * v0 + matrix[0][1] * v1;
            std::complex<double> new_v1 = matrix[1][0] * v0 + matrix[1][1] * v1;

            new_state[i & ~target_mask] = new_v0;
            new_state[i1] = new_v1;
        }
    }

    qreg.set_state(new_state);
}

template <class GateType>
std::vector<std::vector<std::complex<double>>> ControlledGate<GateType>::get_matrix() const
{
    return base_gate_.get_matrix();
}

#endif