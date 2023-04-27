// RkGate.h

#ifndef RK_GATE_H
#define RK_GATE_H
#include "./QuantumGate.h"

// Rk gate subclass
class RkGate : public QuantumGate 
{
public:
    RkGate(int target_qubit, int k) : QuantumGate(target_qubit + 1), target_qubit_(target_qubit), k_(k) {}

    std::vector<std::vector<std::complex<double>>> get_matrix() const override;

    void apply(QubitRegister& qreg) const override;

private:
    int target_qubit_;
    int k_;
};

std::vector<std::vector<std::complex<double>>> RkGate::get_matrix() const
{
    double angle = 2 * M_PI / pow(2, k_);
    std::complex<double> phase(std::cos(angle), std::sin(angle));

    return {
        {1, 0},
        {0, phase}
    };
}

void RkGate::apply(QubitRegister& qreg) const
{
    // Get the state vector
    auto state = qreg.get_state();

    // Apply the R_k gate to the target qubit
    unsigned int target_mask = (1 << target_qubit_);
    for (unsigned int i = 0; i < state.size(); i++) {
        if (i & target_mask) {
            double angle = 2 * M_PI / pow(2, k_);
            std::complex<double> phase(std::cos(angle), std::sin(angle));
            state[i] *= phase;
        }
    }

    // Set the new state vector
    qreg.set_state(state);
}

#endif