// SwapGate.h

#ifndef SWAP_GATE_H
#define SWAP_GATE_H
#include "./QuantumGate.h"

// Swap gate subclass
class SwapGate : public QuantumGate 
{
public:
    SwapGate(int qubit1, int qubit2) : QuantumGate(std::max(qubit1, qubit2) + 1), qubit1_(qubit1), qubit2_(qubit2) {}

    std::vector<std::vector<std::complex<double>>> get_matrix() const override;

    void apply(QubitRegister& qreg) const override;

private:
    int qubit1_;
    int qubit2_;
};

std::vector<std::vector<std::complex<double>>> SwapGate::get_matrix() const 
{ // matrix is included for completeness, but a more efficient method is used to implement SwapGate in SwapGate::apply
    return {
        {1, 0, 0, 0},
        {0, 0, 1, 0},
        {0, 1, 0, 0},
        {0, 0, 0, 1}
    };
}

void SwapGate::apply(QubitRegister& qreg) const 
{


    // Create a copy of the state vector
    std::vector<std::complex<double>> new_state = qreg.get_state();

    unsigned int numStates = new_state.size();
    unsigned int mask1 = (1 << qubit1_);
    unsigned int mask2 = (1 << qubit2_);

    for (unsigned int i = 0; i < numStates; i++) {
        if (((i & mask1) != 0) != ((i & mask2) != 0)) {
            unsigned int swappedIndex = (i ^ mask1) ^ mask2;
            if (i < swappedIndex) {
                std::swap(new_state[i], new_state[swappedIndex]);
            }
        }
    }

    // Set the state vector to the new state
    qreg.set_state(new_state);
}

#endif
