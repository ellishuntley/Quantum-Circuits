// QuantumGate.h

#ifndef QUANTUM_GATE_H
#define QUANTUM_GATE_H
#include "./QubitRegister.h"
#include <vector>
#include <complex>

// Quantum gate base class
class QuantumGate 
{
public:
    QuantumGate(int num_qubits) : num_qubits_(num_qubits) {}

    virtual ~QuantumGate() {}

    virtual std::vector<std::vector<std::complex<double>>> get_matrix() const = 0;

    virtual void apply(QubitRegister& qreg) const = 0;

private:
    int num_qubits_;
};

#endif