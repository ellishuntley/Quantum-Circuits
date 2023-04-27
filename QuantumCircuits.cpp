#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<string> 
#include<limits>
#include<sstream>
#include<vector>
#include<algorithm>
#include <math.h>
#include <cstdlib>
#include <complex>
#include <random>
#include <numeric>
#include <regex>
#include <istream>
#include <cstdio>

// Object Orientated Programming in C++
// Final Project
// Name: Ellis Robert Huntley
// Student ID: 10583262

// Qubit regsiter class
class QubitRegister 
{
    // Overloaded operator to print out the state of the qubit register
    friend std::ostream& operator<<(std::ostream& os, const QubitRegister& qreg);

    // Overloaded operator to input a user-defined state
    friend std::istream& operator>>(std::istream& is, QubitRegister& qreg);
public:
    QubitRegister(unsigned int numQubits) : numQubits(numQubits) 
    {
        unsigned int numStates = 1 << numQubits;
        state = std::vector<std::complex<double>>(numStates);
        state[0] = 1.0;
        for (unsigned int i = 1; i < numStates; i++) {
            state[i] = 0.0;
        }
    }
    ~QubitRegister() {};
    // member function to return the state
    const std::vector<std::complex<double>>& get_state() const;
    // member function to a set a new state in the register
    void set_state(const std::vector<std::complex<double>>& new_state);
    // member function to print a definite state in accordance with the probabilities
    void measure(int qubit_index, std::vector<int>& results); // measures qubit
    // returns number of qubits
    int get_num_qubits() const;
    // print function
    void print() const;
    // reset register to |0*n>
    void reset();
    // sets random state
    void generate_random_state(int target_qubit_index);
private:
    unsigned int numQubits;
    std::vector<std::complex<double>> state;
};

// Generates a random, normalised state
void QubitRegister::generate_random_state(int target_qubit_index)
{
    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);

    // Generate random complex numbers for a single qubit state
    std::complex<double> random_state[2];
    double sum = 0.0;
    for (int i = 0; i < 2; i++) {
        random_state[i] = std::complex<double>(dis(gen), dis(gen));
        sum += std::norm(random_state[i]);
    }

    // Normalize the amplitudes
    double scale_factor = 1.0 / std::sqrt(sum);
    for (int i = 0; i < 2; i++) {
        random_state[i] *= scale_factor;
    }

    // Create a new state vector to store the updated state
    std::vector<std::complex<double>> new_state(state.size(), std::complex<double>(0.0, 0.0));

    // Update the state of the target qubit in the new state vector
    for (int i = 0; i < state.size(); i++) {
        int other_qubits = i & ~(1 << target_qubit_index);
        if (((i >> target_qubit_index) & 1) == 1) {
            new_state[i] = random_state[1] * state[other_qubits];
        } else {
            new_state[i] = random_state[0] * state[other_qubits];
        }
    }

    // Update the state vector with the new state
    state = new_state;
}

// Reverts to |0>
void QubitRegister::reset() 
{
    for (int i = 0; i < numQubits; i++) {
        state[i] = std::complex<double>(0, 0);
    }
    state[0] = std::complex<double>(1, 0);
}

// Prints state 
void QubitRegister::print() const
{
    std::cout << "State: " << *this << "\n"<< std::endl;
}

// Returns number of qubits
int QubitRegister::get_num_qubits() const 
{
    return numQubits;
}

// Collapses state of given qubit
void QubitRegister::measure(int qubit_index, std::vector<int>& results) 
{
    std::vector<double> probabilities;
    double norm = 0.0;
    for (const auto& amplitude : state) {
        double probability = std::norm(amplitude);
        probabilities.push_back(probability);
        norm += probability;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> dist(probabilities.begin(), probabilities.end());
    int index = dist(gen);

    if (qubit_index >= 0 && qubit_index < numQubits) {
        // Collapse the state to the measurement outcome on the specific qubit
        double p0 = 0.0, p1 = 0.0;
        for (unsigned int i = 0; i < state.size(); i++) {
            if ((i & (1 << qubit_index)) == 0) {
                p0 += std::norm(state[i]);
            } else {
                p1 += std::norm(state[i]);
            }
        }
        double normalizationFactor = 1.0 / std::sqrt((index & (1 << qubit_index)) ? p1 : p0);
        for (unsigned int i = 0; i < state.size(); i++) {
            if ((i & (1 << qubit_index)) != (index & (1 << qubit_index))) {
                state[i] = 0.0;
            } else {
                state[i] *= normalizationFactor;
            }
        }
        
        // Save the measurement result
        results.push_back((index & (1 << qubit_index)) != 0);
    } else {
        // Collapse the state to the measurement outcome on the entire register
        std::vector<std::complex<double>> new_state(state.size(), 0);
        new_state[index] = 1.0;
        state = new_state;

        // Save the measurement result
        results.push_back(index);
    }
}

// Sets new state
void QubitRegister::set_state(const std::vector<std::complex<double>>& new_state)
{
    if (new_state.size() != state.size()) {
        throw std::invalid_argument("New state has wrong size");
    }
    state = new_state;
}

// Returns current state
const std::vector<std::complex<double>>& QubitRegister::get_state() const 
{
    return state;
}

// Overloads ostream to print in ket form
std::ostream& operator<<(std::ostream& os, const QubitRegister& qreg)
{
    unsigned int numStates = qreg.state.size();
    unsigned int numQubits = qreg.numQubits;

    // Flag to indicate if the current state is the first non-zero ket in the output
    bool firstKet = true;

    // Iterate through all possible states
    for (unsigned int i = 0; i < numStates; i++) {
        // Extract real and imaginary parts of the complex coefficient
        double re = std::real(qreg.state[i]);
        double im = std::imag(qreg.state[i]);

        // Check if either the real or imaginary part is non-zero
        if (std::abs(re) > 1e-10 || std::abs(im) > 1e-10) {
            // If it's not the first non-zero ket, add a plus or minus sign
            if (!firstKet) {
                if (re < 0) {
                    os << " - ";
                    re = -re;
                } else if (im < 0) {
                    os << " - ";
                    im = -im;
                } else {
                    os << " + ";
                }
            } else {
                firstKet = false;
                if (re < 0) {
                    os << "-";
                    re = -re;
                }
            }

            // Determine if the coefficient is complex (both real and imaginary parts are non-zero)
            bool isComplex = std::abs(re) > 1e-10 && std::abs(im) > 1e-10;

            // If the coefficient is complex, add an opening parenthesis
            if (isComplex) os << "(";

            // If the real part is non-zero, add it to the output
            if (std::abs(re) > 1e-10) {
                os << std::fixed << std::setprecision(3) << re;
            }

            // If the imaginary part is non-zero, add it to the output with the correct sign
            if (std::abs(im) > 1e-10) {
                if (std::abs(re) > 1e-10) {
                    os << (im < 0 ? " - " : " + ");
                }
                os << std::fixed << std::setprecision(3) << std::abs(im) << "i";
            }

            // If the coefficient is complex, add a closing parenthesis
            if (isComplex) os << ")";

            // Add the ket notation for the current state
            os << "|";
            for (int j = numQubits - 1; j >= 0; j--) {
                if ((i >> j) & 1) {
                    os << "1";
                } else {
                    os << "0";
                }
            }
            os << ">";
        }
    }
    return os;
}

// Overload the insertion operator for user-defined states
std::istream& operator>>(std::istream& is, QubitRegister& qreg)
{
    unsigned int numStates = 1 << qreg.get_num_qubits();
    std::vector<std::complex<double>> new_state(numStates);
    double a, b, c, d;
    char ch1, ch2, ch3, ch4, ch5, ch6;
    std::string ket0, ket1;
    
    is >> ch1 >> a >> ch2 >> b >> ch3 >> ket0 >> ch4 >> ch5 >> c >> ch6 >> d >> ket1;
    new_state[0] = std::complex<double>(a, b);
    new_state[1] = std::complex<double>(c, d);

    // Check for normalization
    double total_probability = 0.0;
    for (const auto& amplitude : new_state) {
        total_probability += std::norm(amplitude);
    }

    if (std::abs(total_probability - 1.0) > 1e-3) {
        throw std::runtime_error("The state is not properly normalised.");
    }

    qreg.set_state(new_state);
    return is;
}

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

// QFT pre-built
void QFT(QubitRegister& qreg, int numQubits) 
{
    std::vector<std::shared_ptr<QuantumGate>> gates;

    for (int i = 0; i < numQubits; i++) {
        // Apply Hadamard gate to qubit i
        gates.push_back(std::make_shared<HadamardGate>(i));

        // Apply controlled R_k gates
        for (int j = 1; j <= numQubits - i - 1; j++) {
            gates.push_back(std::make_shared<RkGate>(i + j, j));
        }
    }

    // Swap the order of qubits
    for (int i = 0; i < numQubits / 2; i++) {
        gates.push_back(std::make_shared<SwapGate>(i, numQubits - i - 1));
    }

    for (const auto& gate : gates) {
        gate->apply(qreg);
    }
}

// Inverse QFT pre-built
void invQFT(QubitRegister& qreg, int numQubits) 
{
    std::vector<std::shared_ptr<QuantumGate>> gates;

    // Swap the order of qubits
    for (int i = 0; i < numQubits / 2; i++) {
        gates.push_back(std::make_shared<SwapGate>(i, numQubits - i - 1));
    }

    for (int i = numQubits - 1; i >= 0; i--) {
        // Apply controlled inverse R_k gates
        for (int j = numQubits - i - 1; j >= 1; j--) {
            gates.push_back(std::make_shared<RkGate>(i + j, -j));
        }

        // Apply Hadamard gate to qubit i
        gates.push_back(std::make_shared<HadamardGate>(i));
    }

    for (const auto& gate : gates) {
        gate->apply(qreg);
    }
}

// Quantum teleportation pre-built
void quantum_teleportation(QubitRegister& qreg, int helper_index, int alice_qubit_index, int bob_qubit_index) 
{
    std::vector<std::shared_ptr<QuantumGate>> gates;

    gates.push_back(std::make_shared<HadamardGate>(alice_qubit_index));
    gates.push_back(std::make_shared<ControlledGate<PauliXGate>>(bob_qubit_index, alice_qubit_index));
    gates.push_back(std::make_shared<ControlledGate<PauliXGate>>(alice_qubit_index, helper_index));
    gates.push_back(std::make_shared<HadamardGate>(helper_index));


    for (const auto& gate : gates) {
        gate->apply(qreg);
    }


    // STEP 3: Measure Alice's qubits
    std::vector<int> measurement_results;
    qreg.measure(helper_index, measurement_results);
    qreg.measure(alice_qubit_index, measurement_results);
    int alice_measurement = measurement_results[0];
    int entangled_measurement = measurement_results[1];

    // STEP 4: Perform the necessary operations on Bob's side based on the measurement results
    if (entangled_measurement) {
        PauliXGate pxg3(bob_qubit_index);
        pxg3.apply(qreg);
    }
    if (alice_measurement) {
        PauliZGate pzg(bob_qubit_index);
        pzg.apply(qreg);
    }
}

// Function to get number for prompts and qubit index
int get_number(int options) 
{
    int n;
    bool integer; 
    std::string user_input;
    
    do { // do while loop so that user input is non-zero
        do { // do while loop that loops until the input is of the correct format
            integer = true;
            std::cout << "Enter number: " << std::endl;
            std::getline(std::cin, user_input);

            if (user_input == "") {
                integer = false;
                std::cout << "Invalid input" << std::endl; // hit enter to display the prompt again when running program.
                std::cin.clear();
                std::cin.ignore(1000, '\n');
                continue;
            }

            for (int i=0; i<user_input.size(); i++) { // loops through every character in the string
                if (std::isspace(user_input[i]) || !std::isdigit(user_input[i])) { // checks the current character is an integer
                    integer = false;
                    std::cout << "Invalid input" << std::endl; // hit enter to display the prompt again when running program.
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // ignores the largest possible input
                    break; 
                }
            }
        } while(!integer);

        n = std::stoi(user_input);

        if (n < 0 || n > options - 1) { // only allows the user to enter one of the correct numbers
            std::cout << "Invalid input" << std::endl;
            std::cin.clear();
            std::cin.ignore(1000, '\n');
        }

    } while (n < 0 || n > options -1);

    return n;
}

// Prompts the user to continue once finished reading previous output
void press_any_key_to_continue() 
{
    std::cout << "\nPress any key to continue...";
    std::cin.clear(); // clear any error flags
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // discard any remaining characters
    getchar(); // wait for user to press any key
}

// Function prototypes - necessary so functions can be used in any order
void UI(); 
void info_menu();
void build_circuit();
void add_component(std::vector<QuantumGate*>& gates, QubitRegister& qreg);
void pre_built();

// Contains information about the gates
void info_menu()
{
    std::cout << "\n0) Hadamard Gate \n1) Pauli X Gate \n2) Pauli Y Gate \n3) Pauli Z Gate \n4) Phase Shift Gate \n5) Controlled-Not Gate \n6) Swap Gate \n7) Rk Gate \n8) Main Menu" <<std::endl;
    int instruction = get_number(9);

    switch (instruction){
        case 0:
            std::cout << "\nMatrix: \n1 1 \n1 -1 \n" << std::endl;
            std::cout<< "The Hadamard gate transforms a qubit from the computational basis to the Hadamard basis, which is a superposition of the computational basis states. \n" << std::endl;
            std::cout << "H|0> = |0> + |1>  \nH|1> = |0> - |1>" << std::endl;
            press_any_key_to_continue();
            info_menu();
            break;
        case 1:
            std::cout << "\nMatrix: \n0 1 \n1 0 \n" << std::endl;
            std::cout << "The Pauli X gate flips the state of a qubit. \n" << std::endl;
            std::cout << "X|0> = |1>  \nX|1> = |0>" << std::endl;
            press_any_key_to_continue() ;
            info_menu();
            break;
        case 2:
            std::cout << "\nMatrix: \n0 -i \ni 0 \n" << std::endl;
            std::cout << "The Pauli Y gate flips the state of a qubit and introduces a phase shift of pi/2. \n" << std::endl;
            std::cout << "Y|0> = i|1>  \nY|1> = -i|0>" << std::endl;
            press_any_key_to_continue();
            info_menu();
            break;
        case 3:
            std::cout << "\nMatrix: \n1 0 \n0 -1 \n" << std::endl;
            std::cout << "The Pauli Z gate introduces a phase shift of pi. \n" << std::endl;
            std::cout << "Z|0> = |0>  \nZ|1> = -|1>" << std::endl;
            press_any_key_to_continue();
            info_menu();
            break;
        case 4:
            std::cout << "\nMatrix: \n1 0 \n0 e^(i*pi*x) \n" << std::endl;
            std::cout << "The Phase Shift gate introduces a phase shift of x*pi to a qubit. \n" << std::endl;
            press_any_key_to_continue();
            info_menu();
            break;
        case 5:
            std::cout << "\nMatrix: \n1 0 0 0 \n0 1 0 0 \n0 0 0 1 \n0 0 1 0 \n" << std::endl;
            std::cout << "The Controlled-Not gate flips the target qubit if the control qubit is in the state |1>. Otherwise, it leaves the target qubit unchanged. \n" << std::endl;
            std::cout << "If the control qubit is in state |0>: \nCNOT|00> = |00> \nCNOT|01> = |01> \nCNOT|10> = |11> \nCNOT|11> = |10>" << std::endl;
            std::cout << "If the control qubit is in state |1>: \nCNOT|00> = |01> \nCNOT|01> = |00> \nCNOT|10> = |10> \nCNOT|11> = |11>" << std::endl;
            press_any_key_to_continue();
            info_menu();
            break;
        case 6:
            std::cout << "Matrix: \n1 0 0 0 \n0 0 1 0 \n0 1 0 0 \n0 0 0 1 \n" << std::endl;
            std::cout << "The Swap gate swaps the states of two qubits. \n" << std::endl;
            std::cout << "Swap|00> = |00> \nSwap|01> = |10> \nSwap|10> = |01> \nSwap|11> = |11>" << std::endl;
            press_any_key_to_continue();
            info_menu();
            break;
        case 7:
            std::cout << "\nMatrix: \n1 0 \n0 e^(i*k*pi / 2) \n" << std::endl;
            std::cout << "The Rk gate is a family of gates that introduce a phase shift of k*pi to a single qubit. Similar to the Phase Shift gate, this gate is crucial to implementing the Quantum Fourier Transform, and its inverse." << std::endl;
            press_any_key_to_continue();
            info_menu();
            break;
        case 8:
            UI();
            break;
    }
}

// Contains all gates for the user to implement
void add_component(std::vector<std::shared_ptr<QuantumGate>>& gates, QubitRegister& qreg)
{
    std::cout << "Which component?" << std::endl;
    std::cout << "\n0) Hadamard Gate \n1) Pauli X Gate \n2) Pauli Y Gate \n3) Pauli Z Gate \n4) Phase Shift Gate \n5) Controlled-Not Gate \n6) Swap Gate \n7) Rk Gate" <<std::endl;
    int instruction = get_number(8);
    int qubit;
    int swap;
    int control;
    double phase;
    int k;
    int type;

    switch (instruction){
        case 0:
            std::cout << "\n0) Not controlled \n1) Controlled" << std::endl;
            type = get_number(2);
            switch (type) {
                case 0:
                    std::cout << "Apply to which (0-indexed) qubit?" << std::endl;
                    qubit = get_number(qreg.get_num_qubits());
                    gates.push_back(std::make_shared<HadamardGate>(qubit));
                    break;
                case 1:
                    std::cout << "Apply to which (0-indexed) qubit?" << std::endl;
                    qubit = get_number(qreg.get_num_qubits());
                    std::cout << "Control qubit?" << std::endl;
                    control = get_number(qreg.get_num_qubits());
                    gates.push_back(std::make_shared<ControlledGate<HadamardGate>>(qubit, control));
                    break;
            }
            break;

        case 1:
            std::cout << "\n0) Not controlled \n1) Controlled" << std::endl;
            type = get_number(2);
            switch (type) {
                case 0:
                    std::cout << "Apply to which (0-indexed) qubit?" << std::endl;
                    qubit = get_number(qreg.get_num_qubits());
                    gates.push_back(std::make_shared<PauliXGate>(qubit));
                    break;
                case 1:
                    std::cout << "Apply to which (0-indexed) qubit?" << std::endl;
                    qubit = get_number(qreg.get_num_qubits());
                    std::cout << "Control qubit?" << std::endl;
                    control = get_number(qreg.get_num_qubits());
                    gates.push_back(std::make_shared<ControlledGate<PauliXGate>>(qubit, control));
                    break;
            }
            break;
        case 2:
            std::cout << "\n0) Not controlled \n1) Controlled" << std::endl;
            type = get_number(2);
            switch (type) {
                case 0:
                    std::cout << "Apply to which (0-indexed) qubit?" << std::endl;
                    qubit = get_number(qreg.get_num_qubits());
                    gates.push_back(std::make_shared<PauliYGate>(qubit));
                    break;
                case 1:
                    std::cout << "Apply to which (0-indexed) qubit?" << std::endl;
                    qubit = get_number(qreg.get_num_qubits());
                    std::cout << "Control qubit?" << std::endl;
                    control = get_number(qreg.get_num_qubits());
                    gates.push_back(std::make_shared<ControlledGate<PauliYGate>>(qubit, control));
                    break;
            }
            break;
        case 3:
            std::cout << "\n0) Not controlled \n1) Controlled" << std::endl;
            type = get_number(2);
            switch (type) {
                case 0:
                    std::cout << "Apply to which (0-indexed) qubit?" << std::endl;
                    qubit = get_number(qreg.get_num_qubits());
                    gates.push_back(std::make_shared<PauliZGate>(qubit));
                    break;
                case 1:
                    std::cout << "Apply to which (0-indexed) qubit?" << std::endl;
                    qubit = get_number(qreg.get_num_qubits());
                    std::cout << "Control qubit?" << std::endl;
                    control = get_number(qreg.get_num_qubits());
                    gates.push_back(std::make_shared<ControlledGate<PauliZGate>>(qubit, control));
                    break;
            }
            break;

        case 4:
            std::cout << "Apply to which (0-indexed) qubit?" << std::endl;
            qubit = get_number(qreg.get_num_qubits());
            
            std::cout << "Enter phase as a number in radians: " << std::endl;
            while (true) {
                try {
                    std::cin >> phase;
                    if (std::cin.fail()) {
                        throw std::runtime_error("Invalid input");
                    } else {
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        break;
                    }
                }
                catch (const std::runtime_error& e) {
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::cout << "Invalid input. Please enter a number for the phase: " << std::endl;
                }
            }
            gates.push_back(std::make_shared<PhaseShiftGate>(qubit, phase));
            break;
        case 5:
            std::cout << "Apply to which (0-indexed) qubit?" << std::endl;
            qubit = get_number(qreg.get_num_qubits());
            std::cout << "Control qubit?" << std::endl;
            control = get_number(qreg.get_num_qubits());
            gates.push_back(std::make_shared<CNOTGate>(qubit, control));
            break;
        case 6:
            std::cout << "First qubit." << std::endl;
            qubit = get_number(qreg.get_num_qubits());
            std::cout << "Second qubit." << std::endl;
            swap = get_number(qreg.get_num_qubits());
            gates.push_back(std::make_shared<SwapGate>(qubit, swap));
            break;
        case 7:
            std::cout << "Apply to which (0-indexed) qubit?" << std::endl;
            qubit = get_number(qreg.get_num_qubits());
            std::cout << "Enter k-value (integer): " << std::endl;
            std::string user_input;
            while (true) {
                try {
                    std::getline(std::cin, user_input);

                    if (user_input.empty()) {
                        throw std::runtime_error("Invalid input");
                    } 
                    for (char c : user_input) {
                        if (!std::isdigit(c)) {
                            throw std::runtime_error("Invalid input");
                        }
                    }
                    k = std::stoi(user_input);
                    break;
                }
                catch (const std::runtime_error& e) {
                    std::cout << "Invalid input. Please enter an integer for k: " << std::endl;
                }
            }
            gates.push_back(std::make_shared<RkGate>(qubit, k));
            break;
    }
}

// Lets user build and interact with circuits
void build_circuit()
{
    std::cout << "\nNumber of qubits (generated in |0>): \n  " << std::endl;
    int num_qubits;
    do {
        num_qubits = get_number(std::numeric_limits<int>::max());
        if (num_qubits == 0) {
            std::cout << "Number of qubits must be greater than 0." << std::endl;
        }
    } while (num_qubits == 0);
    std::vector<std::shared_ptr<QuantumGate>> gates;
    QubitRegister qreg(num_qubits);
    
    bool continue_loop = true;
    while (continue_loop){
        std::cout << "\n0) Add gate \n1) Apply gates(s) \n2) Measure qubit \n3) Print state \n4) Clear circuit \n5) Main menu \n" << std::endl;
        int instruction = get_number(6);
        switch(instruction){
            case 0:
                add_component(gates, qreg);
                break;
            case 1:
                qreg.reset(); // Reset the state to the initial so user gates are not reapplied if used more than once
                for (const auto& gate : gates) {
                    gate->apply(qreg);
                }
                std::cout << "\nApplied!" << std::endl;
                press_any_key_to_continue();
                break;
            case 2: {
                std::cout << "Qubit to be measured (0-indexed): " << std::endl;
                int measure = get_number(num_qubits);
                std::vector<int> measurement_results;
                qreg.measure(measure, measurement_results);
                std::cout << "\nQubit " << measure <<  " measured!" << std::endl;
                press_any_key_to_continue();
                break;
            }

            case 3:
                std::cout << "\n" << std::endl;
                qreg.print();
                press_any_key_to_continue();
                break;
            case 4:
                 // clear circuit so user can make a new one
                gates.clear();
                qreg.reset();
                std::cout << "\nCircuit cleared!" << std::endl;
                press_any_key_to_continue();
                build_circuit();
                break;
            case 5:
                continue_loop = false;
                UI();
                break;
        }
    }
}

// Lets user use pre-built circits
void pre_built()
{
    std::cout << "\n0) Quantum Fourier Transform \n1) Inverse Quantum Fourier Transform \n2) Quantum Teleportation" << std::endl;
    int instruction = get_number(3);
    switch (instruction){
        case 0: {
            std::cout << "\n\nThe Quantum Fourier Transform (QFT) is a unitary transformation that maps a quantum state from the computational basis to the Fourier basis.\nThe QFT is a key component in many quantum algorithms, including Shor's algorithm for factoring large numbers and the quantum phase estimation algorithm.\n\nTo perform a QFT on a quantum register, first apply a Hadamard gate to each qubit, which creates a superposition of all possible states. Then, apply a series of controlled R_k gates to introduce relative phases between the basis states. Finally, swap the order of the qubits to obtain the transformed state." << std::endl;
            QubitRegister qreg(2);
            qreg.print();
            QFT(qreg, 2);
            std::cout << "\nApplying QFT circuit...\n" << std::endl;
            std::cout << "Transformed ";
            qreg.print();
            press_any_key_to_continue();
            UI();
            break;
        }
        case 1: {
            std::cout << "\n\nThe Inverse Quantum Fourier Transform (IQFT) is the inverse of the Quantum Fourier Transform. It maps a quantum state from the Fourier basis to the computational basis.\nLike the QFT, the IQFT is necessary to implement many quantum algorithms, including the quantum phase estimation algorithm.\n\nTo perform an inverse QFT on a quantum register, first reverse the order of the qubits. Then, apply a series of controlled R_k gates with negative angles to remove the relative phases between the basis states. Finally, apply a Hadamard gate to each qubit to obtain the transformed state.\n" << std::endl;
            QubitRegister qreg(2);
            QFT(qreg, 2);
            std::cout << "QFT ";
            qreg.print();
            invQFT(qreg, 2);
            std::cout << "Applying inverse QFT circuit...\n" << std::endl;
            std::cout << "Transformed ";
            qreg.print();
            press_any_key_to_continue();
            UI();
            break;
        }
        case 2: {
            std::cout << "\n\nIn quantum teleportation, an unknown quantum state is transmitted from one party (Alice) to another (Bob) using a previously shared entangled state and classical communication. In this example, we'll use a 3-qubit system, where qubit 0 represents Alice's qubit, qubit 1 represents the entangled qubit, and qubit 2 represents Bob's qubit. The process involves creating an entangled state between Alice's qubit and the entangled qubit, applying a series of gates to Alice's qubit and the entangled qubit, measuring Alice's qubit and the entangled qubit, and applying correction gates to Bob's qubit based on the measurement results.\n" << std::endl;
            std::string teleportation_result_text = "After performing quantum teleportation, qubit 2 has successfully taken on the state of the original state, up to a phase, while the other qubits have been measured, collapsing their states. As a result of the protocol, the initial state has been transferred from qubit 0 to qubit 2, effectively 'teleporting' the quantum information while preserving the superposition and coherence. The original states of qubits 0 and 1 are no longer accessible due to the measurements made during the teleportation process.";
            QubitRegister qreg(3);

            int alice_qubit_index = 0;
            int entangled_qubit_index = 1;
            int bob_qubit_index = 2;

            std::cout << "0) Enter your own state \n1) Random state\n" << std::endl;
            int option = get_number(2);

            switch (option){
                case 0: {
                    // Let user enter their own state for a Alice's qubit, keep qubit 1 and 2 in |0>.
                    bool state_valid = false;
                    while (!state_valid) {
                        try {
                            std::cout << "\nEnter the state for Alice's qubit in the format (a+bi)|0> + (c+di)|1>:" << std::endl;
                            std::cin >> qreg;
                            state_valid = true;
                        } catch (std::runtime_error& e) {
                            std::cout << "Error: " << e.what() << std::endl;
                            std::cout << "Please re-enter the state." << std::endl;
                            std::cin.clear(); // Clear error state
                            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore the rest of the line
                            }
                        }
                    std::cout << "\nInputted ";
                    qreg.print();
                    std::cout << "Building circuit...\nTeleporting state..." << std::endl;
                    quantum_teleportation(qreg, alice_qubit_index, entangled_qubit_index, bob_qubit_index);
                    std::cout << "\nAfter quantum teleportation:" << std::endl;
                    qreg.print();
                    std::cout << teleportation_result_text << std::endl;
                    press_any_key_to_continue();
                    UI();
                    break;
                }

                case 1: {
                    // Prepare the initial state of Alice's qubit with a random state
                    qreg.generate_random_state(alice_qubit_index);
                    std::cout << "\nRandom ";
                    qreg.print();
                    std::cout << "Building circuit...\nTeleporting state..." << std::endl;
                    quantum_teleportation(qreg, alice_qubit_index, entangled_qubit_index, bob_qubit_index);
                    std::cout << "\nAfter quantum teleportation:" << std::endl;
                    qreg.print();
                    std::cout << teleportation_result_text << std::endl;
                    press_any_key_to_continue();
                    UI();
                    break;
                }
            }
        break;
        }
    }
}

// Main user interface
void UI()
{
    std::cout << "\nChoose an option: \n \n0) Gate Information \n1) Build your own circuit \n2) Pre-built circuits \n3) Exit program" << std::endl;
    int instruction = get_number(4); // prompts the user to enter a number within allowed range
    switch(instruction){
        case 0:
            // information menu
            info_menu();
            break;
            
        case 1:
            //circuit builder 
            build_circuit();
            break;
        case 2:
            pre_built();
            break;
        case 3:
            std::cout << "\nGoodbye!" << std::endl; 
            break;
    }
}

// Main function
int main() 
{
    UI();
    return 0;
}