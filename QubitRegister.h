// QubitRegsiter.h

#ifndef QUBIT_REGISTER_H
#define QUBIT_REGISTER_H
#include <vector> 
#include<iostream>
#include<iomanip>
#include <complex>

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
        state = std::vector<std::complex<double> >(numStates);
        state[0] = 1.0;
        for (unsigned int i = 1; i < numStates; i++) {
            state[i] = 0.0;
        }
    }
    ~QubitRegister() {};
    // member function to return the state
    const std::vector<std::complex<double> >& get_state() const;

    // member function to a set a new state in the register
    void set_state(const std::vector<std::complex<double> >& new_state);

    // member function to print a definite state in accordance with the probabilities
    void measure(int qubit_index, std::vector<int>& results);

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
    std::vector<std::complex<double> > state;
};

#include <complex>
#include <random>


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

#endif

