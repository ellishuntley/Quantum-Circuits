// main.cpp

// class includes
#include "./QubitRegister.h"
#include "./QuantumGate.h"
#include "./ControlledGate.h"
#include "./HadamardGate.h"
#include "./PauliXGate.h"
#include "./PauliYGate.h"
#include "./PauliZGate.h"
#include "./PhaseShiftGate.h"
#include "./CNOTGate.h"
#include "./SwapGate.h"
#include "./RkGate.h"

// other includes
#include <memory>

// Pre-built Circuits:

// QFT
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

// Inverse QFT
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

// Quantum teleportation
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

// Helper functions:

// Gets number for prompts and qubit index
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
    std::cout << "\nPress enter to continue...";
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::getchar();
}

// User interface functions:

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

// Main menu
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

int main()
{
    UI();
    return 0;
}