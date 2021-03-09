/*
 * simulator.cpp
 *
 *  Created on: 9 mars 2021
 *      Author: Joel Nilsson
 *
 *	Simulates a simple quantum circuit.
 *	Gates supported: X, Y, Z, Hadamard, sqrtNOT, CNOT.
 *	Supports computations with 4 qubits. May be scaled
 *	up by lifting restrictions in start().
 *
 */

#include<iostream>
#include<algorithm>
#include<sstream>
#include<vector>
#include"simulator.h"

typedef std::complex<double> cd;
typedef std::pair<std::complex<double>, std::complex<double>> cc;
typedef std::vector<std::complex<double>> vc;

/*
 * Starts simulations and handles input.
 */
void Simulator::start()
{
	printstart();
	printhelp();

	while(true)
	{

		std::string arg{}, line{};
		std::cin >> arg;

		std::transform(arg.begin(), arg.end(), arg.begin(), ::toupper);

		if(arg == "QUIT")
		{
			break;
		}
		else if(arg == "HELP")
		{
			printhelp();
		}
		else if(arg == "NEW")
		{
			std::cout << "Enter number of qubits: ";
			bool ok = false;
			do
			{
				if(std::cin >> nqubits && nqubits < 5 && nqubits > 0)
				{
					grid.resize(nqubits);
					ok = true;
				}
				else
				{
					std::cout << "Invalid number of qubits (enter a number between 1 and 4)." << std::endl;
					std::cout << "Enter number of qubits: ";
				}

			} while(!ok);

			std::cout << "Enter a grid of " << nqubits << " lines, each 15 characters wide:" << std::endl;

			for(int i{0}; i < nqubits; ++i)
			{
				bool reset = false;
				std::vector<char> gridline(15);
				std::cin >> line;
				if(line.length() < 15)
				{
					line += "---------------";
				}

				for(int j{0}; j < 15; ++j)
				{
					if(validGate(line[j]))
					{
						gridline[j] = line[j];
					}
					else
					{
						std::cout << "Invalid grid line, please enter line again:" << std::endl;
						reset = true;
					}
				}

				grid[i] = gridline;
				if(reset)
				{
					--i;
				}
			}

			// Print grid
			std::cout << "Quantum grid:" << std::endl;
			for(int i{0}; i < nqubits; ++i)
			{
				for(int j{0}; j < 15; ++j)
				{
					std::cout << grid[i][j];
				}
				std::cout << std::endl;
			}

			/*
			 * Setup computation
			 */
			do
			{
				bool invalid_init = false;
				std::cout << "\nEnter an initial state for (" << nqubits << ") qubit(s), or type BACK to go to menu: " << std::flush;
				std::cin >> arg;

				std::transform(arg.begin(), arg.end(), arg.begin(), ::toupper);
				if(arg == "BACK")
				{
					std::cout << "Type NEW to create a new circuit, or QUIT to exit the program. Type HELP for a tutorial." << std::endl;
					break;
				}
				else if(arg.length() < nqubits)
				{
					arg += "0000";
				}

				int statenr{0};
				int incr{1};

				for(int i{nqubits-1}; i >= 0; --i)
				{
					if(arg[i] == '1')
					{
						statenr += incr;
					}
					else if(arg[i] == '0')
					{
						// Add nothing
					}
					else
					{
						std::cout << "Invalid initial state, only 1s and 0s accepted." << std::endl;
						invalid_init = true;
						break;
					}

					incr <<= 1;
				}

				if(invalid_init)
				{
					continue;
				}

				int nstates = 1 << nqubits;
				state = vc(nstates, cd(0, 0));
				state[statenr] = cd(1, 0);

				/*
				 * Need to divide by norm to reduce rounding error?
				 * Divide by norm before presenting probability
				 */
				compute();

			} while (true);

		}

	}

}

/*
 * Performs computation on the given initial state
 * with the grid that the Simulator object has.
 */
void Simulator::compute()
{
	char op{};
	bool do_operation;

	for(int i{0}; i < 15; ++i)
	{
		int CTRLqbit{-1};
		int affectedqbit{-1};
		do_operation = false;
		for(int j{0}; j < nqubits; ++j)
		{
			if(grid[j][i] != '-')
			{
				do_operation = true;
				if(grid[j][i] == 'C')
				{
					affectedqbit = j;
				}
				else if(grid[j][i] == 'O')
				{
					CTRLqbit = j;
				}
			}
		}

		if(do_operation)
		{
			if(CTRLqbit != -1 && affectedqbit != -1)
			{
				CNOT(CTRLqbit, affectedqbit);
			}

			// Matrix for storing operations
			std::vector<std::vector<cd>> matrixop{};

			// Begin at last qubit to get correct multiplication.
			// Tensor product used for parallel operations.
			for(int j{nqubits-1}; j >= 0; --j)
			{
				op = grid[j][i];
				matrixop = tensorprod(op, matrixop);
			}

			int n = matrixop.size();
			vc newstate(n, cd(0, 0));

			// Perform matrix multiplication
			for(int i{0}; i < n; ++i)
			{
				for(int j{0}; j < n; ++j)
				{
					newstate[i] += state[j]*matrixop[j][i];
				}
			}

			state = newstate;
		}

	}

	/*
	 * Print state probabilities.
	 */
	int nstates = 1 << nqubits;

	std::cout << "State probabilities:";

	for(int i{0}; i < nstates; ++i)
	{
		double prob = norm(state[i]);
		if(prob < 1e-6)
		{
			prob = 0;
		}
		std::cout << "\n" << state_translator(i) << " " << prob;
	}
	std::cout << "\n" << std::endl;

	return;
}

/*
 * CNOT-operator. Toggles the affectedqbit of every state where
 * the CTRLqbit is on.
 */
void Simulator::CNOT(int CTRLqbit, int affectedqbit)
{
	vc newstate = state;
	int nstates = 1 << nqubits;

	for(int ix{0}; ix < nstates; ++ix)
	{
		if(ix & (1 << (nqubits - CTRLqbit -1)))
		{
			newstate[ix] = state[( ix ^ (1 << (nqubits - affectedqbit - 1)) )];
		}
	}

	state = newstate;

	return;
}

/*
 * Returns a complex matrix which is the tensor product of the entered
 * matrix 'mat', and the matrix operator represented by 'op'.
 */
std::vector<std::vector<cd>> Simulator::tensorprod(char op, std::vector<std::vector<cd>> mat)
{
	vc coef(4);

	cd i = cd(0, 1);
	cd cmp1 = cd(1, 0);
	cd cmp0 = cd(0, 0);
	// Operations CNOT handled in other function and called before.
	if(op == '-' || op == 'C' || op == 'O')
	{
		coef[0] = cmp1;
		coef[1] = cmp0;
		coef[2] = cmp0;
		coef[3] = cmp1;
	}
	else if(op == 'X')
	{
		coef[0] = cmp0;
		coef[1] = cmp1;
		coef[2] = cmp1;
		coef[3] = cmp0;
	}
	else if(op == 'Y')
	{
		coef[0] = cmp0;
		coef[1] = -i;
		coef[2] = i;
		coef[3] = cmp0;
	}
	else if(op == 'Z')
	{
		coef[0] = cmp1;
		coef[1] = cmp0;
		coef[2] = cmp0;
		coef[3] = -cmp1;
	}
	else if(op == 'N')
	{
		coef[0] = 0.5*cmp1+0.5*i;
		coef[1] = 0.5*cmp1-0.5*i;;
		coef[2] = 0.5*cmp1-0.5*i;;
		coef[3] = 0.5*cmp1+0.5*i;;
	}
	else if(op == 'H')
	{
		coef[0] = cd(1/sqrt(2),0);
		coef[1] = cd(1/sqrt(2),0);
		coef[2] = cd(1/sqrt(2),0);
		coef[3] = -cd(1/sqrt(2),0);
	}

	int sz = mat.size();

	if(sz == 0)
	{
		std::vector<std::vector<cd>> newmat1 = std::vector<std::vector<cd>>(2, std::vector<cd>(2));
		newmat1[0][0] = coef[0];
		newmat1[0][1] = coef[1];
		newmat1[1][0] = coef[2];
		newmat1[1][1] = coef[3];

		return newmat1;
	}

	std::vector<std::vector<cd>> newmat = std::vector<std::vector<cd>>(2*sz, std::vector<cd>(2*sz));

	for(int m{0}; m < 2; ++m)
	{
		for(int k{0}; k < 2; ++k)
		{
			for(int i{0}; i < sz; ++i)
			{
				for(int j{0}; j < sz; ++j)
				{
					newmat[i+k*sz][j+m*sz] = coef[k+2*m]*mat[i][j];
				}
			}
		}
	}

	return newmat;

}

/*
 * Returns true if the entered gate is valid.
 * Does not validate syntax of CNOT gates.
 */
bool Simulator::validGate(char gate)
{
	return (gate == '-' || gate == 'H' || gate == 'X' || gate == 'N' || gate == 'Y' || gate == 'Z' || gate == 'C' || gate == 'O' );
}

/*
 * Prints a tutorial.
 */
void Simulator::printhelp()
{
	std::cout << "User instructions:\n";
	std::cout << "Enter a number of qubits between 1 and 4, and press 'Enter'.\n";
	std::cout << "Enter a quantum grid with gates X, Y, Z, N (square root NOT), H, or CNOT (using C and O combined, see below:\n";
	std::cout << "-X-H-X---O--H--\n";
	std::cout << "-Y-O-H------N--\n";
	std::cout << "-X-C-----C--H--\n";
	std::cout << "Any C or O used by itself will be ignored. Only one O can be used in a given column if the grid.\n";
	std::cout << "As an example, if the qubit at O is in state 1 with probability 1, then the qubit at C will be inverted.\n";
	std::cout << "Finally, enter an initial state (in 0s or 1s) of the qubit(s) in the setup. Empty states will be filled in with zeros.\n";
	std::cout << "The first digit will be the state of the qubit in the topmost line.\n\n";
	std::cout << "After evaluation a new initial state can be provided as input, or the command BACK can be entered to setup a new quantum grid.\n";
	std::cout << "Type QUIT to exit the program.\n";
	std::cout << "Type HELP to view this message again\n";
	std::cout << "Type NEW to begin creating a new grid.\n" << std::endl;

	return;
}

/*
 * Prints a start-up message.
 */
void Simulator::printstart()
{
	std::cout << "......____                    __                     _____ _                 __      __            \n";
	std::cout << "...../ __ \\__  ______ _____  / /___  ______ ___     / ___/(_)___ ___  __  __/ /___ _/ /_____  _____\n";
	std::cout << "..../ / / / / / / __ `/ __ \\/ __/ / / / __ `__ \\    \\__ \\/ / __ `__ \\/ / / / / __ `/ __/ __ \\/ ___/\n";
	std::cout << ".../ /_/ / /_/ / /_/ / / / / /_/ /_/ / / / / / /   ___/ / / / / / / / /_/ / / /_/ / /_/ /_/ / /    \n";
	std::cout << "...\\___\\_\\__,_/\\__,_/_/ /_/\\__/\\__,_/_/ /_/ /_/   /____/_/_/ /_/ /_/\\__,_/_/\\__,_/\\__/\\____/_/     \n\n";
	std::cout << "Welcome to Quantum Simulator version 1.0\n\n";
	return;
}

/*
 * Returns a string with the representation of an integer i in binary.
 */
std::string Simulator::state_translator(int i)
{
	std::string state{};
	int rest{};
	while(i != 0)
	{
		if(i & 1)
		{
			state.append("1");
		}
		else
		{
			state.append("0");
		}
		i >>= 1;
	}

	while(state.length() < nqubits)
	{
		state.append("0");
	}

	std::reverse(state.begin(), state.end());

	return state;
}









