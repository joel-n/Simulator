/*
 * simulator.h
 *
 *  Created on: 9 mars 2021
 *      Author: Joel Nilsson
 *
 *      Quantum circuit simulator class.
 *      See definitions in Simulator.cpp
 */

#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include <vector>
#include <complex>
#include <cmath>

typedef std::complex<double> cd;
typedef std::vector<cd> vc;

class Simulator
{
public:

	void start();

	void compute();

	void CNOT(int CTRLqbit, int affectedqbit);

	std::vector<std::vector<cd>> tensorprod(char op, std::vector<std::vector<cd>> mat);

	void printhelp();
	void printstart();

	bool validGate(char gate);

	std::string state_translator(int i);

private:

	const double PI = acos(-1);

	std::vector<std::vector<char>> grid;

	vc state;

	int nqubits;
};


#endif /* SIMULATOR_H_ */
