#pragma once

#include <iostream>

class Solver
{
public:
	Solver();
	Solver(const char* init);
	Solver(const Solver* init);

	~Solver();

	void print(std::ostream &s);

	bool isSolved();
	bool isAllowed(char val, int x, int y);
	bool solveBackTrack();

	void set(char val, int x, int y);

private:
	char data[9][9];
};

