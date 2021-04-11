
#include "mpi.h"
#include <iostream>
#include <chrono>
#include <random>

int main()
{
	MPI::Init();
	std::cout << "f" << std::endl;
	int rank = MPI::COMM_WORLD.Get_rank();
	int cluster = MPI::COMM_WORLD.Get_size();

	const unsigned int size = 1000000;
	const float vmax = 256;
}
