#include "mpi.h"
#include <iostream>
#include <chrono>
#include <random>

void add(int* invec, int* inoutvec, int* len, MPI_Datatype* dtype)
{
	for (int i = 0; i < *len; i++) {
		inoutvec[i] += invec[i];
	}
}


int main()
{
	MPI::Init();
	std::cout << "f" << std::endl;
	int rank = MPI::COMM_WORLD.Get_rank();
	int cluster = MPI::COMM_WORLD.Get_size();
	const unsigned int size = 1000000;
	const int vmax = 100;
	const int barrier = 50;

	int extra = size % cluster;
	int chunk = size / cluster;
	if (extra > 0) {
		chunk++;
	}

	std::vector<int> vector;
	if (rank == 0) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<int> distr(0, vmax);
		for (size_t index = 0; index < size; ++index)
		{
			vector.push_back(distr(gen));
		}
		for (size_t index = 0; index < extra; ++index)
		{
			vector.push_back(vmax);
		}
	}
	auto start = std::chrono::high_resolution_clock::now();
	std::cout << "g" << std::endl;
	std::vector<int> work(chunk);
	MPI::COMM_WORLD.Scatter(vector.data(), chunk, MPI_INT32_T, work.data(), chunk, MPI_INT32_T, 0);
	std::vector<int> out(chunk);
	for (int i = 0; i < chunk; ++i) {
		if (work[i] < barrier) {
			out[i] = 1;
		} else {
			out[i] = 0;
		}
	}
	std::cout << "h" << std::endl;
	std::vector<int> middle(size + extra);
	MPI::COMM_WORLD.Gather(out.data(), chunk, MPI_INT32_T, middle.data(), chunk, MPI_INT32_T, 0);
	if (rank == 0) {
		std::cout << "i" << std::endl;
		std::vector<int> last(size + extra);
		MPI::Op op;
		op.Init((MPI::User_function*)add, true);
		MPI::COMM_WORLD.Exscan(middle.data(), last.data(), size + extra, MPI_INT32_T, op);
		op.Free();
		std::cout << "j" << std::endl;
		int result = last[last.size() - 1] + 1;
		std::vector<int> compacted(result);
		for (int i = 0; i < last.size(); ++i) {
			if (middle[i]) {
				compacted[last[i]] = vector[i];
			}
		}

		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		std::cout << duration.count() << " milliseconds" << std::endl;
	}
	MPI::Finalize();
	return 0;
}