#include "mpi.h"
#include <iostream>
#include <chrono>
#include <random>

void add(int*, int*, int*, MPI_Datatype*);
void add(int* invec, int* inoutvec, int* len, MPI_Datatype* dtype) {
	std::cout << "1" << std::endl;
	for (int i = 0; i < *len; i++) {
		std::cout << invec[i] << std::endl;
		inoutvec[i] += invec[i];
	}
}


int main()
{
	MPI::Init();
	std::cout << "f" << std::endl;
	int rank = MPI::COMM_WORLD.Get_rank();
	int cluster = MPI::COMM_WORLD.Get_size();
	const unsigned int size = 1000;
	const float vmax = 100;


	std::vector<float> values;
	std::vector<int> column;
	std::vector<int> row;
	std::vector<float> vector(size);
	int sizes[4];

	if (rank == 0) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<float> val(0, vmax);
		std::bernoulli_distribution dist(0.01);
		float temp[size][size];
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				if (dist(gen)) {
					temp[i][j] = val(gen);
				} else {
					temp[i][j] = 0;
				}
			}
		}
		for (int j = 0; j < size; ++j) {
				vector[j] = val(gen);
		}
		for (int i = 0; i < size; ++i) {
			row.push_back(values.size());
			for (int j = 0; j < size; ++j) {
				float var = temp[i][j];
				if (var != 0) {
					values.push_back(var);
					column.push_back(i);
				}
			}
		}
		sizes[0] = values.size();
		sizes[1] = column.size();
		sizes[2] = row.size();
		sizes[3] = vector.size();
	}
	MPI::COMM_WORLD.Bcast(sizes, 4, MPI_INT, 0);
	MPI::COMM_WORLD.Bcast(values.data(), sizes[0], MPI_FLOAT, 0);
	MPI::COMM_WORLD.Bcast(column.data(), sizes[1], MPI_FLOAT, 0);
	MPI::COMM_WORLD.Bcast(row.data(), sizes[2], MPI_FLOAT, 0);
	MPI::COMM_WORLD.Bcast(vector.data(), sizes[3], MPI_FLOAT, 0);

	std::cout << values.size() << std::endl;
	std::cout << column.size() << std::endl;
	std::cout << row.size() << std::endl;
	std::cout << vector.size() << std::endl;
	/*
	auto start = std::chrono::high_resolution_clock::now();
	std::vector<int> work(chunk);
	MPI::COMM_WORLD.Scatter(vector.data(), chunk, MPI_INT32_T, work.data(), chunk, MPI_INT32_T, 0);
	std::vector<int> out(chunk);
	for (int i = 0; i < chunk; ++i) {
		if (work[i] < barrier) {
			out[i] = 1;
		}
		else {
			out[i] = 0;
		}
	}
	std::vector<int> middle(size + extra);
	MPI::COMM_WORLD.Gather(out.data(), chunk, MPI_INT32_T, middle.data(), chunk, MPI_INT32_T, 0);
	MPI::Op* op = new MPI::Op();
	op->Init((MPI::User_function*)add, true);
	std::vector<int> last(size + extra);

	MPI::COMM_WORLD.Exscan(middle.data(), last.data(), middle.size(), MPI_INT32_T, *op);
	op->Free();

	if (rank == 0) {
		std::cout << "m" << std::endl;
		int result = last[last.size() - 1] + 1;
		std::cout << "l" << std::endl;
		std::vector<int> compacted(result);
		for (int i = 0; i < last.size(); ++i) {
			std::cout << i << std::endl;
			if (middle[i]) {
				compacted[last[i]] = vector[i];
			}
		}

		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		std::cout << duration.count() << " milliseconds" << std::endl;
	}
	MPI::Finalize();*/
	return 0;
}