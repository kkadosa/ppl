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
	const unsigned int size = 100000;
	const float vmax = 10;


	std::vector<float> valuesV;
	std::vector<int> columnV;
	std::vector<int> rowV;
	float* vector = new float[size];

	int sizes[3];
	float* values;
	int* column;
	int* row;

	if (rank == 0) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<float> val(0, vmax);
		std::bernoulli_distribution dist(0.01);

		float** temp = new float*[size];
		for (int j = 0; j < size; ++j) {
			temp[j] = new float[size];
		}

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
			rowV.push_back(valuesV.size());
			for (int j = 0; j < size; ++j) {
				float var = temp[i][j];
				if (var != 0) {
					valuesV.push_back(var);
					columnV.push_back(i);
				}
			}
		}
		rowV.push_back(size);

		sizes[0] = valuesV.size();
		sizes[1] = columnV.size();
		sizes[2] = rowV.size();
		values = valuesV.data();
		column = columnV.data();
		row = rowV.data();

		for (int j = 0; j < size; ++j) {
			delete[] temp[j];
		}
		delete[] temp;
	}

	//BEGIN
	
		auto start = std::chrono::high_resolution_clock::now();

		MPI::COMM_WORLD.Bcast(sizes, 3, MPI_INT, 0);
		if (rank != 0) {
			values = new float[sizes[0]];
			column = new int[sizes[1]];
			row = new int[sizes[2]];
		}
		MPI::COMM_WORLD.Bcast(values, sizes[0], MPI_FLOAT, 0);
		MPI::COMM_WORLD.Bcast(column, sizes[1], MPI_INT, 0);
		MPI::COMM_WORLD.Bcast(row, sizes[2], MPI_INT, 0);
		MPI::COMM_WORLD.Bcast(vector, size, MPI_FLOAT, 0);

		float* result = new float[size];
		for (int i = rank; i < size; i += cluster) {
			float accumulator = 0;
			int begin = row[i];
			int end = row[i + 1];
			for (int j = begin; j < end; ++j) {
				int left = values[j];
				int right = vector[column[j]];
				accumulator += left * right;
			}
			result[i] = accumulator;
		}

		float* receive;
		if (rank == 0) {
			receive = new float[size * cluster];
		}

		MPI::COMM_WORLD.Gather(result, size, MPI_FLOAT, receive, size, MPI_FLOAT, 0);

		if (rank == 0) {
			float* out = new float[size];
			for (int i = 0; i < size; i += cluster) {
				for (int j = 0; j < cluster && i + j < size; ++j) {
					int index = j * size + i;
					out[i + j] = receive[index];
				}
			}
			//END
			auto end = std::chrono::high_resolution_clock::now();
			std::cout << cluster << ", " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
			delete[] out;
		}


		delete[] vector;
		delete[] result;
		if (rank != 0) {
			delete[] values;
			delete[] column;
			delete[] row;
		}
		else {
			delete[] receive;
		}
	
	MPI::Finalize();
	return 0;
}