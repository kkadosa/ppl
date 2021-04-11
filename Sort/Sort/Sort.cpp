
#include "mpi.h"
#include <iostream>
#include <chrono>
#include <random>

void merge(int* input, int leftEnd, int middle, int rightEnd) {
	int leftSize = middle - leftEnd + 1;
	int rightSize = rightEnd - middle;

	int* left = new int[leftSize];
	int* right = new int[rightSize];
	for (int i = 0; i < leftSize; ++i) {
		left[i] = input[leftEnd + i];
	}
	for (int i = 0; i < rightSize; ++i) {
		right[i] = input[middle + 1 + i];
	}
	int i = 0;
	int j = 0;
	int k = leftEnd;
	while (i < leftSize && j < rightSize) {
		if (left[i] < right[j]) {
			input[k++] = left[i++];
		} else {
			input[k++] = right[j++];
		}
	}
	while (i < leftSize) {
		input[k++] = left[i++];
	}
	while (j < rightSize) {
		input[k++] = right[j++];
	}
	delete[] left;
	delete[] right;
}

void sort(int* input, int leftEnd, int rightEnd) {
	if (leftEnd < rightEnd) {
		int middle = leftEnd + (rightEnd - leftEnd) / 2;
		sort(input, leftEnd, middle);
		sort(input, middle + 1, rightEnd);
		merge(input, leftEnd, middle, rightEnd);
	}
}

void mergesort(int* input, int size, int rank, int cluster) {
	auto start = std::chrono::high_resolution_clock::now();
	int base = size / cluster;
	int boost = size - (base * (cluster - 1));
	int* sizes = new int[cluster];
	sizes[0] = boost;
	for (int i = 1; i < cluster; ++i) {
		sizes[i] = base;
	}
	int* displacement = new int[cluster];
	displacement[0] = 0;
	for (int i = 1; i < cluster; ++i) {
		displacement[i] = displacement[i - 1] + sizes[i - 1];
	}

	int* mine = new int[sizes[rank]];
	MPI::COMM_WORLD.Scatterv(input, sizes, displacement, MPI_INT, mine, sizes[rank], MPI_INT, 0);
	sort(mine, 0, sizes[rank]);
	MPI::COMM_WORLD.Gatherv(mine, sizes[rank], MPI_INT, input, sizes, displacement, MPI_INT, 0);

	int prevthreads = cluster;
	int threads = cluster / 2;
	int* prevsizes = sizes;
	while (threads > 3) {
		sizes = new int[cluster];
		displacement = new int[cluster];
		if (prevthreads % 2 == 1) {
			sizes[0] = prevsizes[0] + prevsizes[1] + prevsizes[2];
			for (int i = 1; i < threads; ++i) {
				sizes[i] = prevsizes[2 * i + 1] + prevsizes[2 * i + 2];
			}
			for (int i = threads; i < cluster; ++i) {
				sizes[i] = 0;
			}
			displacement[0] = 0;
			for (int i = 1; i < cluster; ++i) {
				displacement[i] = displacement[i - 1] + sizes[i - 1];
			}
			if (sizes[rank] > 0) {
				delete[] mine;
				int* mine = new int[sizes[rank]];
			}
			MPI::COMM_WORLD.Scatterv(input, sizes, displacement, MPI_INT, mine, sizes[rank], MPI_INT, 0);
			if (sizes[rank] > 0) {
				if (rank == 0) {
					merge(mine, 0, prevsizes[0] - 1, prevsizes[0] + prevsizes[1] - 1);
					merge(mine, 0, prevsizes[0] + prevsizes[1] - 1, prevsizes[0] + prevsizes[1] + prevsizes[2] - 1);
				} else {
					merge(mine, 0, prevsizes[2 * rank + 1] - 1, sizes[rank] - 1);
				}
			}
			MPI::COMM_WORLD.Gatherv(mine, sizes[rank], MPI_INT, input, sizes, displacement, MPI_INT, 0);
		} else {
			for (int i = 0; i < threads; ++i) {
				sizes[i] = prevsizes[2 * i] + prevsizes[2 * i + 1];
			}
			for (int i = threads; i < cluster; ++i) {
				sizes[i] = 0;
			}
			displacement[0] = 0;
			for (int i = 1; i < cluster; ++i) {
				displacement[i] = displacement[i - 1] + sizes[i - 1];
			}
			if (sizes[rank] > 0) {
				delete[] mine;
				int* mine = new int[sizes[rank]];
				MPI::Intracomm comm = MPI::COMM_WORLD.Split(0, rank);
				std::cout << rank << " s" << std::endl;
				comm.Scatterv(input, sizes, displacement, MPI_INT, mine, sizes[rank], MPI_INT, 0);
				std::cout << rank << " m" << std::endl;
				merge(mine, 0, prevsizes[2 * rank] - 1, sizes[rank] - 1);
			} else {
				MPI::COMM_WORLD.Split(MPI_UNDEFINED, rank);
				std::cout << rank << " n" << std::endl;
			}
			std::cout << rank << " f" << std::endl;
			MPI::COMM_WORLD.Gatherv(mine, sizes[rank], MPI_INT, input, sizes, displacement, MPI_INT, 0);
		}
		prevthreads = threads;
		threads = threads / 2;
		delete[] prevsizes;
		prevsizes = sizes;
	}
	delete[] displacement;
	delete[] mine;

	if (rank == 0) {
		if (threads == 3) {
			merge(mine, 0, sizes[0] - 1, sizes[0] + sizes[1] - 1);
			merge(mine, 0, sizes[0] + sizes[1] - 1, size - 1);
		} else {
			merge(input, 0, sizes[0] - 1, size - 1);
		}
		auto end = std::chrono::high_resolution_clock::now();
		std::cout << cluster << ", " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << ", ms" << std::endl;
	}
	delete[] sizes;
}

int main()
{
	MPI::Init();
	std::cout << "f" << std::endl;
	int rank = MPI::COMM_WORLD.Get_rank();
	int cluster = MPI::COMM_WORLD.Get_size();

	const unsigned int size1 = 1000000;
	const int vmax = 256;

	int* input;

	if (rank == 0) {
		std::cout << "f1" << std::endl;
		input = new int[size1];
		std::cout << "f2" << std::endl;
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<int> val(0, vmax);

		for (int i = 0; i < size1; ++i) {
			input[i] = val(gen);
		}
	}

	for (int i = 0; i < 40; ++i) {
		mergesort(input, size1, rank, cluster);
	}

	const unsigned int size2 = 100000000;

	if (rank == 0) {
		delete[] input;
		input = new int[size2];
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<int> val(0, vmax);

		for (int i = 0; i < size2; ++i) {
			input[i] = val(gen);
		}
	}

	for (int i = 0; i < 40; ++i) {
		mergesort(input, size2, rank, cluster);
	}
	if (rank == 0) {
		delete[] input;
	}
	MPI::Finalize();
	return 0;
}
