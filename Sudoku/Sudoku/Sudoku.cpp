
#include "mpi.h"
#include <iostream>
#include <chrono>
#include <vector>

#include "Solver.h"

const int base = 2;
const int size = 4;

bool isAllowed(std::vector<int> board, int x, int y, int digit) {
	bool allowed = true;

	if (board[y * size + x] != 0) {
		allowed = false;
	}

	// Azonos sorban vagy oszlopban csak egy 'val' lehet
	for (int i = 0; allowed && i < 9; ++i) {
		if (board[y * size + i] == digit) {
			allowed = false;
		}
		if (board[i*size+x] == digit) {
			allowed = false;
		}
	}

	// Az adott cellában csak egy 'val' lehet
	int upperY = base * (y / base);
	int leftX = base * (x / base);

	for (int i = upperY; allowed && i < upperY + 3; ++i) {
		for (int j = leftX; allowed && j < leftX + 3; ++j) {
			if (board[i * size + j] == digit) {
				allowed = false;
			}
		}
	}

	return allowed;
}

bool isSolved(std::vector<int> board) {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			if (board[i * size + j] == 0) {
				return true;
			}
		}
	}
	return false;
}

void solveBack(std::vector<int> board) {
	if (isSolved(board)) {
		std::cout << "Home" << std::endl;
		MPI::COMM_WORLD.Send(board.data(), size * size, MPI_INT, 0, 0);
	} else {
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				for (int k = 1; k <= size; ++k) {
					if (isAllowed(board, i, j, k)) {
						std::vector<int> t(board);
						t[i * size + j] = k;
						solveBack(t);
					}
				}
			}
		}
	}
}

int main()
{
	MPI::Init();
	std::cout << "f" << std::endl;
	int rank = MPI::COMM_WORLD.Get_rank();
	int cluster = MPI::COMM_WORLD.Get_size();

	std::vector<int> initial = { 0,0,0,8,0,1,0,0,0, 0,0,0,0,0,0,0 };//, 4, 3, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 8, 0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 6, 0, 0, 0, 0, 0, 0, 7, 5, 0, 0, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 6, 0, 0};
	std::vector<std::vector<int> > beginnings;
	std::vector<std::vector<int> > solutions;

	if (rank == 0) {
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				for (int k = 1; k <= size; ++k) {
					if (isAllowed(initial, i, j, k)) {
						std::vector<int> t(initial);
						t[i * size + j] = k;
						beginnings.push_back(t);
					}
				}
			}
		}
		std::cout << "First loop?" << beginnings.size() << std::endl;
		int next = 0;
		int done = 0;
		MPI::Status* status = new MPI::Status();
		while (next < beginnings.size()) {
			std::vector<int> buf(size * size);
			MPI::COMM_WORLD.Recv(buf.data(), size * size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, *status);
			std::cout << "Prompt" << std::endl;
			int tag = status->Get_tag();
			if (tag == 0) {
				std::cout << "Get good." << std::endl;
				solutions.push_back(buf);
			} else {
				std::cout << "Send work." << std::endl;
				MPI::COMM_WORLD.Send(beginnings[next++].data(), size * size, MPI_INT, status->Get_source(), 0);
			}
		}
		while (done < cluster -1) {
			std::vector<int> buf(size * size);
			std::cout << "Beginning of the end." << std::endl;
			MPI::COMM_WORLD.Recv(buf.data(), size * size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, *status);
			int tag = status->Get_tag();
			if (tag == 0) {
				solutions.push_back(buf);
			} else {
				MPI::COMM_WORLD.Send(nullptr, 0, MPI_INT, status->Get_source(), 0);
				++done;
			}
		}
		delete status;

		if (solutions.empty()) {
			std::cout << "No solutions found." << std::endl;
		} else {
			std::cout << "Solutions: ";
			for (std::vector<int> solution : solutions) {
				for (int i = 0; i < size * size; ++i) {
					std::cout << solution[i];
				}
				std::cout << std::endl;
			}
		}
		
	} else {
		std::vector<int> work(size*size);
		bool run = true;
		MPI::Status *status = new MPI::Status();
		while (run) {
			MPI::COMM_WORLD.Send(nullptr, 0, MPI_INT, 0, 1);
			std::cout << rank << "A." << std::endl;
			MPI::COMM_WORLD.Recv(work.data(), size * size, MPI_INT, 0, MPI_ANY_TAG, *status);
			int tag = status->Get_tag();
			if (tag == 0) {
				std::cout << rank << " Work: ";
				for (int i = 0; i < size * size; ++i) {
					std::cout << work[i];
				}
				std::cout << std::endl;
				solveBack(work);
				std::cout << rank << "V" << std::endl;
			} else {
				run = false;
			}
		}
		delete status;
	}

	MPI::Finalize();
    return 0;
}
