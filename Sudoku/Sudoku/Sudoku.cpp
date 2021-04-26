
#include "mpi.h"
#include <iostream>
#include <chrono>
#include <vector>

#include "Solver.h"

const int base = 3;
const int size = 9;

bool isAllowed(const std::vector<int>& board, int x, int y, int digit) {
	// Azonos sorban vagy oszlopban csak egy 'val' lehet
	for (int i = 0; i < size; ++i) {
		if (board[y * size + i] == digit) {
			return false;
		}
		if (board[i*size+x] == digit) {
			return false;
		}
	}

	// Az adott cellában csak egy 'val' lehet
	int upperY = base * (int)(y / base);
	int leftX = base * (int)(x / base);

	for (int i = upperY; i < upperY + base; ++i) {
		for (int j = leftX; j < leftX + base; ++j) {
			if (board[i * size + j] == digit) {
				return false;
			}
		}
	}

	return true;
}

bool isSolved(const std::vector<int>& board) {
	for (int i = 0; i < size * size; ++i) {
		if (board[i] == 0) {
			return false;
		}
	}
	return true;
}

void solveBack(const std::vector<int>& board) {
	
	for (int i = 0; i < size * size; ++i) {
		std::cout << board[i];
	}
	std::cout << std::endl;

	if (isSolved(board)) {
		MPI::COMM_WORLD.Send(board.data(), size * size, MPI_INT, 0, 0);
	} else {
		bool go = true;
		for (int y = 0; go && y < size; ++y) {
			for (int x = 0; go && x < size; ++x) {
				if (board[y * size + x] == 0) {
					std::vector<int> possibilities;

					for (int k = 1; k <= size; ++k) {
						if (isAllowed(board, x, y, k)) {
							possibilities.push_back(k);
						}
					}

					if (possibilities.empty()) {
						go = false;
						std::cout << "up" << std::endl;
					} else {
						for (int k : possibilities) {
							std::vector<int> t(board);
							std::cout << y << " " << x <<" " <<  k << std::endl;
							t[y * size + x] = k;
							solveBack(t);
						}
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

	std::vector<int> initial = { 0,0,0,8,0,1,0,0,0, 0,0,0,0,0,0,0 , 4, 3, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 8, 0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 6, 0, 0, 0, 0, 0, 0, 7, 5, 0, 0, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 6, 0, 0 };
	std::vector<std::vector<int> > beginnings;
	std::vector<std::vector<int> > solutions;

	if (rank == 0) {
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				for (int k = 1; k <= size; ++k) {
					if (isAllowed(initial, j, i, k)) {
						std::vector<int> t(initial);
						t[i * size + j] = k;
						beginnings.push_back(t);
					}
				}
			}
		}
		int next = 0;
		int done = 0;
		MPI::Status* status = new MPI::Status();
		while (next < beginnings.size()) {
			std::vector<int> buf(size * size);
			MPI::COMM_WORLD.Recv(buf.data(), size * size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, *status);
			int tag = status->Get_tag();
			if (tag == 0) {
				solutions.push_back(buf);
			} else {
				MPI::COMM_WORLD.Send(beginnings[next++].data(), size * size, MPI_INT, status->Get_source(), 0);
			}
		}
		while (done < cluster -1) {
			std::vector<int> buf(size * size);
			MPI::COMM_WORLD.Recv(buf.data(), size * size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, *status);
			int tag = status->Get_tag();
			if (tag == 0) {
				solutions.push_back(buf);
			} else {
				MPI::COMM_WORLD.Send(nullptr, 0, MPI_INT, status->Get_source(), 1);
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
			MPI::COMM_WORLD.Recv(work.data(), size * size, MPI_INT, 0, MPI_ANY_TAG, *status);
			int tag = status->Get_tag();
			if (tag == 0) {
				std::cout << "Work for " << rank << ": ";
				for (int i = 0; i < size * size; ++i) {
					std::cout << work[i];
				}
				std::cout << std::endl;
				solveBack(work);
				std::cout << "Work done " << rank << std::endl;
			} else {
				run = false;
			}
		}
		delete status;
	}

	MPI::Finalize();
    return 0;
}
