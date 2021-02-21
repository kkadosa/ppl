// Fractal.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "mpi.h"
#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <complex>
#include <chrono>

void WriteTGA_RGB(const char* filename, unsigned char* data, unsigned int width, unsigned int height)
{
	FILE* f = fopen(filename, "wb");
	if (!f) {
		fprintf(stderr, "Unable to create output TGA image `%s'\n", filename);
		exit(EXIT_FAILURE);
	}

	fputc(0x00, f); /* ID Length, 0 => No ID        */
	fputc(0x00, f); /* Color Map Type, 0 => No color map included   */
	fputc(0x02, f); /* Image Type, 2 => Uncompressed, True-color Image */
	fputc(0x00, f); /* Next five bytes are about the color map entries */
	fputc(0x00, f); /* 2 bytes Index, 2 bytes length, 1 byte size */
	fputc(0x00, f);
	fputc(0x00, f);
	fputc(0x00, f);
	fputc(0x00, f); /* X-origin of Image    */
	fputc(0x00, f);
	fputc(0x00, f); /* Y-origin of Image    */
	fputc(0x00, f);
	fputc(width & 0xff, f); /* Image Width      */
	fputc((width >> 8) & 0xff, f);
	fputc(height & 0xff, f); /* Image Height     */
	fputc((height >> 8) & 0xff, f);
	fputc(0x18, f); /* Pixel Depth, 0x18 => 24 Bits */
	fputc(0x20, f); /* Image Descriptor     */

	for (int y = height - 1; y >= 0; y--) {
		for (size_t x = 0; x < width; x++) {
			const size_t i = (y * width + x) * 3;
			fputc(data[i + 2], f); /* write blue */
			fputc(data[i + 1], f); /* write green */
			fputc(data[i], f); /* write red */
		}
	}
}

int main()
{
	MPI::Init();
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int cluster = MPI::COMM_WORLD.Get_size();
	std::cout << rank << " c" << std::endl;
	const unsigned int size = 1000;
	if (rank == 0) {
		auto start = std::chrono::high_resolution_clock::now();
		unsigned char* data = new unsigned char[size * size * 3];
		std::memset(data, 0, size * size * 3 * sizeof(unsigned char));

		int finished = 1;
		unsigned int* buf = new unsigned int[3];
		MPI::Status* status = new MPI::Status();
		while (finished < cluster) {
			MPI::COMM_WORLD.Recv(buf, 3, MPI_UINT32_T, MPI_ANY_SOURCE, MPI_ANY_TAG, *status);
			int tag = status->Get_tag();
			if (tag == 0) {
				unsigned char* start = data + 3 * (buf[0] * size + buf[1]);
				//std::cout << "teve " << buf[0] << " " << buf[1] << " " << buf[2] << std::endl;
				std::memset(start, -1, 3 * (buf[2] - buf[1]) * sizeof(unsigned char));
			}
			else if (tag == 1) {
				++finished;
				std::cout << "fin: " << status->Get_source() << std::endl;
			}
		}
		WriteTGA_RGB("mandelbrot.tga", data, size, size);
		delete status;
		delete[] data;
		delete[] buf;
		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		std::cout << duration.count() << " milliseconds" << std::endl;
		MPI::Finalize();
	} else {
		std::complex<double> center(-1.68, -1.23);
		double scale = 2.35;
		const unsigned int maxIterations = 100;
		std::cout << rank << " a" << std::endl;
		for (unsigned int y = (rank - 1); y < size; y += (cluster - 1)) {
			bool started = false;
			unsigned int begin;
			std::cout << rank << " l" << std::endl;
			for (unsigned int x = 0; x < size; ++x) {
				bool black = false;

				std::complex<double> c(x / (double)size * scale + center.real(),
					y / (double)size * scale + center.imag());
				std::complex<double> z(c);

				for (unsigned int iteration = 0; iteration < maxIterations && !black; ++iteration) {
					z = z * z + c;
					black = std::abs(z) > 1.0f;
				}

				if (started) {
					if (!black) {
						unsigned int buf[3] = { y, begin, x };
						started = false;
						MPI::COMM_WORLD.Send(buf, 3, MPI_UINT32_T, 0, 0);
					}
				} else {
					if (black) {
						started = true;
						begin = x;
					}
				}
			}
			if (started) {
				unsigned int buf[3] = { y, begin, size };
				MPI::COMM_WORLD.Send(buf, 3, MPI_UINT32_T, 0, 0);
			}
		}
		unsigned int buf[3] = { 1, 2, 3};
		MPI::COMM_WORLD.Send(buf, 3, MPI_UINT32_T, 0, 1);
		MPI::Finalize();
	}
	return 0;
}