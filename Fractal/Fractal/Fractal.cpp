// Fractal.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "mpi.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <complex>

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
	const unsigned int size = 1000;
	if (rank == 0) {
		unsigned char* data = new unsigned char[size * size * 3];
		std::memset(data, 0, size * size * 3 * sizeof(unsigned char));

		int finished = 1;

		while (finished < cluster) {
			unsigned int buf[3];
			MPI::Status status;
			MPI::COMM_WORLD.Recv(buf, 3, MPI_UINT32_T, MPI_ANY_SOURCE, MPI_ANY_TAG, status);
			int tag = status.Get_tag();
			if (tag == 0) {
				unsigned char* start = data + 3 * buf[0] * size + 3 * buf[1];
				std::memset(start, -1, 3 * buf[2] - buf[1] * sizeof(unsigned char));
			}
			else if (tag == 1) {
				++finished;
			}
			WriteTGA_RGB("mandelbrot.tga", data, size, size);
			delete[] data;
		}
	} else {
		std::complex<double> K(0.353, 0.288);
		std::complex<double> center(-1.68, -1.23);
		double scale = 2.35;
		const unsigned int maxIterations = 100;

		for (unsigned int y = rank - 1; y < size; y += rank - 1) {
			bool started = false;
			unsigned int begin;
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
	}
	return 0;
}