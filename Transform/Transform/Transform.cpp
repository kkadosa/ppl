#include <complex>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <chrono>

#include "lodepng.h"
#include "mpi.h"

#ifndef M_PI
constexpr auto  M_PI = 3.14159265358979323846;
#endif

void loadPNG(std::string filename, unsigned& w, unsigned& h, std::vector<unsigned char>& image)
{
	unsigned char* imageData;
	unsigned err = lodepng_decode_file(&imageData, &w, &h, filename.c_str(), LCT_GREY, 8);
	if (err != 0)
	{
		std::cout << "Image decoder error: " << err << std::endl;
		exit(-1);
	}
	image.resize(w * h);
	memcpy(&image[0], imageData, w * h);
}

void savePNG(std::string filename, unsigned& w, unsigned& h, std::vector<unsigned char>& image)
{
	unsigned err = lodepng_encode_file(filename.c_str(), image.data(), w, h, LCT_GREY, 8);
	if (err != 0)
	{
		std::cout << "Image encoder error: " << err << std::endl;
		exit(-1);
	}
}

void RealToComplex(std::vector<unsigned char>& in, std::vector<std::complex<double>>& out)
{
	out.clear();
	out.reserve(in.size());
	for (auto r : in)
	{
		out.push_back(std::complex<double>(static_cast<double>(r), 0.0));
	}
}

void ComplexToReal(std::vector<std::complex<double>>&in, std::vector<unsigned char>& out)
{
	out.clear();
	out.reserve(in.size());
	for (auto c : in)
	{
		double length = sqrt(c.real() * c.real() + c.imag() * c.imag());
		out.push_back(length);
	}
}

/*
void DFT(std::vector<std::complex<double>>& in, std::vector<std::complex<double>>& out, unsigned w, unsigned h, bool horizontal, bool inverse)
{
	out.clear();
	out.resize(in.size());

	for (unsigned i = 0; i < h; ++i)
	{
		for (unsigned k = 0; k < w; ++k)
		{
			std::complex<double> sum(0.0, 0.0);
			for (unsigned j = 0; j < w; ++j)
			{
				size_t addr = horizontal ? j + i * w : i + j * w;
				auto angle = (inverse ? -2.0f : 2.0) * M_PI * j * k / w;
				sum.real(sum.real() + in[addr].real() * cos(angle) - in[addr].imag() * sin(angle));
				sum.imag(sum.imag() + in[addr].real() * sin(angle) + in[addr].imag() * cos(angle));
			}

			if (!inverse)
			{
				sum *= 1.0 / w;
			}

			out[k + i * w] = sum;
		}
		std::cout << i << std::endl;
	}
}
*/

void transform(std::vector<std::complex<double>>& in, std::vector<std::complex<double>>& out, unsigned width, bool inverse, unsigned sendn) {
	out.clear();
	out.resize(in.size());

	for (unsigned i = 0; i < sendn; ++i) {
		for (unsigned k = 0; k < width; ++k) {

			std::complex<double> sum(0.0, 0.0);
			for (unsigned j = 0; j < width; ++j) {
				size_t addr = j + i * width;
				auto angle = (inverse ? -2.0f : 2.0) * M_PI * j * k / width;
				sum.real(sum.real() + in[addr].real() * cos(angle) - in[addr].imag() * sin(angle));
				sum.imag(sum.imag() + in[addr].real() * sin(angle) + in[addr].imag() * cos(angle));
			}
			if (!inverse) {
				sum *= 1.0 / width;
			}
			out[k + i * width] = sum;
		}
		std::cout << i << std::endl;
	}
}

void t(std::vector<std::complex<double>>& in, std::vector<std::complex<double>>& out, unsigned width) {
	for (int i = 0; i < width; ++i) {
		for (int j = 0; j < width; ++j) {
			out[i + j * width] = in[i * width + j];
		}
	}
}

int main()
{
	MPI::Init();
	int rank = MPI::COMM_WORLD.Get_rank();
	int cluster = MPI::COMM_WORLD.Get_size();

	unsigned width = 0, h = 0;
	std::vector<unsigned char> image;
	std::string inName("lena.png");
	std::string dftName("lenaDFT.png");
	std::string outName("lenaOut.png");
	if (rank == 0) {
		loadPNG(inName, width, h, image);
	}

	for (int l = 0; l < 1; ++l) {
		auto start = std::chrono::high_resolution_clock::now();

		MPI::COMM_WORLD.Bcast(&width, 1, MPI_UNSIGNED, 0);
		unsigned sendn = width / cluster;
		MPI::Datatype LINE = MPI::DOUBLE.Create_contiguous(width * 2);

		std::vector<std::complex<double>> full1;
		if (rank == 0) {
			RealToComplex(image, full1);
		}
		std::vector<std::complex<double>> in(width * sendn);
		std::vector<std::complex<double>> out;

		MPI::COMM_WORLD.Scatter(full1.data(), sendn, LINE, in.data(), sendn, LINE, 0);
		transform(in, out, width, false, sendn);
		std::cout << in.size() << " " << out.size() << std::endl;
		MPI::COMM_WORLD.Gather(out.data(), sendn, LINE, full1.data(), sendn, LINE, 0);
		std::vector<std::complex<double>> full2(full1.size());
		if (rank == 0) {
			t(full1, full2, width);
		}
		MPI::COMM_WORLD.Scatter(full2.data(), sendn, LINE, in.data(), sendn, LINE, 0);
		transform(in, out, sendn, false, cluster);
		MPI::COMM_WORLD.Gather(out.data(), sendn, LINE, full2.data(), sendn, LINE, 0);
		if (rank == 0) {
			t(full2, full1, width);
			ComplexToReal(full1, image);
			savePNG(dftName, width, width, image);
			std::cout << "DFT finished" << std::endl;
		}

		MPI::COMM_WORLD.Scatter(full1.data(), sendn, LINE, in.data(), sendn, LINE, 0);
		transform(in, out, width, true, sendn);
		MPI::COMM_WORLD.Gather(out.data(), sendn, LINE, full1.data(), sendn, LINE, 0);
		if (rank == 0) {
			t(full1, full2, width);
		}
		MPI::COMM_WORLD.Scatter(full2.data(), sendn, LINE, in.data(), sendn, LINE, 0);
		transform(in, out, sendn, true, cluster);
		MPI::COMM_WORLD.Gather(out.data(), sendn, LINE, full2.data(), sendn, LINE, 0);
		if (rank == 0) {
			t(full2, full1, width);
			ComplexToReal(full1, image);
			savePNG(dftName, width, width, image);
			std::cout << "IDFT finished" << std::endl;
		}

		if (rank == 0) {
			auto end = std::chrono::high_resolution_clock::now();
			std::cout << cluster << " " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
		}
	}

	/*
	std::vector<std::complex<double>> f1;
	std::vector<std::complex<double>> f2;
	RealToComplex(image, f1);
	DFT(f1, f2, w, h, true, false);
	DFT(f2, f1, w, h, false, false);
	std::cout << "DFT finished" << std::endl;
	ComplexToReal(f1, image);
	savePNG(dftName, w, h, image);

	DFT(f1, f2, w, h, true, true);
	DFT(f2, f1, w, h, false, true);
	std::cout << "IDFT finished" << std::endl;
	ComplexToReal(f1, image);
	savePNG(outName, w, h, image);
	*/
	return 0;
}

