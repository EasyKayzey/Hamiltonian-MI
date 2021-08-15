#ifndef ME_EX_MAIN_H
#define ME_EX_MAIN_H

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <ctime>
#include <complex>
#include <functional>
#include <atomic>
#include <chrono>
#include <unordered_set>
#include <iomanip>
#include <atomic>
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "Eigen/LU"
#include "unsupported/Eigen/MatrixFunctions"
#include "omp.h"

#if !(defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64))
#include <filesystem>
#define HAS_FILESYSTEM 1
#endif

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define path "C:\\EZKZ\\ME-G1\\"
#elif defined(linux) || defined(_linux_) || defined(unix) || defined(_unix_) || defined(__unix__)
#define path "./"
#endif

using namespace std;
using namespace Eigen;

const int DIM = 6;
const int L = (DIM * (DIM - 1)) / 2;
const double HBAR = 1;
const int N_TO = 2;
const int N_OBS = DIM * N_TO;
const int N_H = L;

typedef complex<double> Complex;
typedef Matrix<Complex, DIM, DIM> EMatrix;
typedef DiagonalMatrix<Complex, DIM> EDMatrix;
typedef Matrix<Complex, DIM, 1> EVector;
typedef Matrix<Complex, 1, DIM> ECovector;
typedef Matrix<Complex, Dynamic, 1> TVector;
typedef Matrix<double, Dynamic, 1> RTVector;
typedef array<Complex, N_OBS> OArr;
typedef array<double, DIM> DArr;
typedef array<pair<double, double>, L> FGenome;
typedef mt19937 rng;

int main(int argc, char** argv);

double envelope_funct(double t);

EMatrix to_full_matrix(EMatrix upper);

pair<pair<EMatrix, EMatrix>, EVector> diag_vec(const EMatrix& mu, const EDMatrix& C);

OArr evolve_initial_nonhermitian(const vector<double>& epsilon, const EMatrix& mu, const EVector& psi_i, const array<ECovector, DIM>& anal_pop);

Complex get_only_element(Matrix<Complex, -1, -1> scalar);

pair<int, int> calc_loc(int u_i);

// vector<DArr> gen_pop_graphs(const vector<double>& eps_inter, const EMatrix& CP, const EMatrix& PdC,
//                            const EVector& lambda, const EVector& psi_i, const array<ECovector, DIM>& anal_pop);

vector<double> get_field(const FGenome& epsilon);

void ptime();

template <class T, class F> void print_vec(vector<vector<T>> vec, ofstream& outfile, F lambda);
template <class T, class F, size_t N> void print_arr(vector<array<T, N>> vec, ofstream& outfile, F lambda);

inline double normalize(double rand, double min, double max) {
    return rand * (max - min) + min;
}

// Here's a hash function for arrays:
template<typename T, size_t N>
struct hash<array<T, N>>
{
    typedef array<T, N> argument_type;
    typedef size_t result_type;

    result_type operator()(const argument_type& a) const
    {
        hash<T> hasher;
        result_type h = 0;
        for (result_type i = 0; i < N; ++i)
        {
            h = h * 31 + hasher(a[i]);
        }
        return h;
    }
};

// FFT stuff

#include <valarray>
typedef valarray<Complex> CArray;

void fft(CArray &x)
{
	// DFT
	unsigned int N = x.size(), k = N, n;
	double thetaT = 3.14159265358979323846264338328L / N;
	Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				Complex t = x[a] - x[b];
				x[a] += x[b];
				x[b] = t * T;
			}
			T *= phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			Complex t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}
	// // Normalize (This section make it not working correctly)
	// Complex f = 1.0 / sqrt(N);
	// for (unsigned int i = 0; i < N; i++)
	// 	x[i] *= f;
}
 
// inverse fft (in-place)
void ifft(CArray& x)
{
    // conjugate the complex numbers
    x = x.apply(std::conj);
 
    // forward fft
    fft( x );
 
    // conjugate the complex numbers again
    x = x.apply(std::conj);
 
    // scale the numbers
    x /= x.size();
}

#endif