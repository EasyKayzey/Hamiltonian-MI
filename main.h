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
#include <valarray>
#include <map>
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "Eigen/LU"
#include "unsupported/Eigen/KroneckerProduct"
#include "unsupported/Eigen/MatrixFunctions"
#include "omp.h"

#if !(defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64))
#include <filesystem>
#define HAS_FILESYSTEM 1
#endif

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define path "C:\\EZKZ\\HMI-1\\"
#elif defined(linux) || defined(_linux_) || defined(unix) || defined(_unix_) || defined(__unix__)
#define path "./"
#endif

using namespace std;
using namespace Eigen;

// User constants (change me!)
const int DIM = 4;
const int N_TO = 1;
const int N_FIELDS = 2;
const int BASE = 7;
const string field_dir = "fields/gaf-2/";

// Automatic constants and typedefs
const int L = (DIM * (DIM - 1)) / 2;
const int N_OBS = DIM * N_TO;
const int N_H = L;
const double HBAR = 1;
const double MY_PI = 3.14159265358979323846264338327950288419716939937510582097494459230781;

#define USE_LONG_DOUBLE // comment out to use default double (change me!)
#ifdef USE_LONG_DOUBLE
typedef long double DOUBLE_TYPE;
#else
typedef double DOUBLE_TYPE;
#endif
typedef complex<DOUBLE_TYPE> Complex;
typedef array<DOUBLE_TYPE, DIM> DArr;
typedef array<Complex, N_OBS> OArr;
typedef Matrix<Complex, DIM, DIM> EMatrix;
typedef DiagonalMatrix<Complex, DIM> EDMatrix;
typedef Matrix<Complex, DIM, 1> EVector;
typedef Matrix<Complex, 1, DIM> ECovector;
typedef Matrix<Complex, Dynamic, 1> TVector;
typedef Matrix<double, Dynamic, 1> RTVector;

typedef Matrix<Complex, 2, 2> EMatrix2;
typedef array<vector<double>, N_FIELDS> FieldSet;
typedef array<EMatrix, N_FIELDS> DipoleSet;
typedef valarray<Complex> CArray;
typedef mt19937 rng;

// Hamlitonian (change me!)
DipoleSet dipoles_upper;
EVector H0D;
void gen_hamiltonians() {

    EMatrix2 I2 = EMatrix2::Identity();
    EMatrix2 Sx, Sy, Sz;
    Sx << 0,   1,
          1,   0;
    Sy << 0, -1i,
          1i,  0;
    Sz << 1,   0,
          0,  -1;
    Sx /= 2, Sy /= 2, Sz /= 2;

    double omega_1 = 2* MY_PI *15;
    double omega_2 = -2* MY_PI *21;
    double J12 = 2* MY_PI * 1;
    H0D = (omega_1 * kroneckerProduct(Sz, I2).eval() + omega_2 * kroneckerProduct(I2, Sz).eval() + J12 * kroneckerProduct(Sz, Sz).eval()).diagonal();
    dipoles_upper[0] = (kroneckerProduct(Sx, I2) + kroneckerProduct(I2, Sx)).triangularView<Eigen::Upper>();
    dipoles_upper[1] = (kroneckerProduct(Sy, I2) + kroneckerProduct(I2, Sy)).triangularView<Eigen::Upper>();

//        double omega_1 = 0;
//    H0D = (-1 * omega_1 * 2 * Sz).diagonal();
//    dipoles_upper[0] = (2 * Sx).triangularView<Eigen::Upper>();
//    dipoles_upper[1] = (2 * Sy).triangularView<Eigen::Upper>();

    // dipoles_upper[0] <<  0,   0.061,  -.013,
    //                    0,   0,  .083,
    //                    0,  0,  0;
}

// Encoding (change me!)
enum enc_scheme { other, order, partial, full };
enum enc_type { hermitian, antihermitian, nonhermitian };

const enc_scheme cur_scheme = partial;
const enc_type cur_type = nonhermitian;

EMatrix get_partial_encoding_integers() {
	EMatrix partial_encoding_integers;
//	partial_encoding_integers <<
//		0, 0,
//		1, 0;
    // partial_encoding_integers <<
    //         0, 2, 3,
    //         0, 0, 4,
    //         1, 0, 0;
//        partial_encoding_integers <<
//                0, 0, 0, 0,
//                0, 0, 0, 0,
//                1, 0, 0, 0,
//                0, 0, 0, 0;
    partial_encoding_integers <<
            0, 2, 3, 0,
            0, 0, 0, 4,
            1, 0, 0, 5,
            0, 0, 0, 0;

	return partial_encoding_integers;
}

// Other headers and utility functions

int main(int argc, char** argv);

int autorun_states(int argc, char** argv);

void run_prompts(string& message);

int read_field_file(string filename, FieldSet& fields);

EMatrix gen_encoding_integers();

double envelope_funct(double t);

EMatrix to_full_matrix_hermitian(EMatrix upper);
EMatrix to_full_matrix_antihermitian(EMatrix upper);

pair<pair<EMatrix, EMatrix>, EVector> diag_vec(const EMatrix& mu, const EDMatrix& C);

// OArr evolve_initial_hermitian(const FieldSet& fields, const DipoleSet& dipoles, const EVector& psi_i);
OArr evolve_initial_nonhermitian(const FieldSet& fields, const DipoleSet& dipoles, const EVector& psi_i);

vector<CArray> run_order_analysis(bool prints, const FieldSet& fields, const EVector& psi_i, bool hermitian,
                                  const EMatrix& encoding_integers);

Complex get_only_element(Matrix<Complex, -1, -1> scalar);

pair<int, int> calc_loc(int u_i);

pair<vector<DArr>, EVector> gen_pop_graphs(const FieldSet& fields, const DipoleSet& dipoles, const EVector& psi_i);

void ptime();

template <class T, class F> void print_vec(vector<vector<T>> vec, ofstream& outfile, F lambda);
template <class T, class F, size_t N> void print_arr(vector<array<T, N>> vec, ofstream& outfile, F lambda);

inline double normalize(double rand, double min, double max) {
    return rand * (max - min) + min;
}

// inline EMatrix kroneckerProduct3(const EMatrix2& A, const EMatrix2& B, const EMatrix2& C) {
// 	return kroneckerProduct(A, kroneckerProduct(B, C).eval()).eval();
// }

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