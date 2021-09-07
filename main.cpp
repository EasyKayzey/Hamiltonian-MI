#include "main.h"
#include "omp.h"

double T = 15, DELTA_T, N_T_double = 250;
int N_T;
int ORD = 0;
int BASE = 6;
int main_start_time;

EMatrix V;
// EVector H0D;
// EDMatrix C;
array<ECovector, DIM> anal_pop;

// #define USE_FIELD_FILE
int main(int argc, char** argv) {
    { // this will only work until 2038 so be careful
        time_t now;
        main_start_time = time(&now);
        assert(now == main_start_time);
        ptime();
    }
#ifndef USE_FIELD_FILE
    T = 15;
    N_T = 250;
#endif

    N_T = (int) round(N_T_double);
    DELTA_T = T / N_T;

    string message = "AT1";
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            message += '_';
            message += argv[i];
        }
    }


    V   <<  0, 0.2i, 0,
            0.2i, 0, 0.1i,
            0, 0.1i, 0;

    EVector psi_i = EVector::Zero();
    psi_i[0] = 1;

    for (int i = 0; i < DIM; ++i) {
        ECovector cur = ECovector::Zero();
        cur(i) = 1;
        anal_pop[i] = cur;
    }

#ifdef USE_FIELD_FILE
    FGenome field_genome{};
    // string ffn = string(argv[1]);
    string ffn;
    cout << "Field file name?" << endl;
    cin >> ffn;
    if (ffn.empty() || ffn == "n")
        ffn = "field";

    {
        ifstream field_file(path + ffn + ".txt");
        if (field_file.good()) {
            cout << "Using field from " << ffn << ".txt" << endl;
            try {
                int al;
                field_file >> al;
                int fn = al;
                if (fn == 1) {
                    for (int j = 0; j < L; ++j)
                        field_file >> field_genome[j].first;
                    for (int j = 0; j < L; ++j)
                        field_file >> field_genome[j].second;
                } else if (fn == 0)
                    cout << "Not using fields." << endl;
                else
                    throw runtime_error("Unexpected value of input file fn.");
            } catch (...) {
                cout << "Reading fields failed..." << endl;
                exit(0);
            }
        } else {
            cout << "Reading fields failed..." << endl;
            exit(0);
        }
    }
    vector<double> field = get_field(field_genome);
#else
    cout << "Using hardcoded field..." << endl;
    vector<double> field(N_T);
    fill(field.begin(), field.end(), 1);
    string ffn = "ONES";
#endif

    EMatrix upper_triangle_ones = EMatrix::Zero();
    for (int i = 1; i < DIM; ++i)
        for (int j = 0; j < i; ++j)
            if (V(i, j) != 0.il)
                upper_triangle_ones(i, j) = 1.;
    
    EMatrix encoding_integers;
    enum enc_scheme { other, partial, full };
    enum enc_type { hermitian, antihermitian, nonhermitian };
    const enc_scheme cur_scheme = other;
    const enc_type cur_type = hermitian;

    if (cur_scheme == other) {
        ORD = 256;
        encoding_integers << 0, 1, 0, 0, 0, 8, 0, 0, 0;

    } else if (cur_scheme == partial) {
        encoding_integers << 
                0, 1, 2, 3, 0, 0,
                0, 0, 4, 5, 6, 0,
                0, 0, 0, 7, 8, 0,
                0, 0, 0, 0, 9, 0,
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0;
    } else if (cur_scheme == full) {
        if (cur_type == nonhermitian) {
            encoding_integers = upper_triangle_ones + upper_triangle_ones.transpose();
        } else {
            encoding_integers = upper_triangle_ones;
        }
        int ctr = 0;
        for (Complex& d : upper_triangle_ones.reshaped())
            if (d.real() != 0)
                d = (double) ++ctr;
    } else {
        throw runtime_error("Unsupported encoding scheme.");
    }

    if (cur_type == hermitian) {
        encoding_integers = (-encoding_integers + encoding_integers.transpose()).eval();
    } else if (cur_type == antihermitian) {
        encoding_integers =  (encoding_integers + encoding_integers.transpose()).eval();
    } else if (cur_type == nonhermitian) {
        // do nothing
    } else {
        throw runtime_error("Unsupported encoding type.");
    }

    if (cur_scheme != other) {
        long double max_ord = 0;
        for (Complex& c : encoding_integers.reshaped()) {
            if (c.real() != 0) {
                max_ord = max(max_ord, c.real());
                c = round(pow(BASE, (double) c.real() - 1));
            }
        }
        ORD = 1 << (int) ceil(log2(pow(BASE, (double) max_ord)));
    } else {
        if (ORD == 0)
            throw runtime_error("ORD not set...");
    }

    vector<CArray> anal_res = run_order_analysis(field, psi_i, false, encoding_integers);


    ofstream outfile(string(path) + "HMI_" + to_string(main_start_time) + "_" + ffn + (message == "#" ? "" : "_" + message) + ".txt");

    int out_ints[] = {DIM, N_T, main_start_time, L, N_H, N_TO, N_OBS, ORD, BASE};
    double out_doubles[] = {T, HBAR};

    if (message.empty())
        message = "#";
    outfile << "HMI2 " << 0 << ' ' << message << ' ' << time(nullptr) - main_start_time << endl;
    for (int o : out_ints)
        outfile << ' ' << o;
    outfile << endl;
    for (double o : out_doubles)
        outfile << ' ' << o;
    outfile << endl;

    outfile << "Preliminaries:" << endl;

    // outfile << H0D.real().transpose() << endl;
    outfile << -1 << endl;
    outfile << V.imag() << endl;

    // outfile << mu_t_upper.real() + mu_t_upper.real().transpose() << endl;

    outfile << psi_i.real().transpose() << endl;

    outfile << "Field:" << endl;

    for (double d : field)
        outfile << d << ' ';
    outfile << endl;

    outfile << encoding_integers.real() << endl;

    for (CArray& arr : anal_res) {
        for (Complex d : arr)
            outfile << d.real() << ' ';
        outfile << endl;
        for (Complex d : arr)
            outfile << d.imag() << ' ';
        outfile << endl;
    }   

    if (!outfile.good())
        cerr << "Writing failed." << endl;

    outfile.close();

}


double envelope_funct(double t) {
    // static_assert(N_TO == 2, "The current envelope function is a double bell curve...\n");
    // return exp(-30 * (2 * t / T - .5) * (2 * t / T - .5)) + exp(-30 * ((2 * t - T) / T - .5) * ((2 * t - T) / T - .5));
    static_assert(N_TO == 1, "The current envelope function is a single bell curve...\n");
    return exp(-30 * (t / T - .5) * (t / T - .5));
}

EMatrix to_full_matrix_hermitian(EMatrix upper) {
    return (upper + upper.adjoint()).transpose(); // transpose is because upper triangle is backwards transitions
}

EMatrix to_full_matrix_antihermitian(EMatrix upper) {
    return upper + upper.transpose();
}

pair<pair<EMatrix, EMatrix>, EVector> diag_vec(const EMatrix& mu, const EDMatrix& C) {
    ComplexSchur<EMatrix> schur;
    schur.compute(mu, true);
    if (schur.info() != Success) {
        cout << "Schur computation failed." << endl;
        exit(1);
    }

//    cout << "C: " << endl << C.toDenseMatrix() << endl;

    EVector lambda;
    lambda = schur.matrixT().diagonal();
    EMatrix P;
    P = schur.matrixU();
//    cout << "Lambda diagonal: " << endl << lambda.transpose() << endl;
//    cout << "P: " << endl << P << endl;
//    cout << "P Lambda P^d:" << endl << P * lambda.asDiagonal() * P.adjoint() << endl;

    EMatrix CP, PdC;
    CP = C * P;
    PdC = P.adjoint() * C;
    return {{CP, PdC}, lambda};
}

// OArr evolve_initial_hermitian(const vector<double>& epsilon, const EMatrix& mu, const EVector& psi_i) {
//     auto diag_ret = diag_vec(mu, C);
//     EVector lambda = diag_ret.second;
//     EMatrix CP = diag_ret.first.first;
//     EMatrix PdC = diag_ret.first.second;
//     EMatrix PdCCP = PdC * CP;

//     vector<EDMatrix, aligned_allocator<EDMatrix>> E(N_T);
//     for (int i = 0; i < N_T; ++i)
//         E[i] = exp(1i * DELTA_T / HBAR * lambda.array() * epsilon[i]).matrix().asDiagonal();

//     vector<EVector, aligned_allocator<EVector>> it(N_T + 1);
//     it[0] = PdC * psi_i;
//     for (int i = 1; i < N_T; ++i)
//         it[i] = PdCCP * (E[i - 1] * it[i - 1]);
//     it[N_T] = CP * (E[N_T - 1] * it[N_T - 1]);

//     OArr samples{};
//     for (int i = 0; i < N_TO; ++i)
//         for (int j = 0; j < DIM; ++j)
//             samples[i * DIM + j] = it[(i + 1) * N_T / N_TO][j];
//     return samples;
// }

OArr evolve_initial_nonhermitian(const vector<double>& epsilon, const EMatrix& mu, const EVector& psi_i) {
    // auto diag_ret = diag_vec(mu, C);
    // EVector lambda = diag_ret.second;
    // EMatrix CP = diag_ret.first.first;
    // EMatrix PdC = diag_ret.first.second;
    // EMatrix PdCCP = PdC * CP;

    vector<EMatrix, aligned_allocator<EMatrix>> Ex(N_T);
    for (int i = 0; i < N_T; ++i)
        Ex[i] = (1i * DELTA_T / HBAR * mu * epsilon[i]).exp();

    vector<EVector, aligned_allocator<EVector>> it(N_T + 1);
    it[0] = psi_i;
    for (int i = 1; i < N_T; ++i)
        it[i] = (Ex[i - 1] * it[i - 1]);
    it[N_T] =  (Ex[N_T - 1] * it[N_T - 1]);

    OArr samples{};
    for (int i = 0; i < N_TO; ++i)
        for (int j = 0; j < DIM; ++j)
            samples[i * DIM + j] = it[(i + 1) * N_T / N_TO][j];
    return samples;
}

vector<CArray> run_order_analysis(const vector<double>& epsilon, const EVector& psi_i, bool hermitian, const EMatrix& encoding_integers) {
    cout << time(nullptr) << endl;
    int ord = ORD;
    vector<OArr> order_results(ord);
    cout << "Running analysis with ord=" << ord << endl;
#pragma omp parallel for default(shared)
    for (int s = 0; s < ord; ++s) {
        if (s % 1000 == 0)
            cout << "Doing s=" << s << endl;
        EMatrix mu = V;
        double g = 2 * M_PI * s / ord;
        EMatrix encoding = encoding_integers;
        for (Complex& d : encoding.reshaped())
            d = polar(1.l, d.real() * g);
        EMatrix encoded = (mu.array() * encoding.array());
        if (hermitian)
            order_results[s] = evolve_initial_nonhermitian(epsilon, encoded, psi_i);
        else
            order_results[s] = evolve_initial_nonhermitian(epsilon, encoded, psi_i);
    }
    cout << time(nullptr) << endl;

    vector<CArray> ffts;
    for (int i = 0; i < DIM; ++i) {
        int ii = i + DIM * (N_TO - 1);
        CArray tfft(ord);
        for (int j = 0; j < ord; ++j) {
            tfft[j] = order_results[j][ii];
        }
        fft(tfft);
        tfft /= tfft.size();
        cout << "\nFFT for 1 to " << i + 1 << ':' << endl;
        // for (auto& d : tfft)
        //     cout << abs(d) << ' ';
        for (int k = 0; k < ord; ++k) {
            if (abs(tfft[k]) > 0.01)
                cout << k << ": " << abs(tfft[k]) << endl;
        }
        cout << endl << "Sum of values: " << tfft.sum() << "; magnitude " << abs(tfft.sum()) << "; prob " << norm(tfft.sum()) << endl;
        ffts.push_back(tfft);
    }
    double sum_of_probs = 0;
    for (CArray& fft : ffts)
        sum_of_probs += norm(fft.sum());
    cout << "\nSum of calculated sum-of-fftval probablities is " << sum_of_probs << endl;

    return ffts;
}


complex<double> get_only_element(Matrix<complex<double>, -1, -1> scalar) {
    if (scalar.rows() > 1 || scalar.cols() > 1) {
        cout << scalar << endl;
        throw runtime_error("Tried to get single element from matrix, see cout for matrix");
    }
    return scalar(0, 0);
}

pair<int, int> calc_loc(int u_i) { // I could binary search this but I'm too lazy
    for (int i = DIM - 1; i > 0; --i)
        if ((i * (i - 1)) / 2 <= u_i)
            return {i, u_i - ((i * (i - 1)) / 2)};
    throw runtime_error("calc_loc failed");
}

/*
vector<DArr> gen_pop_graphs(const vector<double>& eps_inter, const EMatrix& CP, const EMatrix& PdC,
                            const EVector& lambda, const EVector& psi_i, const array<ECovector, DIM>& anal_pop) {
    vector<EMatrix, aligned_allocator<EMatrix>> E(N_T);
    for (int i = 0; i < N_T; ++i)
        E[i] = exp(1i * DELTA_T / HBAR * lambda.array() * eps_inter[i]).matrix().asDiagonal();

    vector<EVector, aligned_allocator<EVector>> it(N_T + 1);
    it[0] = psi_i;
    for (int i = 1; i <= N_T; ++i)
        it[i] = CP * (E[i - 1] * (PdC * it[i - 1]));

    vector<DArr> pops(N_T + 1);
    for (int i = 0; i <= N_T; ++i)
        for (int o = 0; o < DIM; ++o)
            pops[i][o] = (anal_pop[o] * it[i]).squaredNorm();
    return pops;
}
*/

// vector<double> get_field(const FGenome& epsilon) {
//     // throw runtime_error("Not implemented");
//     array<double, L> omega;

//     for (int l = 1, c = 0; l < DIM; ++l)
//         for (int m = 0; m < l; ++m)
//             omega[c++] = (H0D(l) - H0D(m)).real() / HBAR;

//     vector<double> eps_inter(N_T);
//     for (int i = 0; i < N_T; ++i) {
//         eps_inter[i] = 0;
//         double t_i = (i + 0.5) * DELTA_T;
//         int c = 0;
//         for (int l = 1; l < DIM; ++l) {
//             for (int m = 0; m < l; ++m) {
//                 eps_inter[i] += epsilon[c].first * sin(omega[c] * t_i + epsilon[c].second);
//                 ++c;
//             }
//         }
//         eps_inter[i] *= envelope_funct(t_i);
//     }
//     return eps_inter;
// }

template <class T, class F> void print_vec(vector<vector<T>> vec, ofstream& outfile, F lambda) {
    for (vector<T> i : vec) {
        for (T j : i)
            outfile << lambda(j) << ' ';
        outfile << endl;
    }
}

template <class T, class F, size_t N> void print_arr(vector<array<T, N>> vec, ofstream& outfile, F lambda) {
    for (array<T, N> i : vec) {
        for (T j : i)
            outfile << lambda(j) << ' ';
        outfile << endl;
    }
}

void ptime() {
    time_t t = time(nullptr);
    cout << "Unix time " << (long) t << ", runtime " << (long) (t - main_start_time) << ", date " << ctime(&t);
}
