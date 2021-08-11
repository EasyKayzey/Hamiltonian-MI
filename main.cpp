#include "main.h"

double T = 4000, DELTA_T, N_T_double = 250;
int N_T;
int ORD_R = 128;
EMatrix mu_t_upper;
EVector H0D;
EDMatrix C;
int main_start_time;

int main(int argc, char** argv) {
    { // this will only work until 2038 so be careful
        time_t now;
        main_start_time = time(&now);
        assert(now == main_start_time);
        ptime();
    }
    N_T = (int) round(N_T_double);
    DELTA_T = T / N_T;

    string message = "AT1";
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            message += '_';
            message += argv[i];
        }
    }

    H0D << 0, 0.00820226918, 0.01558608386, 0.02215139847, 0.02789825857, 0.03282661861;
    C = exp(H0D.array() * -1i * DELTA_T / 2 / HBAR).matrix().asDiagonal();

    mu_t_upper   <<  0, 0.06116130402, -0.01272999623, 0.003275382148,  0,               0,
                     0, 0,              0.0834968862, -0.02150346764,   0.006457586337,  0,
                     0, 0,              0,             0.09850502453,  -0.02960110398,   0.01003668145,
                     0, 0,              0,             0,               0.1093531798,   -0.0371225703,
                     0, 0,              0,             0,               0,               0.1171966409,
                     0, 0,              0,             0,               0,               0;

    EVector psi_i = EVector::Zero();
    psi_i[0] = 1;

    array<ECovector, DIM> anal_pop;
    for (int i = 0; i < DIM; ++i) {
        ECovector cur = ECovector::Zero();
        cur(i) = 1;
        anal_pop[i] = cur;
    }

    FGenome field_genome{};
    {
        ifstream field_file(path + string(argv[1]) + ".txt");
        if (field_file.good()) {
            cout << "Using field from " << argv[1] << ".txt" << endl;
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

    vector<OArr> order_results(ORD_R);
#pragma omp parallel for default(shared)
    for (int s = 0; s < ORD_R; ++s) {
        EMatrix mu_upper = mu_t_upper;
        double g = 2 * M_PI * s / ORD_R;
        Complex m = polar(1., g);
        mu_upper *= m;
        order_results[s] = evolve_initial(field, to_full_matrix(mu_upper), psi_i, anal_pop);
    }


}


double envelope_funct(double t) {
    static_assert(N_TO == 2, "The current envelope function is a double bell curve...\n");
    return exp(-30 * (2 * t / T - .5) * (2 * t / T - .5)) + exp(-30 * ((2 * t - T) / T - .5) * ((2 * t - T) / T - .5));
}

EMatrix to_full_matrix(EMatrix upper) {
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

OArr evolve_initial(const vector<double>& epsilon, const EMatrix& mu, const EVector& psi_i, const array<ECovector, DIM>& anal_pop) {
    auto diag_ret = diag_vec(mu, C);
    EVector lambda = diag_ret.second;       
    EMatrix CP = diag_ret.first.first;
    EMatrix PdC = diag_ret.first.second;
    EMatrix PdCCP = PdC * CP;

    vector<EDMatrix, aligned_allocator<EDMatrix>> E(N_T);
    for (int i = 0; i < N_T; ++i)
        E[i] = exp(1i * DELTA_T / HBAR * lambda.array() * epsilon[i]).matrix().asDiagonal();

    vector<EVector, aligned_allocator<EVector>> it(N_T + 1);
    it[0] = PdC * psi_i;
    for (int i = 1; i < N_T; ++i)
        it[i] = PdCCP * (E[i - 1] * it[i - 1]);
    it[N_T] = CP * (E[N_T - 1] * it[N_T - 1]);

    OArr samples{};
    for (int i = 0; i < N_TO; ++i)
        for (int j = 0; j < DIM; ++j)
            samples[i * DIM + j] = (anal_pop[j] * it[(i + 1) * N_T / N_TO]).squaredNorm();
    return samples;
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

vector<double> get_field(const FGenome& epsilon) {
    throw runtime_error("Not implemented");
    // vector<double> eps_inter(N_T);
    // for (int i = 0; i < N_T; ++i) {
    //     eps_inter[i] = 0;
    //     double t_i = (i + 0.5) * DELTA_T;
    //     int c = 0;
    //     for (int l = 1; l < DIM; ++l) {
    //         for (int m = 0; m < l; ++m) {
    //             eps_inter[i] += epsilon[c].first * sin(omega[c] * t_i + epsilon[c].second);
    //             ++c;
    //         }
    //     }
    //     eps_inter[i] *= envelope_funct(t_i);
    // }
    // return eps_inter;
}

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
