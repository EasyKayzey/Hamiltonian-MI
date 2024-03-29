#include "main.h"
#include "omp.h"

double T;
int N_T;
int ORD = 0;
int time_scale = 1;
int main_start_time;
double amp_scale = 1;
bool ask_any_prompts = true;

bool use_t_arr = false;
bool ask_t_scale = false;
bool ask_amp_scale = false;
bool print_ffts = false;
bool rerun = true;

int autostate = -1;
string autofield = "";
string automessage = "";

time_t now = 0;
array<ECovector, DIM> anal_pop;

int main(int argc, char** argv) {
    if (now == 0) {
        main_start_time = time(&now);  // this will only work until 2038 so be careful
        assert(now == main_start_time);
        ptime();

        gen_hamiltonians();

        cout << "5 initial prompts: tarr, tscale, ascale, printffts, allstates" << endl;

        char c_tmp;
        cout << "If T arrays are desired, enter \"y\" at the prompt." << endl;
        cin >> c_tmp;
        if (c_tmp == 'Y' || c_tmp == 'y') {
            use_t_arr = true;
            cout << "Using T arrays." << endl;
        } else {
            use_t_arr = false;
        }
        
        cout << "If you'd like to be prompted for time point scale, enter \"y\" at the prompt." << endl;
        cin >> c_tmp;
        if (c_tmp == 'Y' || c_tmp == 'y') {
            ask_t_scale = true;
        } else {
            ask_t_scale = false;
        }
                
        cout << "If you'd like to be prompted for amplitude scale, enter \"y\" at the prompt." << endl;
        cin >> c_tmp;
        if (c_tmp == 'Y' || c_tmp == 'y') {
            ask_amp_scale = true;
        } else {
            ask_amp_scale = false;
        }

        cout << "If you'd like to print FFTs, enter \"y\" at the prompt." << endl;
        cin >> c_tmp;
        if (c_tmp == 'Y' || c_tmp == 'y') {
            print_ffts = true;
        } else {
            print_ffts = false;
        }

        cout << "If you'd like to autorun all initial states, enter \"y\" at the prompt." << endl;
        cin >> c_tmp;
        if (c_tmp == 'Y' || c_tmp == 'y') {
            return autorun_states(argc, argv);
        }
    }

    string message = "";
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            message += '_';
            message += argv[i];
        }
    }
    string message_backup = message;

    FieldSet fields{};
    // string ffn = string(argv[1]);
    string ffn = autofield;
    if (autofield.empty()) {
        cout << "Reading fields from " << field_dir << ". Field file name?" << endl;
        cin >> ffn;
        if (ffn.empty() || ffn == "#")
            throw runtime_error("Need field file name!");
        else if (ffn == "exit")
            return 0;
    }

    ffn = field_dir + ffn;

    string init_state_str = "";
    int init_state;
    if (autostate == -1) {
        cout << "Initial state index?" << endl;
        cin >> init_state_str;
        if (!init_state_str.empty() && init_state_str != "0") {
            init_state = stoi(init_state_str);
        } else {
            init_state = 0;
        }
        message += (message.length() == 0 ? "S" : "_S") + to_string(init_state);
    } else {
        init_state = autostate;
        message += (message.length() == 0 ? "AS" : "_AS") + to_string(init_state);
    }

    if (ask_any_prompts)
        run_prompts(message);

    string message_append = automessage;
    if (message_append.empty()) {
        cout << "Message append? (can use # for no)" << endl;
        cin >> message_append;
    }
    if (message_append.length() != 0 && message_append != "#")
        message += (message.length() == 0 ? "" : "_") + message_append;

    EVector psi_i = EVector::Zero();
    psi_i[init_state] = 1;
    // EVector psi_i;
    // psi_i << 2.25974e-05+1.39985e-05i, 0.47164-0.881791i, -5.03996e-06-1.80181e-05i;

    for (int i = 0; i < DIM; ++i) {
        ECovector cur = ECovector::Zero();
        cur(i) = 1;
        anal_pop[i] = cur;
    }


    if (read_field_file(path + ffn + ".txt", fields) == -1) {
        cout << "Field file not found. Exiting..." << endl;
        exit(-1);
    }

    cout << "Successfully read fields." << endl;

    EMatrix encoding_integers = gen_encoding_integers();

    cout << "Running end analysis..." << endl;
    vector<CArray> anal_res_end = run_order_analysis(print_ffts, fields, psi_i, cur_type == hermitian, encoding_integers);

    DipoleSet dipoles = dipoles_upper;
    for (EMatrix &dipole : dipoles) dipole = (dipole + dipole.adjoint()).eval(); 
    auto PGR = gen_pop_graphs(fields, dipoles, psi_i);


    ofstream outfile(string(path) + "HMI_" + to_string(main_start_time) + "_" + ffn.substr(ffn.find_last_of("/\\") + 1) 
                     + (message == "#" || message.empty() ? "" : "_" + message) + ".txt");

    int out_ints[] = {DIM, N_T, main_start_time, L, N_H, N_TO, time_scale, N_FIELDS, ORD, BASE, 
                      cur_scheme, cur_type, (use_t_arr ? N_T : 1), 0,0,0,0};
    double out_doubles[] = {T, HBAR, amp_scale, 0,0,0,0};

    if (message.empty())
        message = "#";
    outfile << "HMR2 " << 1 << ' ' << message << ' ' << time(nullptr) - main_start_time << endl;
    for (int o : out_ints)
        outfile << ' ' << o;
    outfile << endl;
    for (double o : out_doubles)
        outfile << ' ' << o;
    outfile << endl;

    outfile << "Preliminaries:" << endl;

    outfile << H0D.real().transpose() << endl;

    for (EMatrix& m : dipoles_upper)
        outfile << m + m.adjoint() << endl;

    outfile << psi_i.real().transpose() << endl;
    outfile << psi_i.imag().transpose() << endl;

    outfile << PGR.second.real().transpose() << endl;
    outfile << PGR.second.imag().transpose() << endl;

    for (int i = 0; i < DIM; ++i) {
        for (DArr &a : PGR.first) {
            outfile << a[i] << ' ';
        }
        outfile << endl;
    }

    outfile << "Fields:" << endl;

    for (vector<double>& field : fields) {
        for (double d : field)
            outfile << d << ' ';
        outfile << endl;
    }

    outfile << encoding_integers.real() << endl;

    if (use_t_arr) {
        cout << "Running time analysis on " << N_T << " points..." << endl;
        int N_T_back = N_T;
        for (int i = 0; i < N_T_back; ++i) {
            outfile << i << endl;
            if (N_T_back < 100 || i % (N_T_back / 100) == 0)
                cout << "Running time point " << i << endl;
            N_T = i + 1;
            auto anal_res = run_order_analysis(false, fields, psi_i, cur_type == hermitian, encoding_integers);
            for (CArray& arr : anal_res) {
                for (Complex d : arr)
                    outfile << d.real() << ' ';
                outfile << endl;
                for (Complex d : arr)
                    outfile << d.imag() << ' ';
                outfile << endl;
            }   
        }
        cout << "Finished time analysis." << endl;
        N_T = N_T_back;
    } else {
        outfile << N_T << endl;
        for (CArray& arr : anal_res_end) {
            for (Complex d : arr)
                outfile << d.real() << ' ';
            outfile << endl;
            for (Complex d : arr)
                outfile << d.imag() << ' ';
            outfile << endl;
        }   
    }

    if (!outfile.good())
        cerr << "Writing failed." << endl;
    else {
        outfile.close();
        cout << "Finished writing." << endl;
    }

    message = message_backup;
    if (rerun)
        return main(argc, argv);
    return 0;
}

int autorun_states(int argc, char** argv) {
    bool rerunback = rerun, aapback = ask_any_prompts;
    rerun = false;
    ask_any_prompts = false;

    cout << "Reading fields from " << field_dir << ". Field file name?" << endl;
    cin >> autofield;
    if (autofield.empty() || autofield == "#")
        throw runtime_error("Need field file name!");
    else if (autofield == "exit")
        return 0;

    automessage = "";
    run_prompts(automessage);
    
    string automessage_append = "";
    cout << "Message append? (can use # for no)" << endl;
    cin >> automessage_append;
    if (automessage_append.length() != 0 && automessage_append != "#")
        automessage += (automessage.length() == 0 ? "" : "_") + automessage_append;
    if (automessage.empty())
        automessage = "#";

    for (int i = 0; i < DIM; ++i) {
        autostate = i;
        cout << "Running initial state " << i << "..." << endl;
        main(argc, argv);
    }

    autofield = "";
    autostate = -1;
    rerun = rerunback;
    ask_any_prompts = aapback;
    if (rerun)
        return autorun_states(argc, argv);
    return 0;
}

void run_prompts(string& message) {
    if (ask_t_scale) {
        string t_scale_str = "";
        cout << "Time point scale?" << endl;
        cin >> t_scale_str;
        if (!t_scale_str.empty() && t_scale_str != "1") {
            time_scale = stoi(t_scale_str);
            message += (message.length() == 0 ? "" : "_") + t_scale_str;
        } else {
            time_scale = 1;
        }
    }
    
    if (ask_amp_scale) {
        string amp_scale_str = "";
        cout << "Amplitude scale?" << endl;
        cin >> amp_scale_str;
        if (!amp_scale_str.empty() && amp_scale_str != "1") {
            amp_scale = stod(amp_scale_str);
            message += (message.length() == 0 ? "" : "_") + amp_scale_str;
        } else {
            amp_scale = 1;
        }
    }
}

int read_field_file(string filename, FieldSet& fields) {
    ifstream field_file(filename);
    if (field_file.good()) {
        cout << "Using field from " << filename << endl;
        try {
            double d;
            int n_fields, n_skip;
            string s_tmp;

            getline(field_file, s_tmp);
            if (s_tmp[s_tmp.length() - 1] == '\r')
                s_tmp = s_tmp.substr(0, s_tmp.length() - 1);
            if (s_tmp != "TIME N_T N_FIELDS SKIP")
                throw runtime_error("Field file header " + s_tmp + " incorrect.");
            
            field_file >> T;
            field_file >> N_T;
            field_file >> n_fields;
            if (n_fields != N_FIELDS)
                throw runtime_error("Field file N_FIELDS=" + to_string(n_fields)
                                        + " does not match header " + to_string(N_FIELDS) + ".");
            field_file >> n_skip;
            cout << "Read T=" << T << ", N_T=" << N_T << ", SKIP=" << n_skip << endl;
            if (time_scale != 1) {
                N_T = N_T * time_scale;
                cout << "With time scale " << time_scale << " applied, N_T=" << N_T << endl;
            }

            getline(field_file, s_tmp); // needed to get to end of line 2
            getline(field_file, s_tmp);
            if (s_tmp[s_tmp.length() - 1] == '\r')
                s_tmp = s_tmp.substr(0, s_tmp.length() - 1);
            if (s_tmp != "FIELDS")
                throw runtime_error("Field file line 3 " + s_tmp + " incorrect.");

            fill(fields.begin(), fields.end(), vector<double>(N_T));
            for (int i = 0; i < N_T / time_scale; ++i) {
                for (int k = 0; k < n_skip; ++k)
                    field_file >> d;
                for (int j = 0; j < N_FIELDS; ++j) {
                    field_file >> d;
                    for (int k = 0; k < time_scale; ++k) {
                        fields[j][i * time_scale + k] = d * amp_scale;
                    }
                }
            }
            if (!field_file.eof()) {
                char c;
                while (field_file.get(c))
                    if (!std::isspace(c))
                        throw runtime_error("Field file too long...");
            }
            return 0;
        } catch (runtime_error& e) {
            cout << "Reading fields failed... Error: " << e.what() << endl;
            cerr << "Reading fields failed... Error: " << e.what() << endl;
            exit(-1);
        }
    } else {
        return -1;
    }
}

EMatrix gen_encoding_integers() {
	EMatrix upper_triangle_ones = EMatrix::Zero();
    for (int i = 1; i < DIM; ++i)
        for (int j = 0; j < i; ++j)
            for (EMatrix& upper : dipoles_upper)
                if (upper(i, j) != upper(j, i))
                    upper_triangle_ones(i, j) = 1.;

    EMatrix encoding_integers;

    if (cur_scheme == other) {

    } else if (cur_scheme == order) {
        encoding_integers = upper_triangle_ones + (cur_type == nonhermitian) * upper_triangle_ones.transpose();
    } else if (cur_scheme == partial) {
        encoding_integers = get_partial_encoding_integers();
    } else if (cur_scheme == full) {
        if (cur_type == nonhermitian) {
            encoding_integers = upper_triangle_ones + upper_triangle_ones.adjoint();
        } else {
            encoding_integers = upper_triangle_ones;
        }
        int ctr = 0;

        for (Complex& d : encoding_integers.reshaped())
            if (d.real() != 0)
                d = (double) ++ctr;
    } else {
        throw runtime_error("Unsupported encoding scheme.");
    }

    if (cur_type == hermitian) {
        encoding_integers = (encoding_integers - encoding_integers.adjoint()).eval();
    } else if (cur_type == antihermitian) {
        encoding_integers = (encoding_integers + encoding_integers.adjoint()).eval();
    } else if (cur_type == nonhermitian) {
        // do nothing
    } else {
        throw runtime_error("Unsupported encoding type.");
    }

    if (cur_scheme != other) {
        DOUBLE_TYPE max_ord = 0;
        for (Complex& c : encoding_integers.reshaped()) {
            if (c.real() != 0) {
                max_ord = max(max_ord, abs(c.real()));
                c = copysign(round(pow(BASE, (double) abs(c.real()) - 1)), c.real());
            }
        }
        ORD = 1 << (int) ceil(log2(pow(BASE, (double) max_ord)));
    } else {
        if (ORD == 0)
            throw runtime_error("ORD not set...");
    }

    return encoding_integers;
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
        exit(-1);
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

// OArr evolve_initial_hermitian(const FieldSet& fields, const DipoleSet& dipoles, const EVector& psi_i) {
//     throw runtime_error("Not implemented right now");
//     // auto diag_ret = diag_vec(mu, C);
//     // EVector lambda = diag_ret.second;
//     // EMatrix CP = diag_ret.first.first;
//     // EMatrix PdC = diag_ret.first.second;
//     // EMatrix PdCCP = PdC * CP;

//     // vector<EDMatrix, aligned_allocator<EDMatrix>> E(N_T);
//     // for (int i = 0; i < N_T; ++i)
//     //     E[i] = exp(1i * T / N_T / HBAR * lambda.array() * epsilon[i]).matrix().asDiagonal();

//     // vector<EVector, aligned_allocator<EVector>> it(N_T + 1);
//     // it[0] = PdC * psi_i;
//     // for (int i = 1; i < N_T; ++i)
//     //     it[i] = PdCCP * (E[i - 1] * it[i - 1]);
//     // it[N_T] = CP * (E[N_T - 1] * it[N_T - 1]);

//     // OArr samples{};
//     // for (int i = 0; i < N_TO; ++i)
//     //     for (int j = 0; j < DIM; ++j)
//     //         samples[i * DIM + j] = it[(i + 1) * N_T / N_TO][j];
//     // return samples;
// }

OArr evolve_initial_nonhermitian(const FieldSet& fields, const DipoleSet& dipoles, const EVector& psi_i) {

    vector<EMatrix, aligned_allocator<EMatrix>> Hs(N_T);
    for (int i = 0; i < N_T; ++i) {
        Hs[i] = H0D.asDiagonal();
        for (int j = 0; j < N_FIELDS; ++j)
            Hs[i] += fields[j][i] * dipoles[j];
    }

    vector<EVector, aligned_allocator<EVector>> it(N_T + 1);
    it[0] = psi_i;
    for (int i = 1; i <= N_T; ++i)
        it[i] = (-1.i * Hs[i-1] * T / N_T).exp() * it[i - 1];

    OArr samples{};
    for (int i = 0; i < N_TO; ++i)
        for (int j = 0; j < DIM; ++j)
            samples[i * DIM + j] = it[(i + 1) * N_T / N_TO][j];
    return samples;
}

vector<CArray> run_order_analysis(bool prints, const FieldSet& fields, const EVector& psi_i, bool hermitian, const EMatrix& encoding_integers) {
    // cout << time(nullptr) << endl;
    int ord = ORD;
    vector<OArr> order_results(ord);
    if (prints)
        cout << "Running analysis with ord=" << ord << endl;

#pragma omp parallel for default(shared)
    for (int s = 0; s < ord; ++s) {
        if (prints && s % 1000 == 0)
            cout << "Doing s=" << s << endl;
        double g = 2 * MY_PI * s / ord;
        EMatrix encoding = encoding_integers;
        for (Complex& d : encoding.reshaped())
            d = polar((DOUBLE_TYPE) 1, d.real() * g);
        DipoleSet encoded;
        for (int i = 0; i < N_FIELDS; ++i)
            encoded[i] = ((dipoles_upper[i] + dipoles_upper[i].adjoint()).array() * encoding.array()).matrix();
        // if (hermitian)
        //     order_results[s] = evolve_initial_hermitian(fields, encoded, psi_i);
        // else
        order_results[s] = evolve_initial_nonhermitian(fields, encoded, psi_i);
    }
    if (hermitian && !hermitian)
        cout << "This print exists to remove an unused variable warning.\n";
    // cout << time(nullptr) << endl;

    vector<CArray> ffts;
    for (int i = 0; i < DIM; ++i) {
        int ii = i + DIM * (N_TO - 1);
        CArray tfft(ord);
        for (int j = 0; j < ord; ++j) {
            tfft[j] = order_results[j][ii];
        }
        fft(tfft);
        tfft /= tfft.size();
        
        if (prints) {
            cout << "\nFFT for 1 to " << i + 1 << ":\n";
            // for (auto& d : tfft)
            //     cout << abs(d) << ' ';
            for (int k = 0; k < ord; ++k) {
                if (abs(tfft[k]) > 0.01)
                    cout << k << ": " << abs(tfft[k]) << '\n';
            }
            cout << endl << "Sum of values: " << tfft.sum() << "; magnitude " << abs(tfft.sum()) << "; prob " << norm(tfft.sum()) << '\n';
        }

        ffts.push_back(tfft);
    }
    double sum_of_probs = 0;
    for (CArray& fft : ffts)
        sum_of_probs += norm(fft.sum());

    if (prints) {
        cout << "\nFinal state when unmodulated: ";
        for (const Complex& c : order_results[0])
            cout << c << ' ';
        // auto default_prec = cout.precision();
        // cout << std::setprecision(8);
        cout << "\nSum of calculated sum-of-fftval probablities is " << sum_of_probs;
        cout << "; error is " << sum_of_probs - 1 << '\n';
        // cout << std::setprecision(default_prec);
    }

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

pair<vector<DArr>, EVector> gen_pop_graphs(const FieldSet& fields, const DipoleSet& dipoles, const EVector& psi_i) {

    vector<EMatrix, aligned_allocator<EMatrix>> Hs(N_T);
    for (int i = 0; i < N_T; ++i) {
        Hs[i] = H0D.asDiagonal();
        for (int j = 0; j < N_FIELDS; ++j)
            Hs[i] += fields[j][i] * dipoles[j];
    }

    vector<EVector, aligned_allocator<EVector>> it(N_T + 1);
    it[0] = psi_i;
    for (int i = 1; i <= N_T; ++i)
        it[i] = (-1.i * Hs[i-1] * T / N_T).exp() * it[i - 1];

    vector<DArr> samples;
    for (EVector& iv : it) {
        DArr a;
        for (int i = 0; i < DIM; ++i)
            a[i] = norm(iv(i));
        samples.push_back(a);
    }
    return {samples, it[N_T]};
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
