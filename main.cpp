#include "main.h"
#include "omp.h"

double T = HEADER_TIME, DELTA_T, N_T_double = HEADER_NTD * TIME_POINT_SCALE;
int N_T;
int ORD = 0;
int BASE = 7;
int main_start_time;
double field_scale_factor = 1;

DipoleSet dipoles_upper;
EVector H0D;
// EDMatrix C;
array<ECovector, DIM> anal_pop;

#define USE_FIELD_FILE
#define USE_GOTO
int main(int argc, char** argv) {
    { // this will only work until 2038 so be careful
        time_t now;
        main_start_time = time(&now);
        assert(now == main_start_time);
        ptime();
    }
#ifndef USE_FIELD_FILE
    T = 20000;
    N_T = 1200;
#endif

    N_T = (int) round(N_T_double);
    DELTA_T = T / N_T;

    string message = "";
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            message += '_';
            message += argv[i];
        }
    }
    string message_backup = message;

    H0D << 0,   0.0082,  .016;

    dipoles_upper[0] <<  0,   0.061,  -.013,
                       0,   0,  .083,
                       0,  0,  0;


    start_label:

#ifdef USE_FIELD_FILE
    FieldSet fields{};
    fill(fields.begin(), fields.end(), vector<double>(N_T));
    // string ffn = string(argv[1]);
    string ffn;
    cout << "Field file name?" << endl;
    cin >> ffn;
    if (ffn.empty() || ffn == "n")
        ffn = "field_nmr";

    ffn = "fields/gaurav_algo_fields/" + ffn;

    string amul;
    // cout << "Amplitude multiplier?" << endl;
    // cin >> amul;
    if (amul.empty() || amul == "1")
        field_scale_factor = 1;
    else {
        message += (message.length() == 0 ? "" : "_") + amul;
        field_scale_factor = stod(amul);
    }

    string init_state_str;
    int init_state;
    // cout << "Initial state index?" << endl;
    // cin >> init_state_str;
    if (!init_state_str.empty() && init_state_str != "0") {
        init_state = stoi(init_state_str);
        message += (message.length() == 0 ? "" : "_") + init_state_str;
    } else {
        init_state = 0;
    }

    string t_str;
    cout << "T?" << endl;
    cin >> t_str;
    if (!t_str.empty() && t_str != "0") {
        T = stod(t_str);
        message += (message.length() == 0 ? "" : "_") + t_str;
        DELTA_T = T / N_T;
    }

    string n_t_str;
    cout << "N_T?" << endl;
    cin >> n_t_str;
    if (!n_t_str.empty() && n_t_str != "0") {
        // cout << 'a';
        N_T = stoi(n_t_str);
        message += (message.length() == 0 ? "" : "_") + n_t_str;
        DELTA_T = T / N_T;
    }
        

    string message_append;
    cout << "Message append? (can use # for no)" << endl;
    cin >> message_append;
    if (message_append.length() != 0 && message_append != "#")
        message += (message.length() == 0 ? "" : "_") + message_append;


    EVector psi_i = EVector::Zero();
    psi_i[init_state] = 1;

    for (int i = 0; i < DIM; ++i) {
        ECovector cur = ECovector::Zero();
        cur(i) = 1;
        anal_pop[i] = cur;
    }


    {
        ifstream field_file(path + ffn + ".txt");
        if (field_file.good()) {
            cout << "Using field from " << ffn << ".txt" << endl;
            try {
                double d;
                for (int i = 0; i < N_T / TIME_POINT_SCALE; ++i) {
                    for (int j = 0; j < N_FIELDS; ++j) {
                        field_file >> d; // doubled up because of the stupid time input
                        field_file >> d;
                        for (int k = 0; k < TIME_POINT_SCALE; ++k) {
                            fields[j][i * TIME_POINT_SCALE + k] = d;
                        }
                    }
                }
                if (!field_file.eof()) {
                    char c;
                    while (field_file.get(c))
                        if (!std::isspace(c))
                            throw runtime_error("Field file too long...");
                }
            } catch (runtime_error& e) {
                cout << "Reading fields failed... Error: " << e.what() << endl;
                exit(0);
            }
        } else {
            cout << "Reading fields failed..." << endl;
#ifdef USE_GOTO
            message = message_backup;
            goto start_label;
#else
            exit(0);
#endif
        }
    }
#else
    cout << "Using hardcoded field..." << endl;
    vector<double> field = {-0.000311776, -0.000320269, -0.000321936, -0.000315862, -0.000301276, -0.000277617, -0.000244615, -0.000202348, -0.000151293, -9.23613e-05, -2.69033e-05, 4.33064e-05, 0.000116116, 0.000189071, 0.000259515, 0.000324709, 0.000381975, 0.000428838, 0.000463172, 0.000483332, 0.000488261, 0.000477569, 0.000451574, 0.000411291, 0.000358382, 0.000295064, 0.000223972, 0.000147999, 7.01187e-05, -6.80498e-06, -8.01877e-05, -0.000147879, -0.000208282, -0.000260426, -0.000303987, -0.000339263, -0.000367085, -0.000388686, -0.000405531, -0.000419121, -0.000430782, -0.000441468, -0.000451578, -0.000460831, -0.000468181, -0.000471818, -0.000469233, -0.000457377, -0.000432878, -0.000392333, -0.00033264, -0.000251345, -0.000146996, -1.94452e-05, 0.000129906, 0.000297955, 0.000479858, 0.000669122, 0.000857826, 0.00103697, 0.00119694, 0.00132804, 0.00142112, 0.00146819, 0.00146299, 0.00140151, 0.00128239, 0.00110717, 0.000880379, 0.000609401, 0.000304183, -2.32333e-05, -0.000359332, -0.000689855, -0.00100058, -0.00127815, -0.0015108, -0.00168905, -0.00180623, -0.00185886, -0.0018468, -0.0017732, -0.0016443, -0.00146892, -0.00125791, -0.00102339, -0.000777921, -0.000533719, -0.000301814, -9.13622e-05, 9.09272e-05, 0.000241176, 0.000358508, 0.000444999, 0.000505391, 0.000546604, 0.000577067, 0.000605923, 0.000642161, 0.000693737, 0.000766761, 0.000864806, 0.0009884, 0.00113474, 0.0012977, 0.001468, 0.00163382, 0.00178144, 0.00189624, 0.0019637, 0.00197052, 0.00190571, 0.00176154, 0.00153433, 0.00122501, 0.000839381, 0.00038804, -0.000113975, -0.000647871, -0.00119204, -0.0017232, -0.00221768, -0.00265273, -0.00300787, -0.00326601, -0.0034145, -0.00344586, -0.00335824, -0.00315554, -0.00284721, -0.00244769, -0.0019756, -0.00145264, -0.000902351, -0.00034875, 0.000185051, 0.000678143, 0.00111297, 0.00147623, 0.00175953, 0.00195968, 0.0020786, 0.002123, 0.00210359, 0.00203415, 0.00193034, 0.00180841, 0.00168389, 0.00157034, 0.00147829, 0.00141443, 0.00138104, 0.00137583, 0.00139217, 0.00141953, 0.00144444, 0.00145149, 0.0014247, 0.00134883, 0.00121078, 0.0010008, 0.000713537, 0.000348797, -8.80539e-05, -0.000586064, -0.00112905, -0.00169635, -0.00226385, -0.00280533, -0.00329396, -0.00370389, -0.00401179, -0.00419832, -0.00424935, -0.00415697, -0.00392001, -0.00354435, -0.00304266, -0.00243385, -0.0017421, -0.000995573, -0.000224913, 0.000538447, 0.00126399, 0.00192378, 0.00249397, 0.002956, 0.00329755, 0.00351303, 0.00360364, 0.00357702, 0.00344651, 0.00323004, 0.0029487, 0.00262524, 0.00228239, 0.00194128, 0.00162002, 0.00133255, 0.00108777, 0.000889101, 0.000734517, 0.000616889, 0.000524804, 0.000443676, 0.000357111, 0.000248416, 0.000102148, -9.44166e-05, -0.000349993, -0.000668284, -0.00104733, -0.00147924, -0.00195044, -0.0024422, -0.00293167, -0.00339321, -0.00379993, -0.00412536, -0.00434522, -0.00443903, -0.00439155, -0.00419399, -0.00384481, -0.00335006, -0.00272339, -0.00198549, -0.00116315, -0.000287898, 0.000605593, 0.00148141, 0.00230435, 0.00304183, 0.00366568, 0.00415369, 0.00449073, 0.00466951, 0.00469079, 0.00456312, 0.00430213, 0.00392927, 0.00347035, 0.00295377, 0.00240855, 0.00186257, 0.00134074, 0.00086354, 0.000445932, 9.66104e-05, -0.000182212, -0.00039457, -0.000550056, -0.000662738, -0.000749758, -0.000829731, -0.000921034, -0.00104014, -0.00120008, -0.0014092, -0.00167022, -0.00197974, -0.00232823, -0.0027004, -0.00307614, -0.00343167, -0.00374119, -0.00397855, -0.00411913, -0.0041416, -0.00402959, -0.00377302, -0.00336915, -0.00282314, -0.00214813, -0.00136484, -0.000500688, 0.000411546, 0.00133545, 0.00223292, 0.00306632, 0.00380057, 0.00440517, 0.00485595, 0.00513642, 0.00523872, 0.00516389, 0.00492176, 0.00453018, 0.00401379, 0.00340236, 0.0027289, 0.00202755, 0.00133147, 0.000670895, 7.13522e-05, -0.00044765, -0.000873564, -0.00120121, -0.0014327, -0.00157681, -0.00164795, -0.00166479, -0.0016485, -0.00162101, -0.0016031, -0.00161271, -0.00166341, -0.00176326, -0.00191403, -0.00211099, -0.00234305, -0.00259354, -0.00284131, -0.00306224, -0.00323099, -0.0033229, -0.00331589, -0.00319224, -0.00294013, -0.00255485, -0.00203953, -0.00140538, -0.000671447, 0.00013619, 0.000985716, 0.00184127, 0.00266496, 0.00341904, 0.00406819, 0.00458162, 0.00493497, 0.00511181, 0.00510466, 0.00491546, 0.00455542, 0.00404429, 0.00340914, 0.00268263, 0.00190101, 0.00110193, 0.00032216, -0.000404543, -0.0010493, -0.00158967, -0.00201072, -0.00230559, -0.00247552, -0.00252938, -0.00248264, -0.00235602, -0.00217376, -0.00196171, -0.00174538, -0.00154807, -0.00138909, -0.00128247, -0.00123588, -0.00125019, -0.00131941, -0.00143124, -0.00156802, -0.00170811, -0.00182761, -0.00190224, -0.00190928, -0.00182949, -0.00164874, -0.0013594, -0.000961244, -0.000461887, 
        0.000123303, 0.000771858, 0.00145537, 0.00214107, 0.0027938, 0.00337819, 0.00386085, 0.00421262, 0.00441048, 0.00443913, 0.00429211, 0.00397235, 0.00349209, 0.00287225, 0.00214123, 0.00133321, 0.000486167, -0.000360368, -0.00116762, -0.0018998, -0.00252607, -0.00302222, -0.00337191, -0.00356734, -0.00360945, -0.00350754, -0.0032784, -0.00294499, -0.00253471, -0.00207758, -0.00160411, -0.00114335, -0.000720973, -0.0003577, -6.80577e-05, 0.000140406, 0.000267401, 0.000319652, 0.000310162, 0.000257013, 0.000181816, 0.000107887, 5.83123e-05, 5.40142e-05, 0.000111964, 0.000243674, 0.000454075, 0.00074086, 0.00109435, 0.00149792, 0.0019289, 0.00236001, 0.0027611, 0.00310123, 0.00335081, 0.00348374, 0.00347936, 0.00332408, 0.00301258, 0.00254846, 0.00194431, 0.00122123, 0.000407709, -0.000461874, -0.0013492, -0.00221419, -0.00301726, -0.00372159, -0.00429521, -0.00471272, -0.00495669, -0.00501851, -0.00489879, -0.0046071, -0.00416131, -0.00358635, -0.00291264, -0.00217421, -0.00140661, -0.000644858, 7.86343e-05, 0.000735851, 0.00130467, 0.00176994, 0.0021241, 0.00236721, 0.00250657, 0.00255578, 0.00253334, 0.00246105, 0.00236209, 0.00225908, 0.00217221, 0.00211757, 0.00210578, 0.00214112, 0.00222106, 0.00233637, 0.00247178, 0.00260704, 0.00271849, 0.00278089, 0.00276941, 0.00266166, 0.00243968, 0.00209155, 0.0016127, 0.00100675, 0.000285723, -0.00053021, -0.00141361, -0.00233113, -0.00324522, -0.00411625, -0.00490467, -0.00557327, -0.00608928, -0.00642628, -0.0065657, -0.00649788, -0.00622269, -0.0057495, -0.0050967, -0.00429065, -0.00336426, -0.00235511, -0.00130344, -0.000249926, 0.000766459, 0.0017105, 0.00255249, 0.00326964, 0.00384699, 0.00427785, 0.0045637, 0.00471352, 0.00474272, 0.00467158, 0.00452349, 0.00432299, 0.00409379, 0.00385693, 0.00362919, 0.00342188, 0.0032401, 0.00308246, 0.00294139, 0.00280392, 0.00265292, 0.00246871, 0.00223091, 0.00192036, 0.00152108, 0.00102194, 0.000418174, -0.000287717, -0.0010854, -0.00195676, -0.00287654, -0.00381339, -0.00473143, -0.00559219, -0.00635667, -0.00698763, -0.00745169, -0.00772139, -0.00777684, -0.00760709, -0.00721084, -0.00659683, -0.00578345, -0.00479799, -0.00367525, -0.00245584, -0.00118402, 9.44905e-05, 0.00133492, 0.0024957, 0.0035405, 0.00443996, 0.00517297, 0.00572741, 0.00610032, 0.00629754, 0.00633275, 0.00622608, 0.00600234, 0.00568898, 0.00531402, 0.00490398, 0.00448201, 0.0040664, 0.00366943, 0.00329687, 0.00294785, 0.0026154, 0.00228743, 0.00194818, 0.00157989, 0.00116476, 0.000686899, 0.000134161, -0.000500222, -0.00121642, -0.00200725, -0.00285789, -0.00374625, -0.0046438, -0.00551697, -0.00632882, -0.00704121, -0.00761694, -0.00802205, -0.00822796, -0.00821338, -0.0079658, -0.00748261, -0.00677161, -0.00585092, -0.00474837, -0.00350024, -0.00214961, -0.000744145, 0.000666228, 0.00203191, 0.00330615, 0.00444742, 0.00542144, 0.00620281, 0.00677602, 0.0071359, 0.00728737, 0.00724466, 0.00702985, 0.00667107, 0.00620029, 0.00565106, 0.00505616, 0.00444543, 0.00384397, 0.00327068, 0.0027374, 0.00224856, 0.00180147, 0.00138709, 0.000991372, 0.000596931, 0.000185013, -0.000262445, -0.000760771, -0.00132072, -0.00194697, -0.00263703, -0.00338074, -0.00416028, -0.00495081, -0.00572167, -0.00643797, -0.00706269, -0.00755885, -0.00789195, -0.00803225, -0.00795684, -0.00765146, -0.00711175, -0.00634405, -0.00536554, -0.00420371, -0.00289531, -0.00148464, -2.14247e-05, 0.00144169, 0.00285186, 0.00415877, 0.0053173, 0.00628981, 0.00704808, 0.00757457, 0.00786308, 0.00791867, 0.00775691, 0.00740244, 0.00688707, 0.00624745, 0.00552255, 0.00475108, 0.00396908, 0.00320777, 0.00249187, 0.00183841, 0.00125616, 0.000745703, 0.000300063, -9.40797e-05, -0.000454734, -0.00080261, -0.00115894, -0.00154327, -0.00197139, -0.00245357, -0.00299313, -0.0035857, -0.00421886, -0.00487258, -0.00552018, -0.00612978, -0.00666629, -0.00709363, -0.00737719, -0.00748616, -0.00739589, -0.00708978, -0.0065608, -0.00581248, -0.00485924, -0.0037261, -0.00244771, -0.00106681, 0.000367826, 0.00180389, 0.00318823, 0.00446971, 0.00560194, 0.00654573, 0.00727121, 0.00775928, 0.00800245, 0.00800497, 0.00778213, 0.00735902, 0.0067686, 0.0060494, 0.00524292, 0.00439097, 0.00353308, 0.00270424, 0.00193294, 0.00123986, 0.000637003, 0.000127592, -0.000293495, -0.000638686, -0.00092617, -0.00117799, -0.00141785, -0.00166889, -0.00195143, -0.00228105, -0.00266695, -0.00311084, -0.00360638, -0.00413925, -0.00468783, -0.00522449, -0.00571726, -0.00613203, -0.00643485, -0.00659434, -0.00658405, -0.00638449, -0.00598481, -0.00538397, -0.00459128, -0.00362636, -0.00251839, -0.0013048, -2.93978e-05, 0.00125989, 0.00251362, 0.00368349, 0.00472508, 0.00560021, 0.00627913, 0.00674202, 0.00697999, 0.00699533, 0.00680107, 0.00641989, 0.00588243, 0.00522514, 
        0.00448786, 0.00371125, 0.00293426, 0.00219182, 0.00151287, 0.0009189, 0.000422986, 2.94396e-05, -0.000265927, -0.000474974, -0.000615898, -0.000711479, -0.000787041, -0.000868257, -0.000978952, -0.00113907, -0.0013629, -0.00165778, -0.00202328, -0.00245094, -0.00292467, -0.00342169, -0.00391401, -0.00437025, -0.00475788, -0.00504546, -0.00520495, -0.00521378, -0.00505663, -0.00472676, -0.00422676, -0.00356875, -0.00277394, -0.0018716, -0.000897497, 0.000108092, 0.00110268, 0.00204404, 0.00289266, 0.00361404, 0.00418074, 0.0045739, 0.00478424, 0.00481248, 0.00466908, 0.00437333, 0.00395198, 0.00343731, 0.00286498, 0.00227168, 0.00169277, 0.00116006, 0.000699983, 0.000332008, 6.76885e-05, -8.98212e-05, -0.000145659, -0.000112603, -9.92521e-06, 0.000138177, 0.000304467, 0.000460734, 0.000579926, 0.000638178, 0.000616626, 0.000502865, 0.000291947, -1.31391e-05, -0.000401561, -0.000855295, -0.0013503, -0.00185811, -0.00234782, -0.00278815, -0.0031497, -0.00340695, -0.00354013, -0.00353665, -0.00339203, -0.0031104, -0.00270427, -0.00219385, -0.00160581, -0.00097159, -0.000325391, 0.000297992, 0.000865434, 0.0013476, 0.00172084, 0.00196868, 0.00208287, 0.00206382, 0.0019205, 0.00166977, 0.00133515, 0.000945176, 0.000531474, 0.000126596, -0.000238144, -0.000534731, -0.000740274, -0.000838483, -0.000820684, -0.000686351, -0.000443089, -0.000106102, 0.000302834, 0.000756814, 0.00122559, 0.00167753, 0.00208175, 0.0024101, 0.00263904, 0.00275121, 0.00273653, 0.00259294, 0.00232646, 0.00195082, 0.00148649, 0.000959343, 0.000398835, -0.000163926, -0.000698287, -0.0011761, -0.00157356, -0.0018728, -0.00206296, -0.00214089, -0.00211119, -0.00198579, -0.00178301, -0.00152616, -0.00124174, -0.00095758, -0.000700687, -0.000495313, -0.000361135, -0.000311788, -0.000353822, -0.000486155, -0.00070008, -0.000979799, -0.00130347, -0.00164468, -0.00197421, -0.0022621, -0.00247963, -0.00260137, -0.00260695, -0.00248253, -0.00222187, -0.0018269, -0.0013078, -0.000682518, 2.41561e-05, 0.000781995, 0.00155706, 0.00231369, 0.00301657, 0.00363286, 0.00413409, 0.00449787, 0.00470914, 0.00476105, 0.00465519, 0.00440143, 0.00401707, 0.0035257, 0.00295547, 0.0023373, 0.00170275, 0.00108206, 0.000502164, -1.48829e-05, -0.000453224, -0.000803997, -0.00106564, -0.00124366, -0.00134994, -0.0014015, -0.00141902, -0.00142498, -0.00144175, -0.00148973, -0.0015855, -0.00174043, -0.00195953, -0.00224085, -0.00257536, -0.00294737, -0.00333545, -0.00371378, -0.00405383, -0.00432631, -0.00450316, -0.00455955, -0.00447573, -0.00423849, -0.00384236, -0.00329026, -0.00259361, -0.00177203, -0.000852394, 0.000132499, 0.00114565, 0.00214785, 0.00309987, 0.00396461, 0.00470912, 0.00530647, 0.00573715, 0.00599001, 0.00606278, 0.00596188, 0.0057018, 0.00530398, 0.0047952, 0.00420577, 0.00356748, 0.00291156, 0.00226672, 0.00165743, 0.00110255, 0.00061442, 0.000198365, -0.0001472, -0.00043021, -0.000663915, -0.00086547, -0.00105424, -0.00124995, -0.00147082, -0.00173181, -0.00204311, -0.00240899, -0.0028271, -0.00328823, -0.00377664, -0.00427088, -0.004745, -0.00517022, -0.00551677, -0.00575592, -0.00586199, -0.00581416, -0.00559815, -0.00520738, -0.0046438, -0.00391808, -0.00304938, -0.00206451, -0.000996664, 0.000116277, 0.00123376, 0.00231472, 0.00331981, 0.00421354, 0.00496616, 0.00555521, 0.00596657, 0.00619504, 0.00624431, 0.00612636, 0.00586041, 0.00547139, 0.00498812, 0.00444131, 0.0038615, 0.00327714, 0.00271282, 0.00218792, 0.00171557, 0.00130222, 0.000947598, 0.000645183, 0.000383189, 0.000145866, -8.48908e-05, -0.000327764, -0.0006002, -0.00091669, -0.00128727, -0.00171634, -0.00220194, -0.00273544, -0.00330179, -0.00388021, -0.00444538, -0.00496893, -0.00542125, -0.00577343, -0.00599918, -0.0060767, -0.00599026, -0.00573141, -0.00529987, -0.00470375, -0.00395942, -0.00309077, -0.00212802, -0.00110617, -6.30616e-05, 0.000962623, 0.00193352, 0.00281559, 0.00357997, 0.00420444, 0.00467446, 0.00498377, 0.00513434, 0.00513588, 0.00500485, 0.00476304, 0.00443587, 0.00405056, 0.00363415, 0.00321179, 0.0028051, 0.00243094, 0.00210048, 0.00181883, 0.00158502, 0.00139248, 0.00122992, 0.00108254, 0.000933479, 0.000765417, 0.0005622, 0.000310365, 4.63257e-07, -0.000371913, -0.000805487, -0.00129307, -0.00182179, -0.00237373, -0.00292694, -0.00345678, -0.00393748, -0.0043438, -0.00465266, -0.00484473, -0.00490571, -0.00482736, -0.00460815, -0.00425345, -0.00377532, -0.00319189, -0.00252632, -0.0018055, -0.00105848, -0.000314845, 0.000396988, 0.00105133, 0.00162661, 0.00210642, 0.00248029, 0.00274397, 0.00289938, 0.00295405, 0.00292032, 0.00281418, 0.00265397, 0.00245892, 0.00224778, 0.00203753, 0.00184228, 0.00167247, 0.00153433, 0.00142974, 0.0013564, 0.00130828, 0.00127633, 0.00124949, 0.00121568, 0.00116298, 0.00108066, 0.000960166, 0.000795882, 0.000585698, 0.000331252, 3.79156e-05, -0.000285517, -0.000627387, -0.000973941, -0.00131025, -0.00162118, -0.00189243, -0.00211144, -0.00226818, -0.0023558, -0.00237101, -0.00231425, -0.00218961, -0.00200446, -0.00176899, -0.00149548, -0.0011975, -0.000889091, -0.000583882, -0.000294349, -3.11357e-05, 0.000197456, 0.000385791, 0.000531031, 0.000633076, 0.000694319, 0.000719252, 0.000713963, 0.000685549, 0.000641525, 0.000589236, 0.00053535, 0.000485437, 0.00044367, 0.000412666, 0.000393451, 0.000385554, 0.000387214, 0.000395664, 0.000407472, 0.000418908, 0.000426311, 0.000426412, 0.000416612, 0.000395185, 0.000361393, 0.00031552, 0.000258819, 0.000193392, 0.000122, 4.78486e-05, -2.56627e-05, -9.51793e-05, -0.00015762, -0.000210368, -0.000251426, -0.000279509, -0.000294088, -0.00029539, -0.000284332, -0.000262432, -0.000231679, -0.00019439, -0.000153053, -0.000110177, -6.81516e-05, -2.91231e-05, 5.10032e-06, 3.31185e-05, 5.39786e-05, 6.71857e-05};
    string ffn = "HARDCODED";
#endif

    cout << "Successfully read fields." << endl;

    EMatrix upper_triangle_ones = EMatrix::Zero();
    for (int i = 1; i < DIM; ++i)
        for (int j = 0; j < i; ++j)
            for (EMatrix& upper : dipoles_upper)
                if (upper(i, j) != upper(j, i))
                    upper_triangle_ones(i, j) = 1.;

 EMatrix encoding_integers;
    enum enc_scheme { other, order, partial, full };
    enum enc_type { hermitian, antihermitian, nonhermitian };
    const enc_scheme cur_scheme = partial;
    const enc_type cur_type = nonhermitian;

    if (cur_scheme == other) {

    } else if (cur_scheme == order) {
        encoding_integers = upper_triangle_ones + (cur_type == nonhermitian) * upper_triangle_ones.transpose();
    } else if (cur_scheme == partial) {
        encoding_integers << 
                0, 2, 3,
                0, 0, 4,
                1, 0, 0;


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
        long double max_ord = 0;
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

    cout << "Running analysis..." << endl;
    vector<CArray> anal_res = run_order_analysis(fields, psi_i, cur_type == hermitian, encoding_integers);

    DipoleSet dipoles = dipoles_upper;
    for (EMatrix &dipole : dipoles) dipole = (dipole + dipole.adjoint()).eval(); 
    auto PGR = gen_pop_graphs(fields, dipoles, psi_i);


    ofstream outfile(string(path) + "HMI_" + to_string(main_start_time) + "_" + ffn.substr(ffn.find_last_of("/\\") + 1) 
                     + (message == "#" || message.empty() ? "" : "_" + message) + ".txt");

    int out_ints[] = {DIM, N_T, main_start_time, L, N_H, N_TO, TIME_POINT_SCALE, N_FIELDS, ORD, BASE, 10 * cur_scheme + cur_type};
    double out_doubles[] = {T, HBAR, field_scale_factor};

    if (message.empty())
        message = "#";
    outfile << "HMR1 " << 1 << ' ' << message << ' ' << time(nullptr) - main_start_time << endl;
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

#ifdef USE_GOTO
    message = message_backup;
    goto start_label; // Yes, this is horrible. Yes, I know I should use a loop, or run the program multiple times, or *anything* else. 
    // But right now, I just want to have main_start_time stay the same without all the work involved, so too bad.
#endif
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

// OArr evolve_initial_hermitian(const FieldSet& fields, const DipoleSet& dipoles, const EVector& psi_i) {
//     throw runtime_error("Not implemented right now");
//     // auto diag_ret = diag_vec(mu, C);
//     // EVector lambda = diag_ret.second;
//     // EMatrix CP = diag_ret.first.first;
//     // EMatrix PdC = diag_ret.first.second;
//     // EMatrix PdCCP = PdC * CP;

//     // vector<EDMatrix, aligned_allocator<EDMatrix>> E(N_T);
//     // for (int i = 0; i < N_T; ++i)
//     //     E[i] = exp(1i * DELTA_T / HBAR * lambda.array() * epsilon[i]).matrix().asDiagonal();

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
        it[i] = (-1.i * Hs[i-1] * DELTA_T).exp() * it[i - 1];

    OArr samples{};
    for (int i = 0; i < N_TO; ++i)
        for (int j = 0; j < DIM; ++j)
            samples[i * DIM + j] = it[(i + 1) * N_T / N_TO][j];
    return samples;
}

vector<CArray> run_order_analysis(const FieldSet& fields, const EVector& psi_i, bool hermitian, const EMatrix& encoding_integers) {
    cout << time(nullptr) << endl;
    int ord = ORD;
    vector<OArr> order_results(ord);
    cout << "Running analysis with ord=" << ord << endl;
#pragma omp parallel for default(shared)
    for (int s = 0; s < ord; ++s) {
        if (s % 1000 == 0)
            cout << "Doing s=" << s << endl;
        double g = 2 * MY_PI * s / ord;
        EMatrix encoding = encoding_integers;
        for (Complex& d : encoding.reshaped())
            d = polar(1.l, d.real() * g);
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
        cout << "\nFFT for 1 to " << i + 1 << ":\n";
        // for (auto& d : tfft)
        //     cout << abs(d) << ' ';
        for (int k = 0; k < ord; ++k) {
            if (abs(tfft[k]) > 0.01)
                cout << k << ": " << abs(tfft[k]) << '\n';
        }
        cout << endl << "Sum of values: " << tfft.sum() << "; magnitude " << abs(tfft.sum()) << "; prob " << norm(tfft.sum()) << '\n';
        ffts.push_back(tfft);
    }
    double sum_of_probs = 0;
    for (CArray& fft : ffts)
        sum_of_probs += norm(fft.sum());
    cout << "\nFinal state when unmodulated: ";
    for (const Complex& c : order_results[0])
        cout << c << ' ';
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
        it[i] = (-1.i * Hs[i-1] * DELTA_T).exp() * it[i - 1];

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
