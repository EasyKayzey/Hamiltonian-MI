#include "main.h"
#include "omp.h"

double T = 4000, DELTA_T, N_T_double = 1200;
int N_T;
int ORD_R = 512;
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

#define USE_FIELD_FILE true
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
    vector<double> field = {-0.000182329, -0.000146157, -0.000106726, -6.36898e-05, -1.65071e-05, 3.54378e-05, 9.26904e-05, 0.000155568, 0.000224013, 0.000297471, 0.000374819, 0.000454342, 0.000533782, 0.000610428, 0.000681264, 0.000743136, 0.000792926, 0.000827717, 0.000844919, 0.000842371, 0.000818397, 0.00077184, 0.000702086, 0.000609081, 0.000493372, 0.000356165, 0.000199394, 2.5804e-05, -0.000160984, -0.000356447, -0.000555198, -0.000751116, -0.000937574, -0.00110778, -0.0012552, -0.00137406, -0.00145981, -0.00150956, -0.00152246, -0.00149975, -0.00144473, -0.00136245, -0.00125911, -0.00114138, -0.00101554, -0.000886697, -0.000758048, -0.000630421, -0.000502113, -0.000369076, -0.000225478, -6.457e-05, 0.000120225, 0.000334119, 0.000579827, 0.000856635, 0.0011598, 0.00148037, 0.00180557, 0.00211951, 0.00240452, 0.00264259, 0.00281705, 0.00291415, 0.00292432, 0.00284312, 0.00267155, 0.00241579, 0.00208648, 0.00169749, 0.00126451, 0.000803524, 0.00032949, -0.000144662, -0.000608614, -0.00105474, -0.00147771, -0.00187378, -0.00223983, -0.0025724, -0.00286689, -0.00311711, -0.0033152, -0.00345216, -0.00351871, -0.00350652, -0.00340965, -0.00322584, -0.00295757, -0.00261266, -0.00220417, -0.00174964, -0.00126967, -0.000786005, -0.000319259, 0.00011332, 0.000499783, 0.000834818, 0.00112025, 0.00136458, 0.00158162, 0.00178826, 0.00200174, 0.00223669, 0.00250241, 0.00280063, 0.0031242, 0.00345695, 0.0037746, 0.00404703, 0.00424128, 0.00432526, 0.00427162, 0.00406119, 0.00368574, 0.00314951, 0.0024694, 0.00167374, 0.000799762, -0.000109895, -0.00101119, -0.00186275, -0.00262944, -0.00328516, -0.00381432, -0.00421217, -0.0044837, -0.00464159, -0.00470326, -0.00468768, -0.00461221, -0.00448999, -0.00432826, -0.00412776, -0.00388332, -0.00358556, -0.0032233, -0.00278651, -0.00226915, -0.00167158, -0.00100208, -0.000277313, 0.00047854, 0.00123575, 0.00196195, 0.00262569, 0.00319991, 0.00366509, 0.00401146, 0.00423989, 0.00436138, 0.00439512, 0.0043654, 0.00429774, 0.00421483, 0.00413276, 0.00405826, 0.00398713, 0.0039043, 0.00378562, 0.00360097, 0.00331863, 0.00291012, 0.00235504, 0.00164505, 0.000786676, -0.000197777, -0.00127112, -0.00238408, -0.00347957, -0.00449801, -0.00538301, -0.00608676, -0.0065745, -0.00682754, -0.00684461, -0.00664121, -0.0062472, -0.00570297, -0.00505451, -0.00434812, -0.00362541, -0.00291921, -0.00225079, -0.00162892, -0.00105059, -0.000503474, 3.04353e-05, 0.000570516, 0.00113315, 0.00172798, 0.00235526, 0.00300476, 0.00365628, 0.00428168, 0.00484826, 0.00532281, 0.00567595, 0.005886, 0.00594179, 0.00584412, 0.00560547, 0.00524813, 0.0048009, 0.00429491, 0.00375918, 0.00321653, 0.00268053, 0.00215384, 0.00162839, 0.0010872, 0.000507874, -0.000132894, -0.000854282, -0.00166657, -0.00256728, -0.00353894, -0.00454867, -0.00554996, -0.00648632, -0.00729648, -0.00792051, -0.00830609, -0.00841439, -0.00822454, -0.00773643, -0.00697143, -0.00597078, -0.00479196, -0.00350331, -0.00217755, -0.000884776, 0.000313976, 0.00137215, 0.00226117, 0.00297176, 0.00351283, 0.00390826, 0.00419199, 0.00440218, 0.00457513, 0.00474007, 0.00491509, 0.00510508, 0.00530158, 0.00548466, 0.00562644, 0.00569564, 0.0056626, 0.0055039, 0.00520608, 0.00476775, 0.00419997, 0.00352476, 0.00277204, 0.00197545, 0.00116773, 0.000376323, -0.000380078, -0.00109354, -0.00176708, -0.00241303, -0.00304969, -0.00369684, -0.00437061, -0.00507869, -0.00581631, -0.00656399, -0.00728729, -0.00793868, -0.00846161, -0.00879617, -0.00888573, -0.00868382, -0.00816047, -0.0073071, -0.0061395, -0.00469833, -0.00304708, -0.00126746, 0.000547337, 0.00230078, 0.0039018, 0.00527313, 0.0063582, 0.00712581, 0.00757181, 0.00771787, 0.00760727, 0.00729853, 0.00685764, 0.00634991, 0.00583255, 0.00534877, 0.00492404, 0.00456492, 0.00426049, 0.00398584, 0.00370742, 0.00338908, 0.00299823, 0.00251118, 0.00191687, 0.0012187, 0.00043406, -0.000408067, -0.00127119, -0.00211634, -0.00290727, -0.00361509, -0.00422164, -0.00472113, -0.00511974, -0.00543343, -0.00568408, -0.00589471, -0.00608437, -0.00626363, -0.00643129, -0.00657281, -0.00666096, -0.00665852, -0.00652282, -0.00621153, -0.00568903, -0.00493242, -0.00393648, -0.00271676, -0.00131042, 0.000225384, 0.00181768, 0.00338375, 0.00483878, 0.00610388, 0.0071139, 0.00782366, 0.00821213, 0.00828367, 0.00806646, 0.00760817, 0.00696956, 0.00621687, 0.00541418, 0.00461648, 0.00386448, 0.00318165, 0.00257367, 0.00203042, 0.00152982, 0.00104328, 0.000541584, 7.24796e-07, -0.000593424, -0.00124205, -0.00193243, -0.00263916, -0.00332728, -0.00395694, -0.00448877, -0.0048893, -0.00513555, -0.00521824, -0.00514304, -0.00492983, -0.00460994, -0.00422178, -0.0038055, -0.00339747, -0.00302517, -0.00270354, -0.00243307, -0.00220005, -0.00197899, -0.00173689, -0.00143866, -0.00105315, -0.000558712, 5.23386e-05, 0.000771046, 0.00157139, 0.0024122, 0.00324117, 0.00400068, 0.00463444, 0.00509438, 0.00534658, 0.00537562, 0.00518664, 0.00480475, 0.00427199, 0.003642, 0.00297341, 0.00232254, 0.00173663, 0.00124827, 0.000871866, 0.000602607, 0.000417904, 0.000281288, 0.000148145, -2.73964e-05, -0.000285361, -0.000653556, -0.0011428, -0.00174442, -0.0024304, -0.00315608, -0.00386521, -0.0044964, -0.00499048, -0.00529741, -0.00538229, -0.00522942, -0.0048441, -0.00425203, -0.00349622, -0.00263202, -0.00172074, -0.000822799, 8.90973e-06, 0.000734351, 0.00132989, 0.00178951, 0.00212359, 0.00235554, 0.00251669, 0.00264038, 0.00275595, 0.00288363, 0.00303103, 0.00319168, 0.00334601, 0.00346439, 0.00351203, 0.00345489, 0.00326582, 0.00292988, 0.00244819, 0.00183958, 0.00113972, 0.00039783, -0.000328761, -0.000981374, -0.00150734, -0.00186692, -0.00203887, -0.0020239, -0.00184554, -0.00154808, -0.00119186, -0.000846232, -0.000581076, -0.000457833, -0.000521177, -0.000792444, -0.0012657, -0.00190705, -0.00265728, -0.00343763, -0.00415795, -0.00472628, -0.0050588, -0.00508903, -0.00477523, -0.0041054, -0.00309928, -0.00180722, -0.000306025, 0.00130776, 0.00292654, 0.00444138, 0.00575172, 0.00677418, 0.00744924, 0.00774533, 0.00765984, 0.0072171, 0.00646387, 0.00546307, 0.00428661, 0.00300846, 0.00169856, 0.000418287, -0.000782305, -0.00186537, -0.0028046, -0.00358334, -0.00419204, -0.00462571, -0.00488182, -0.00495908, -0.00485739, -0.00457889, -0.00412988, -0.00352325, -0.0027808, -0.0019349, -0.00102887, -0.000115867, 0.000744135, 0.00148885, 0.0020595, 0.0024071, 0.00249858, 0.00232165, 0.00188789, 0.00123329, 0.000416212, -0.00048731, -0.00139032, -0.00220315, -0.00284207, -0.00323759, -0.00334141, -0.00313134, -0.00261367, -0.00182274, -0.00081782, 0.000322604, 0.0015083, 0.00264599, 0.00364823, 0.00444155, 0.00497294, 0.00521387, 0.00516153, 0.00483712, 0.00428162, 0.00354969, 0.00270268, 0.00180162, 0.000901212, 4.54301e-05, -0.000734802, -0.00142181, -0.0020094, -0.0024997, -0.00289907, -0.00321403, -0.00344775, -0.00359779, -0.00365545, -0.00360686, -0.00343561, -0.00312664, -0.00267054, -0.00206781, -0.00133209, -0.000491885, 0.000409704, 0.00131746, 0.00216856, 0.00289865, 0.00344884, 0.00377265, 0.00384188, 0.00365077, 0.00321748, 0.00258291, 0.00180672, 0.000961151, 0.00012333, -0.00063304, -0.00124555, -0.00166962, -0.00188303, -0.0018877, -0.00170879, -0.00139098, -0.000992576, -0.000577929, -0.000209344, 6.06994e-05, 0.000195908, 0.000180332, 1.97593e-05, -0.00025982, -0.000616446, -0.000998512, -0.00135228, -0.00162953, -0.00179432, -0.00182786, -0.00173094, -0.00152352, -0.00124157, -0.000931507, -0.000643136, -0.000421851, -0.000301361, -0.000297888, -0.000406741, -0.000601785, -0.000837967, -0.0010566, -0.00119278, -0.00118394, -0.000978591, -0.000543953, 0.000128292, 0.00101906, 0.00208206, 0.00324733, 0.00442768, 0.00552749, 0.00645271, 0.00712105, 0.00747095, 0.00746814, 0.0071091, 0.00642079, 0.00545685, 0.00429089, 0.00300769, 0.00169362, 0.000427455, -0.000727374, -0.00172819, -0.00255564, -0.00321331, -0.00372438, -0.00412608, -0.00446255, -0.00477703, -0.00510462, -0.00546636, -0.00586553, -0.00628653, -0.00669668, -0.00705029, -0.00729484, -0.00737802, -0.00725499, -0.00689471, -0.00628469, -0.00543344, -0.00437046, -0.00314366, -0.00181464, -0.000452331, 0.000874113, 0.00210238, 0.00318305, 0.00408384, 0.00479162, 0.00531186, 0.00566582, 0.00588578, 0.00600921, 0.00607272, 0.00610653, 0.0061304, 0.00615135, 0.00616351, 0.00614993, 0.00608614, 0.00594469, 0.00570008, 0.00533335, 0.00483554, 0.00420955, 0.0034702, 0.00264242, 0.00175793, 0.000850907, -4.665e-05, -0.00090911, -0.00172036, -0.00247531, -0.0031796, -0.00384734, -0.00449732, -0.00514818, -0.00581314, -0.00649532, -0.00718417, -0.0078538, -0.00846355, -0.00896077, -0.00928579, -0.00937836, -0.0091849, -0.00866575, -0.00780153, -0.00659767, -0.00508665, -0.00332754, -0.00140254, 0.000588999, 0.00253879, 0.00433913, 0.00589274, 0.00712186, 0.00797544, 0.0084335, 0.00850822, 0.00824145, 0.00769911, 0.00696325, 0.00612275, 0.00526393, 0.00446212, 0.00377514, 0.00323923, 0.00286776, 0.00265258, 0.00256765, 0.00257429, 0.00262727, 0.00268088, 0.00269425, 0.0026352, 0.00248237, 0.00222528, 0.00186284, 0.00140047, 0.000846674, 0.000209761, -0.0005047, -0.00129503, -0.00216299, -0.00311177, -0.00414263, -0.00525087, -0.00642177, -0.00762726, -0.00882415, -0.00995424, -0.0109467, -0.0117226, -0.0122015, -0.0123086, -0.0119832, -0.0111864, -0.00990696, -0.00816608, -0.00601852, -0.00355119, -0.000878314, 0.00186634, 0.00453966, 0.00700091, 0.00912383, 0.0108077, 0.0119861, 0.0126319, 0.0127588, 0.0124174, 0.0116884, 0.0106722, 0.00947802, 0.00821188, 0.00696679, 0.00581499, 0.00480335, 0.00395217, 0.0032571, 0.00269382, 0.00222463, 0.00180599, 0.0013959, 0.00096027, 0.000477396, -5.99684e-05, -0.000645026, -0.00125966, -0.00187886, -0.00247613, -0.0030289, -0.00352306, -0.00395596, -0.00433726, -0.00468759, -0.00503507, -0.00541022, -0.00583997, -0.00634172, -0.00691825, -0.00755446, -0.00821633, -0.00885258, -0.00939867, -0.00978307, -0.0099347, -0.00979102, -0.00930554, -0.00845397, -0.00723835, -0.00568843, -0.00386037, -0.00183264, 0.000300383, 
    0.00243636, 0.00447382, 0.00632068, 0.00790173, 0.00916397, 0.0100795, 0.0106455, 0.0108813, 0.0108239, 0.0105211, 0.0100251, 0.00938593, 0.0086469, 0.00784124, 0.00699134, 0.00610984, 0.00520237, 0.00427136, 0.00332005, 0.00235608, 0.00139397, 0.000455887, -0.000429224, -0.00122834, -0.00190835, -0.0024408, -0.00280669, -0.00300054, -0.00303302, -0.0029317, -0.00273951, -0.00251102, -0.00230686, -0.00218681, -0.00220255, -0.00239088, -0.00276837, -0.00332837, -0.00404053, -0.00485333, -0.00569899, -0.00650048, -0.00717948, -0.00766463, -0.00789884, -0.00784508, -0.00748985, -0.00684403, -0.00594116, -0.00483336, -0.00358528, -0.00226717, -0.000947548, 0.000313355, 0.00146872, 0.00248821, 0.00335862, 0.00408224, 0.0046732, 0.0051524, 0.005542, 0.00586017, 0.00611706, 0.00631247, 0.00643568, 0.00646741, 0.00638362, 0.00616055, 0.0057802, 0.00523551, 0.00453413, 0.00370033, 0.00277466, 0.00181106, 0.00087199, 2.20098e-05, -0.000679374, -0.00118465, -0.00146421, -0.00151037, -0.00133908, -0.000988949, -0.000517612, 4.17137e-06, 0.000500088, 0.000896654, 0.00113139, 0.00115978, 0.000960128, 0.000535669, -8.62668e-05, -0.000857944, -0.00171718, -0.00259477, -0.00342234, -0.00413974, -0.00470085, -0.00507749, -0.00526069, -0.00525954, -0.00509789, -0.0048093, -0.00443122, -0.00399912, -0.0035415, -0.00307643, -0.00261015, -0.00213788, -0.00164662, -0.00111964, -0.000541633, 9.59017e-05, 0.000790311, 0.00152496, 0.00226857, 0.00297703, 0.00359741, 0.00407397, 0.00435532, 0.00440169, 0.00419155, 0.0037263, 0.00303245, 0.00216084, 0.00118285, 0.000183898, -0.000744761, -0.00151531, -0.00205243, -0.00230113, -0.0022327, -0.00184781, -0.00117676, -0.000276602, 0.00077443, 0.00188471, 0.00295831, 0.0039045, 0.00464643, 0.00512816, 0.00531912, 0.00521548, 0.00483853, 0.00423023, 0.00344686, 0.00255153, 0.00160673, 0.000667883, -0.000221432, -0.00103203, -0.00174931, -0.00237132, -0.00290516, -0.00336229, -0.00375381, -0.00408623, -0.00435858, -0.00456132, -0.00467717, -0.00468382, -0.00455812, -0.00428104, -0.0038427, -0.00324664, -0.00251255, -0.00167701, -0.000791829, 7.98568e-05, 0.00086954, 0.00151056, 0.00194577, 0.00213468, 0.00205905, 0.00172614, 0.00116911, 0.00044444, -0.000373353, -0.00119911, -0.0019456, -0.00253229, -0.00289332, -0.00298386, -0.00278412, -0.00230074, -0.00156553, -0.000631618, 0.000432376, 0.00154948, 0.0026421, 0.00363931, 0.00448307, 0.00513251, 0.00556599, 0.00578071, 0.00579019, 0.00562002, 0.00530273, 0.00487248, 0.00436039, 0.00379123, 0.00318168, 0.00254049, 0.00187019, 0.00117013, 0.000440042, -0.000316311, -0.00108827, -0.00185663, -0.00259356, -0.00326426, -0.00383019, -0.00425334, -0.0045012, -0.00455152, -0.00439632, -0.00404448, -0.00352251, -0.00287323, -0.00215237, -0.00142351, -0.000751764, -0.000196989, 0.000192636, 0.000385896, 0.000371978, 0.000161608, -0.000214246, -0.00070794, -0.00126052, -0.00180819, -0.00228905, -0.00264935, -0.00284851, -0.00286244, -0.00268481, -0.00232632, -0.00181217, -0.00117813, -0.000465734, 0.000282612, 0.00102812, 0.00173908, 0.00239292, 0.0029767, 0.00348605, 0.00392291, 0.00429242, 0.00459967, 0.00484672, 0.00503053, 0.00514217, 0.00516736, 0.00508836, 0.004887, 0.00454818, 0.00406362, 0.00343484, 0.00267525, 0.00181065, 0.000878172, -7.63754e-05, -0.0010025, -0.00184976, -0.00257286, -0.00313624, -0.00351782, -0.00371133, -0.00372691, -0.00358994, -0.00333813, -0.00301731, -0.00267625, -0.00236128, -0.00211116, -0.00195304, -0.00189985, -0.00194943, -0.00208553, -0.00228045, -0.00249895, -0.00270278, -0.00285542, -0.00292617, -0.00289337, -0.00274617, -0.00248496, -0.00212025, -0.00167048, -0.00115905, -0.000610899, -4.94131e-05, 0.000506213, 0.0010426, 0.00155273, 0.00203529, 0.0024931, 0.00293087, 0.00335254, 0.00375899, 0.00414609, 0.00450385, 0.0048165, 0.00506376, 0.0052231, 0.00527261, 0.00519422, 0.00497671, 0.00461807, 0.00412687, 0.00352234, 0.00283315, 0.00209499, 0.00134726, 0.000629195, -2.39005e-05, -0.000584128, -0.00103356, -0.0013657, -0.00158572, -0.0017094, -0.00176094, -0.00176996, -0.00176794, -0.00178465, -0.00184483, -0.00196562, -0.00215492, -0.00241091, -0.00272271, -0.00307192, -0.00343507, -0.00378625, -0.00409997, -0.00435351, -0.00452886, -0.0046138, -0.00460223, -0.00449374, -0.00429254, -0.00400616, -0.00364387, -0.00321536, -0.00272969, -0.00219472, -0.0016171, -0.0010027, -0.000357384, 0.000311961, 0.000996497, 0.00168469, 0.002362, 0.003011, 0.00361217, 0.00414505, 0.00458986, 0.00492925, 0.00515008, 0.00524484, 0.0052127, 0.00505988, 0.00479931, 0.00444966, 0.00403377, 0.00357671, 0.00310362, 0.00263775, 0.0021986, 0.00180074, 0.001453, 0.00115835, 0.000914413, 0.000714286, 0.000547852, 0.000403184, 0.000268, 0.000130964, -1.72571e-05, -0.000183281, -0.00037066, -0.000579893, -0.000808731, -0.00105271, -0.00130579, -0.00156107, -0.00181134, -0.00204965, -0.00226964, -0.00246571, -0.00263313, -0.00276801, -0.00286719, -0.00292822, -0.00294926, -0.00292914, -0.00286736, -0.00276425, -0.00262105, -0.00244005, -0.00222467, -0.00197942, -0.00170983, -0.00142227, -0.00112368, -0.00082124, -0.000521962, -0.000232362, 4.19017e-05, 0.000296291, 0.000527559};
#endif


    cout << time(nullptr) << endl;
    vector<OArr> order_results(ORD_R);
#pragma omp parallel for default(shared)
    for (int s = 0; s < ORD_R; ++s) {
        EMatrix mu_upper = mu_t_upper;
        double g = -2 * M_PI * s / ORD_R;
        Complex m = polar(1., g);
        mu_upper *= m;
        order_results[s] = evolve_initial_nonhermitian(field, to_full_matrix(mu_upper), psi_i, anal_pop);
    }
    cout << time(nullptr) << endl;

    vector<CArray> ffts;
    for (int i = 0; i < DIM; ++i) {
        int ii = i + DIM;
        CArray tfft(ORD_R);
        for (int j = 0; j < ORD_R; ++j) {
            tfft[j] = order_results[j][ii];
        }
        ifft(tfft);
        cout << "\nFFT for 1 to " << i + 1 << ':' << endl;
        // for (auto& d : tfft)
        //     cout << abs(d) << ' ';
        for (int k = 0; k < ORD_R; ++k) {
            if (abs(tfft[k]) > 0.01)
                cout << k << ": " << abs(tfft[k]) << endl;
        }
        cout << endl << "Sum of values: " << tfft.sum() << "; magnitude " << abs(tfft.sum()) << endl;
        ffts.push_back(tfft);
    }

    ofstream outfile(string(path) + "HMI_" + to_string(main_start_time) + "_" + ffn + (message == "#" ? "" : "_" + message) + ".txt");

    int out_ints[] = {DIM, N_T, main_start_time, L, N_H, N_TO, N_OBS, ORD_R};
    double out_doubles[] = {T, HBAR};

    if (message.empty())
        message = "#";
    outfile << "HMI1 " << 1 << ' ' << message << ' ' << time(nullptr) - main_start_time << endl;
    for (int o : out_ints)
        outfile << ' ' << o;
    outfile << endl;
    for (double o : out_doubles)
        outfile << ' ' << o;
    outfile << endl;

    outfile << "Preliminaries:" << endl;

    outfile << H0D.real().transpose() << endl;

    outfile << mu_t_upper + mu_t_upper.transpose() << endl;

    outfile << psi_i.real().transpose() << endl;

    outfile << "Field:" << endl;

    for (double d : field)
        outfile << d << ' ';
    outfile << endl;

    outfile << "Data:" << endl;

    for (CArray& arr : ffts) {
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

OArr evolve_initial_hermitian(const vector<double>& epsilon, const EMatrix& mu, const EVector& psi_i, const array<ECovector, DIM>& anal_pop) {
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
            samples[i * DIM + j] = it[(i + 1) * N_T / N_TO][j];
    return samples;
}

OArr evolve_initial_nonhermitian(const vector<double>& epsilon, const EMatrix& mu, const EVector& psi_i, const array<ECovector, DIM>& anal_pop) {
    // auto diag_ret = diag_vec(mu, C);
    // EVector lambda = diag_ret.second;       
    // EMatrix CP = diag_ret.first.first;
    // EMatrix PdC = diag_ret.first.second;
    // EMatrix PdCCP = PdC * CP;
    EDMatrix CC = C.toDenseMatrix().diagonal().array().pow(2).matrix().asDiagonal();

    vector<EMatrix, aligned_allocator<EMatrix>> Ex(N_T);
    for (int i = 0; i < N_T; ++i)
        Ex[i] = (1i * DELTA_T / HBAR * mu * epsilon[i]).exp();
    
    vector<EVector, aligned_allocator<EVector>> it(N_T + 1);
    it[0] = C * psi_i;
    for (int i = 1; i < N_T; ++i)
        it[i] = CC * (Ex[i - 1] * it[i - 1]);
    it[N_T] = C * (Ex[N_T - 1] * it[N_T - 1]);

    OArr samples{};
    for (int i = 0; i < N_TO; ++i)
        for (int j = 0; j < DIM; ++j)
            samples[i * DIM + j] = it[(i + 1) * N_T / N_TO][j];
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
    // throw runtime_error("Not implemented");
    array<double, L> omega;

    for (int l = 1, c = 0; l < DIM; ++l)
        for (int m = 0; m < l; ++m)
            omega[c++] = (H0D(l) - H0D(m)).real() / HBAR;
    
    vector<double> eps_inter(N_T);
    for (int i = 0; i < N_T; ++i) {
        eps_inter[i] = 0;
        double t_i = (i + 0.5) * DELTA_T;
        int c = 0;
        for (int l = 1; l < DIM; ++l) {
            for (int m = 0; m < l; ++m) {
                eps_inter[i] += epsilon[c].first * sin(omega[c] * t_i + epsilon[c].second);
                ++c;
            }
        }
        eps_inter[i] *= envelope_funct(t_i);
    }
    return eps_inter;
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
