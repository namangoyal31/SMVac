#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/FubiniLipatov.hpp>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;
using namespace SMVacuumDecay;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: benchmark_scan <Mh> <Mt>" << endl;
        return 1;
    }
    
    double Mh = std::stod(argv[1]);
    double Mt = std::stod(argv[2]);
    
    cout << std::scientific << std::setprecision(12);
    cout << "--- Benchmark Scan ---" << endl;
    cout << "Mh: " << Mh << " GeV" << endl;
    cout << "Mt: " << Mt << " GeV" << endl;
    
    auto fl_res = classify_buttazzo(Mh, Mt);
    cout << "\n--- Analytical (Fubini-Lipatov) ---" << endl;
    cout << "Status: " << std::get<0>(fl_res) << endl;
    cout << "Action: " << std::get<1>(fl_res) << endl;
    
    auto num_res = classify_stability(Mh, Mt);
    cout << "\n--- Numerical (Canonical) ---" << endl;
    cout << "Status: " << std::get<0>(num_res) << endl;
    cout << "Action: " << std::get<1>(num_res) << endl;
    cout << "S_approx: " << std::get<2>(num_res) << endl;
    
    return 0;
}
