#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/FubiniLipatov.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
using namespace SMVacuumDecay;

int main(int argc, char* argv[]) {
    bool use_analytical = false;
    size_t start_idx = 0;
    size_t end_idx = 0;
    
    bool use_dense = false;
    
    int arg_offset = 1;
    if (argc > 1) {
        std::string arg1(argv[1]);
        if (arg1 == "--analytical") {
            use_analytical = true;
            arg_offset = 2;
        } else if (arg1 == "--numerical") {
            use_analytical = false;
            arg_offset = 2;
        }
    }
    
    if (argc > arg_offset) {
        std::string next_arg(argv[arg_offset]);
        if (next_arg == "--dense") {
            use_dense = true;
            arg_offset++;
        }
    }
    
    if (argc >= arg_offset + 2) {
        start_idx = std::stoull(argv[arg_offset]);
        end_idx = std::stoull(argv[arg_offset + 1]);
    }
    
    vector<pair<double,double>> points;
    double step = use_dense ? 0.03125 : 0.25;
    double Mt_min = use_dense ? 160.0 : 0.0;
    double Mt_max = use_dense ? 185.0 : 250.0;
    double Mh_min = use_dense ? 110.0 : 0.0;
    double Mh_max = use_dense ? 140.0 : 250.0;
    
    for (double Mt = Mt_min; Mt <= Mt_max; Mt += step) {
        for (double Mh_input = Mh_min; Mh_input <= Mh_max; Mh_input += step) {
            points.emplace_back(Mt, Mh_input);
        }
    }
    
    if (end_idx == 0 || end_idx > points.size()) {
        end_idx = points.size();
    }
    
    std::string prefix = use_analytical ? "analytical" : "numerical";
    std::string filename;
    
    if (use_dense) {
        filename = "results/analytical_data_dense_chunk_" + std::to_string(start_idx) + ".csv";
        if (!use_analytical) filename = "results/numerical_data_dense_chunk_" + std::to_string(start_idx) + ".csv";
    } else {
        if (use_analytical) {
            filename = "results/analytical_data_chunk_" + std::to_string(start_idx) + ".csv";
        } else {
            filename = "results/data_numerical_chunk_" + std::to_string(start_idx) + ".csv";
        }
    }
    
    ofstream file(filename);
    if (use_analytical) {
        file << "Mt,Mh_calc,Stability,S_approx\n";
    } else {
        file << "Mt,Mh_calc,Stability,S_exact,S_approx\n";
    }
    
    cout << "Analyzing " << prefix << " points " << start_idx << " to " << end_idx << "..." << endl;
    
    for (size_t i = start_idx; i < end_idx; ++i) {
        double Mt = points[i].first;
        double Mh_input = points[i].second;
        
        if (use_analytical) {
            auto res = classify_buttazzo(Mh_input, Mt);
            int stability = std::get<0>(res);
            double S_approx = std::get<1>(res);
            double Mh_calc = std::get<2>(res);
            file << Mt << "," << Mh_calc << "," << stability << "," << S_approx << "\n";
        } else {
            auto res = classify_stability(Mh_input, Mt);
            int stability = std::get<0>(res);
            double S_exact = std::get<1>(res);
            double S_approx = std::get<2>(res);
            double Mh_calc = std::get<3>(res);
            file << Mt << "," << Mh_calc << "," << stability << "," << S_exact << "," << S_approx << "\n";
        }
        
        if (use_analytical && i % 5000 == 0) cout << "Processed " << i << std::endl;
        else if (!use_analytical && i % 1000 == 0) cout << "Processed " << i << std::endl;
    }
    
    file.close();
    cout << "Done. Saved to " << filename << endl;
    return 0;
}
