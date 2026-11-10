#include <SMVacuumDecay/EffectivePotential.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace SMVacuumDecay;

struct RefPoint {
    double phi, lambda_eff, V;
};

int main() {
    std::ifstream in("SMVacuumDecay/reference/v1.0/potential_reference.csv");
    if (!in.is_open()) {
        std::cerr << "Could not open potential_reference.csv" << std::endl;
        return 1;
    }
    
    std::string line;
    std::getline(in, line); // header
    
    std::vector<RefPoint> ref_data;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string token;
        RefPoint pt;
        std::getline(ss, token, ','); pt.phi = std::stod(token);
        std::getline(ss, token, ','); pt.lambda_eff = std::stod(token);
        std::getline(ss, token, ','); pt.V = std::stod(token);
        ref_data.push_back(pt);
    }
    
    double Mh = 5.0;
    double Mt = 105.0;
    double MPlanck = 1.22e19;
    
    StandardModelParameters y_match = get_nnlo_matching(Mh, Mt);
    double t_match = 2 * std::log(Mt);
    
    double t_start = 2 * std::log(1.0);
    double tPlanck = 2 * std::log(MPlanck);
    double dt = 0.1;
    
    RGEHelper rge;
    
    StandardModelParameters y_rev = y_match;
    double t_rev = t_match;
    rge.add_point(t_rev, y_rev);
    while (t_rev > t_start) {
        if (t_rev - dt < t_start) dt = t_rev - t_start;
        y_rev = rk4_single_step(y_rev, t_rev, -dt);
        t_rev -= dt;
        rge.add_point(t_rev, y_rev);
    }
    
    StandardModelParameters y_fwd = y_match;
    double t_fwd = t_match;
    dt = 0.1;
    while (t_fwd < tPlanck) {
        if (t_fwd + dt > tPlanck) dt = tPlanck - t_fwd;
        y_fwd = rk4_single_step(y_fwd, t_fwd, dt);
        t_fwd += dt;
        rge.add_point(t_fwd, y_fwd);
    }
    
    double max_err = 0;
    
    for (const auto& ref : ref_data) {
        double t = 2.0 * std::log(ref.phi);
        StandardModelParameters p = rge.get_params(t);
        double lam_eff = get_lambda_eff(p);
        double V = V_eff(ref.phi, rge);
        
        double err_lam = std::abs(lam_eff - ref.lambda_eff);
        double err_V;
        if (std::abs(ref.V) > 1e-30) {
             err_V = std::abs((V - ref.V) / ref.V);
        } else {
             err_V = std::abs(V - ref.V);
        }
        
        double cur_max = std::max({err_lam, err_V});
        if (cur_max > max_err) {
            max_err = cur_max;
        }
    }
    
    std::cout << "Potential Validation Max Error: " << std::scientific << max_err << std::endl;
    
    if (max_err < 1e-12) {
        std::cout << "PASS" << std::endl;
        return 0;
    } else {
        std::cout << "FAIL" << std::endl;
        return 1;
    }
}
