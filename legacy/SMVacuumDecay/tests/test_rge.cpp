#include <SMVacuumDecay/RGE.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace SMVacuumDecay;

struct RefPoint {
    double mu, g1, g2, g3, yt, lambda;
};

int main() {
    std::ifstream in("SMVacuumDecay/reference/v1.0/rge_reference.csv");
    if (!in.is_open()) {
        std::cerr << "Could not open rge_reference.csv" << std::endl;
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
        std::getline(ss, token, ','); pt.mu = std::stod(token);
        std::getline(ss, token, ','); pt.g1 = std::stod(token);
        std::getline(ss, token, ','); pt.g2 = std::stod(token);
        std::getline(ss, token, ','); pt.g3 = std::stod(token);
        std::getline(ss, token, ','); pt.yt = std::stod(token);
        std::getline(ss, token, ','); pt.lambda = std::stod(token);
        ref_data.push_back(pt);
    }
    
    double Mh = 125.1;
    double Mt = 173.1;
    double MPlanck = 1.22e19;
    
    StandardModelParameters y_match = get_nnlo_matching(Mh, Mt);
    double t_match = 2 * std::log(Mt);
    
    double t_start = 2 * std::log(1.0);
    double tPlanck = 2 * std::log(MPlanck);
    double dt = 0.1;
    
    RGEHelper rge;
    
    // Backward
    StandardModelParameters y_rev = y_match;
    double t_rev = t_match;
    rge.add_point(t_rev, y_rev);
    while (t_rev > t_start) {
        if (t_rev - dt < t_start) dt = t_rev - t_start;
        y_rev = rk4_single_step(y_rev, t_rev, -dt);
        t_rev -= dt;
        rge.add_point(t_rev, y_rev);
    }
    
    // Forward
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
        double t = 2.0 * std::log(ref.mu);
        StandardModelParameters p = rge.get_params(t);
        
        double err_g1 = std::abs(p.g1 - ref.g1);
        double err_g2 = std::abs(p.g2 - ref.g2);
        double err_g3 = std::abs(p.g3 - ref.g3);
        double err_yt = std::abs(p.yt - ref.yt);
        double err_lam = std::abs(p.lambda - ref.lambda);
        
        double cur_max = std::max({err_g1, err_g2, err_g3, err_yt, err_lam});
        if (cur_max > 1e-10 && max_err < 1e-10) {
            std::cout << "First divergence at mu = " << ref.mu << " (t=" << t << ")\n";
            std::cout << "Ref: g1=" << ref.g1 << " g2=" << ref.g2 << " g3=" << ref.g3 << " yt=" << ref.yt << " lam=" << ref.lambda << "\n";
            std::cout << "Cur: g1=" << p.g1 << " g2=" << p.g2 << " g3=" << p.g3 << " yt=" << p.yt << " lam=" << p.lambda << "\n";
        }
        if (cur_max > max_err) {
            max_err = cur_max;
        }
    }
    
    std::cout << "RGE Validation Max Error: " << std::scientific << max_err << std::endl;
    
    if (max_err < 1e-12) {
        std::cout << "PASS" << std::endl;
        return 0;
    } else {
        std::cout << "FAIL" << std::endl;
        return 1;
    }
}
