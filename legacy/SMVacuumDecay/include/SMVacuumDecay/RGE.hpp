#pragma once

#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>

namespace SMVacuumDecay {

// Constants
const double pi = 3.14159265358979323846;
const double PI2 = pi * pi;
const double LOOP1 = 16.0 * PI2;
const double LOOP2 = LOOP1 * LOOP1;
const double LOOP3 = LOOP2 * LOOP1;
const double LOOP4 = LOOP2 * LOOP2;
const double v = 246.22;
const double Mtau = 1.777;
const double Mb = 4.0;
const double alpha3_at_Mz = 0.1184;

struct StandardModelParameters { 
    double g1, g2, g3, yt, yb, ytau, lambda, phi = 0; 
    
    StandardModelParameters operator+(const StandardModelParameters& other) const {
        return {g1+other.g1, g2+other.g2, g3+other.g3, yt+other.yt, yb+other.yb, ytau+other.ytau, lambda+other.lambda, phi+other.phi};
    }
    StandardModelParameters operator*(double scalar) const {
        return {g1*scalar, g2*scalar, g3*scalar, yt*scalar, yb*scalar, ytau*scalar, lambda*scalar, phi*scalar};
    }
};

StandardModelParameters get_nnlo_matching(double Mh, double Mt);

double betaG1sq(const StandardModelParameters& p);
double betaG2sq(const StandardModelParameters& p);
double betaG3sq(const StandardModelParameters& p);
double betaLambda(const StandardModelParameters& p);
double betaYt2(const StandardModelParameters& p);
double betaYb2(const StandardModelParameters& p);
double betaYtau2(const StandardModelParameters& p);

StandardModelParameters rk4_single_step(const StandardModelParameters& y, double t, double dt);
bool rk4_adaptive_step(StandardModelParameters& y, double& t, double& dt);

class RGEHelper {
private:
    std::vector<std::pair<double, StandardModelParameters>> table;
    bool sorted = false;
public:
    void clear();
    void add_point(double t, const StandardModelParameters& p);
    void sort_table();
    StandardModelParameters get_params(double t);
};

} // namespace SMVacuumDecay
