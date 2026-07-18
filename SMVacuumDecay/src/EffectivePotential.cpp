#include <SMVacuumDecay/EffectivePotential.hpp>
#include <cmath>

namespace SMVacuumDecay {

double V1_loop(double phi, double yt) {
    double g1init = 0.46;
    double g2init = 0.65;
    
    double Mt2 = 0.5 * yt * yt * phi * phi;
    double Mw2 = 0.25 * g2init * g2init * phi * phi;
    double Mz2 = 0.25 * (g1init * g1init + g2init * g2init) * phi * phi;
    
    double mu = 173.34;
    double mu2 = mu * mu;
    
    double term_t = -12.0 * Mt2 * Mt2 * (std::log(Mt2 / mu2) - 1.5);
    double term_W = 6.0 * Mw2 * Mw2 * (std::log(Mw2 / mu2) - 5.0/6.0);
    double term_Z = 3.0 * Mz2 * Mz2 * (std::log(Mz2 / mu2) - 5.0/6.0);
    
    return (term_t + term_W + term_Z) / (64.0 * PI2);
}

double get_lambda_eff(const StandardModelParameters& p) {
    double yt2 = p.yt * p.yt;
    double yt4 = yt2 * yt2;
    double term_t = -3.0 * yt4 * (std::log(0.5 * yt2) - 1.5);
    
    double g2_2 = p.g2 * p.g2;
    double g2_4 = g2_2 * g2_2;
    double term_W = 0.375 * g2_4 * (std::log(0.25 * g2_2) - 5.0/6.0);
    
    double g12 = p.g1 * p.g1 + g2_2;
    double g12_2 = g12 * g12;
    double term_Z = 0.1875 * g12_2 * (std::log(0.25 * g12) - 5.0/6.0);
    
    double delta_lambda = (term_t + term_W + term_Z) / (16.0 * PI2);
    return p.lambda + delta_lambda;
}

double V_eff(double phi, RGEHelper& rge) {
    if (phi <= 0) return 0.0;
    double t = 2.0 * std::log(phi);
    StandardModelParameters p = rge.get_params(t);
    double lam = get_lambda_eff(p);
    return 0.25 * lam * phi * phi * phi * phi;
}

}
