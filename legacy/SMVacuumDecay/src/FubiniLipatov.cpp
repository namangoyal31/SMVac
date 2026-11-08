#include <SMVacuumDecay/FubiniLipatov.hpp>
#include <SMVacuumDecay/RGE.hpp>
#include <SMVacuumDecay/EffectivePotential.hpp>
#include <cmath>

namespace SMVacuumDecay {

std::tuple<int, double, double> classify_buttazzo(double Mh_input, double Mt) {
    StandardModelParameters y = get_nnlo_matching(Mh_input, Mt);
    double Mh_calc = Mh_input;
    double t0 = 2*std::log(Mt);
    // Standard analytic plots only evaluate up to Planck scale
    double tMax = 2*std::log(1.22e19); 
    
    double min_lambda_eff = 1.0;
    double mu1 = -1.0;
    bool went_negative = false;
    double t = t0;
    double dt = 0.1; // Increased dt for 10x faster integration
    
    double t_v = 2 * std::log(v);
    
    if (t >= t_v) {
        min_lambda_eff = get_lambda_eff(y);
        if (min_lambda_eff <= 0.0) {
            went_negative = true;
            mu1 = std::exp(t / 2.0);
        }
    }
    
    while (t < tMax) {
        if (t + dt > tMax) dt = tMax - t;
        y = rk4_single_step(y, t, dt);
        t += dt;
        
        if (t >= t_v) {
            double lam_eff = get_lambda_eff(y);
            if (!went_negative && lam_eff <= 0.0) {
                mu1 = std::exp(t / 2.0);
                went_negative = true;
            }
            if (lam_eff < min_lambda_eff) min_lambda_eff = lam_eff;
        }
        if (std::abs(y.lambda) > 4*pi || y.yt > 4*pi || y.g1 > 2.0) break;
    }
    
    if (t < 2*std::log(1.22e19)) return std::make_tuple(4, -1.0, Mh_calc);
    if (min_lambda_eff >= 0.0) return std::make_tuple(1, -1.0, Mh_calc);
    
    double S_approx = 8.0 * pi * pi / (3.0 * std::abs(min_lambda_eff));
    double Tv = 1.179e44 / v;
    double S_threshold = 4.0 * std::log(Tv * mu1);
    
    bool is_metastable_S = (S_approx > S_threshold);
    if (is_metastable_S) return std::make_tuple(2, S_approx, Mh_calc);
    else return std::make_tuple(3, S_approx, Mh_calc);
}

}
