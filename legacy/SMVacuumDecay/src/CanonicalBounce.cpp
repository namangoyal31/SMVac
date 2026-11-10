#include <SMVacuumDecay/CanonicalBounce.hpp>
#include <SMVacuumDecay/EffectivePotential.hpp>
#include <SMVacuumDecay/numerics/GoldenSectionSearch.hpp>
#include <SMVacuumDecay/numerics/Simpson.hpp>
#include <cmath>
#include <algorithm>

namespace SMVacuumDecay {

static constexpr double CLAUDE_M_PLANCK = 1.22e19;
static constexpr double ALPHA      = 4.0;
static constexpr double SIGMA_PLAN = 8.0;
static constexpr int    N_QUAD     = 2048;

double get_pure_sm_lambda(RGEHelper& rge, double phi_dimless, double mu_inst) {
    double h = phi_dimless * mu_inst;
    double h_rge = std::min(h, CLAUDE_M_PLANCK);
    if (h_rge <= 246.22) return 0.0;

    double t      = 2.0 * std::log(h_rge);
    StandardModelParameters p      = rge.get_params(t);
    return get_lambda_eff(p);
}

double integrand_u(RGEHelper& rge, double u, double R, double mu_inst, double prefactor) {
    if (std::abs(u) >= 1.0) return 0.0;
    double ratio = (1.0 + u) / (1.0 - u);
    double e2x   = std::pow(ratio, ALPHA);
    double e4x   = e2x * e2x;
    double phi_x  = prefactor * 2.0 / (R * (e2x + 1.0));
    
    double lam_phi = get_pure_sm_lambda(rge, phi_x, mu_inst);
    double phi4 = std::pow(phi_x, 4);
    
    double V_SM = 0.25 * lam_phi * phi4;
    
    double lambda_6 = 1.0;
    double V_6 = (lambda_6 / 6.0) * std::pow(phi_x * mu_inst / CLAUDE_M_PLANCK, 2) * phi4;
    
    double V_x = V_SM + V_6;
    
    double jacobian = ALPHA / (1.0 - u * u);
    // Use pi*pi instead of LOOP1 (which is 16*pi^2)
    return 2.0 * pi * pi * std::pow(R, 4) * e4x * V_x * jacobian;
}

double evaluate_action_at_R(RGEHelper& rge, double mu_inst, double R) {
    if (R <= 0.0) return 1e100;
    double mu_R = mu_inst / R;
    if (mu_R <= 0.0) return 1e100;
    double t_R    = 2.0 * std::log(mu_R);
    StandardModelParameters p_R    = rge.get_params(t_R);
    double lambda_R = p_R.lambda;
    if (lambda_R >= 0.0) return 1e100;
    double abs_lam  = std::abs(lambda_R);
    double prefactor = std::sqrt(2.0 / abs_lam);
    // The kinetic integral for the conformal bounce is exactly 16 pi^2 / (3 |\lambda|)
    double kinetic_term = (16.0 * pi * pi) / (3.0 * abs_lam);
    if (kinetic_term > 5e5) return kinetic_term;
    auto f = [&](double u) {
        return integrand_u(rge, u, R, mu_inst, prefactor);
    };
    double potential_integral = numerics::simpson_integrate(-1.0, 1.0, N_QUAD, f);
    
    return kinetic_term + potential_integral;
}

double find_minimum_action(RGEHelper& rge, double mu_inst, double t_min_lambda) {
    double mu_opt = std::exp(t_min_lambda / 2.0);
    double R_opt = mu_inst / mu_opt;
    double logR_opt = std::log(R_opt);
    double a = logR_opt - 8.0;  
    double b = logR_opt + 8.0;  
    double tol = 1e-13;
    
    auto f = [&](double logR) {
        return evaluate_action_at_R(rge, mu_inst, std::exp(logR));
    };
    
    return numerics::golden_section_search(a, b, tol, f);
}

std::tuple<int, double, double, double> classify_stability(double Mh, double Mt) {
    StandardModelParameters y = get_nnlo_matching(Mh, Mt);
    double Mh_calc = Mh; 
    
    double t0 = 2*std::log(Mt);
    double tPlanck = 2*std::log(CLAUDE_M_PLANCK);
    
    RGEHelper rge;
    bool is_unstable = false;
    double mu_inst = CLAUDE_M_PLANCK;
    
    double min_lambda_eff = 1.0;
    double t_min_lambda = t0;
    
    double t = t0;
    double dt = 0.1;
    rge.add_point(t, y);
    
    double t_v = 2 * std::log(246.22);
    
    if (t >= t_v) {
        min_lambda_eff = get_lambda_eff(y);
        if (min_lambda_eff <= 0.0) {
            mu_inst = std::exp(t/2.0);
            is_unstable = true;
        }
    }
    
    while (t < tPlanck) {
        if (t + dt > tPlanck) dt = tPlanck - t;
        
        y = rk4_single_step(y, t, dt);
        t += dt;
        rge.add_point(t, y);
        
        if (std::abs(y.lambda) > 4*pi || y.yt > 4*pi) 
            return std::make_tuple(4, -1.0, -1.0, Mh_calc);
            
        if (!std::isfinite(y.lambda) || !std::isfinite(y.yt))
            return std::make_tuple(4, -1.0, -1.0, Mh_calc);
            
        if (t >= t_v) {
            double lam_eff = get_lambda_eff(y);
            if (lam_eff < min_lambda_eff) {
                min_lambda_eff = lam_eff;
                t_min_lambda = t;
            }
            if (!is_unstable && lam_eff < 0) {
                mu_inst = std::exp(t/2.0);
                is_unstable = true;
            }
        }
    }
    
    if (!is_unstable) return std::make_tuple(1, -1.0, -1.0, Mh_calc); // Stable
    
    double S_approx = 8.0 * pi * pi / (3.0 * std::abs(min_lambda_eff));
    double S_exact = find_minimum_action(rge, mu_inst, t_min_lambda);
    
    // Exact Metastability Threshold using Universe Age
    double Tv = 1.179e44 / 246.22;
    double S_threshold = 4.0 * std::log(Tv * mu_inst);
    
    if (S_exact > S_threshold) return std::make_tuple(2, S_exact, S_approx, Mh_calc); // Metastable
    return std::make_tuple(3, S_exact, S_approx, Mh_calc); // Unstable
}

}
