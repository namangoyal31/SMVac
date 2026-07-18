#pragma once
#include <SMVacuumDecay/RGE.hpp>
#include <tuple>

namespace SMVacuumDecay {

double get_pure_sm_lambda(RGEHelper& rge, double phi_dimless, double mu_inst);

double integrand_u(RGEHelper& rge, double u, double R, double mu_inst, double prefactor);

double evaluate_action_at_R(RGEHelper& rge, double mu_inst, double R);

double find_minimum_action(RGEHelper& rge, double mu_inst, double t_min_lambda);

std::tuple<int, double, double, double> classify_stability(double Mh, double Mt);

}
