// SMVac — vacuum decay probability.
#include "smvac/decay.hpp"
#include "smvac/constants.hpp"

#include <cmath>

namespace smvac {

double hubble_rate_GeV(double H0_km_s_Mpc) { return H0_km_s_Mpc / kMpcKm * kHbarGeVs; }

double action_threshold(double lnLambdaB, const DecayParameters& d) {
    return 4.0 * (lnLambdaB - std::log(hubble_rate_GeV(d.H0_km_s_Mpc))) + std::log(d.light_cone);
}

double log10_decay_probability(double S, double lnLambdaB, const DecayParameters& d) {
    return (action_threshold(lnLambdaB, d) - S) / std::log(10.0);
}

}  // namespace smvac
