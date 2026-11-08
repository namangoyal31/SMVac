// SMVac — vacuum decay probability in our past light cone (Buttazzo et al. 2013, eqs. 80-82):
//   d p = dt dV Lambda_B^4 exp(-S),   p = 0.15 Lambda_B^4 exp(-S) / H0^4,
// with Lambda_B = 1/R the bounce scale. The vacuum is metastable if p < 1.
#pragma once

namespace smvac {

struct DecayParameters {
    double H0_km_s_Mpc = 67.4;   // Hubble rate used by Buttazzo et al. (2013)
    double light_cone = 0.15;    // past-light-cone 4-volume in units of H0^-4 (their eq. 82)
};

double hubble_rate_GeV(double H0_km_s_Mpc);

// log10 p for action S and ln(Lambda_B / GeV).
double log10_decay_probability(double S, double lnLambdaB, const DecayParameters& d = DecayParameters());

// Action at which p = 1 for the given bounce scale: S_thr = 4 ln(Lambda_B / H0) + ln(0.15).
double action_threshold(double lnLambdaB, const DecayParameters& d = DecayParameters());

}  // namespace smvac
