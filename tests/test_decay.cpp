// Tests of the decay-probability criterion (Buttazzo et al. 2013, eqs. 80-82).
#include "smvac/decay.hpp"
#include "test_util.hpp"

#include <cmath>

using namespace smvac;
using smvac_test::check;
using smvac_test::close;

int main() {
    // 1/H0 for H0 = 67.4 km/s/Mpc is 14.5 Gyr (1 Gyr = 3.15576e16 s); hbar = 6.582119569e-25 GeV s.
    const double invH0_s = 6.582119569e-25 / hubble_rate_GeV(67.4);
    close(invH0_s / 3.15576e16, 14.507, 2e-4, 0, "1/H0 = 14.5 Gyr for H0 = 67.4 km/s/Mpc");
    // Threshold S_thr = 4 ln(Lambda_B/H0) + ln 0.15, checked for Lambda_B = 1e17 GeV against the explicit formula.
    const double lnL = std::log(1e17);
    close(action_threshold(lnL), 4.0 * (lnL - std::log(hubble_rate_GeV(67.4))) + std::log(0.15), 1e-15, 0,
          "threshold formula");
    close(action_threshold(lnL), 540.08, 0, 0.05, "threshold at Lambda_B = 1e17 GeV is ~540");
    close(log10_decay_probability(action_threshold(lnL), lnL), 0.0, 0, 1e-12, "p = 1 exactly at the threshold");
    check(log10_decay_probability(2000.0, lnL) < -600.0, "S = 2000 is deeply metastable");
    // Each factor e in Lambda_B shifts the threshold by 4.
    close(action_threshold(lnL + 1.0) - action_threshold(lnL), 4.0, 1e-13, 0, "dS_thr/dln Lambda_B = 4");
    return smvac_test::finish();
}
