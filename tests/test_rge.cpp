// Tests of the RG evolution.
//  1. Couplings at M_Pl reproduce Buttazzo et al. (2013) eq. (61) to its quoted precision (literature reference).
//  2. Classical RK4 shows fourth-order convergence (error ratio ~16 per step halving).
//  3. Cubic Hermite dense output reproduces an independent finer integration between nodes.
//  4. Non-perturbative inputs are reported as such (no silent failure).
#include "smvac/rge.hpp"
#include "test_util.hpp"

#include <cmath>

using namespace smvac;
using smvac_test::check;
using smvac_test::close;

int main() {
    Inputs in;   // M_h = 125.15, M_t = 173.34, alpha_s = 0.1184, M_W = 80.384 (central values of eq. 61)
    const RGSolution rg(in, 0.025);
    check(rg.perturbative(), "central point runs perturbatively to M_Pl");
    const Couplings c = rg.at(rg.t_end());
    // eq. (61) quotes four decimals; allow one unit in the last place (plus the y_b, y_tau input difference).
    close(c.g1, 0.6154, 0.0, 1e-4, "g1(M_Pl) vs eq. (61)");
    close(c.g2, 0.5055, 0.0, 1e-4, "g2(M_Pl) vs eq. (61)");
    close(c.g3, 0.4873, 0.0, 1e-4, "g3(M_Pl) vs eq. (61)");
    close(c.yt, 0.3825, 0.0, 1e-4, "yt(M_Pl) vs eq. (61)");
    close(c.lambda, -0.0143, 0.0, 1e-4, "lambda(M_Pl) vs eq. (61)");
    // Independent Python reference (DOP853, rtol 1e-13, separate transcription): lambda(M_Pl) = -0.0143688...
    close(c.lambda, -0.014368795660, 0.0, 1e-10, "lambda(M_Pl) vs independent Python reference (DOP853, separate transcription)");
    close(c.yt, 0.3824989140, 0.0, 1e-10, "yt(M_Pl) vs independent Python reference");

    // Fourth-order convergence: errors of dt = 0.2 and 0.1 relative to dt = 0.0125.
    const double ref = RGSolution(in, 0.0125).at(rg.t_end()).lambda;
    const double e1 = std::abs(RGSolution(in, 0.2).at(rg.t_end()).lambda - ref);
    const double e2 = std::abs(RGSolution(in, 0.1).at(rg.t_end()).lambda - ref);
    const double ratio = e1 / e2;
    check(ratio > 12.0 && ratio < 20.0, "RK4 error ratio for step halving is ~16 (got " + std::to_string(ratio) + ")");
    check(e2 < 1e-12, "RK4 with dt = 0.1 is already accurate to < 1e-12 in lambda(M_Pl)");

    // Dense output: midpoints of the dt = 0.1 grid against the dt = 0.025 solution at the same t.
    const RGSolution coarse(in, 0.1);
    double worst = 0.0;
    for (double t = coarse.t0() + 0.05; t < coarse.t_end(); t += 1.0)
        worst = std::fmax(worst, std::abs(coarse.at(t).lambda - rg.at(t).lambda));
    check(worst < 1e-9, "cubic Hermite dense output accurate to < 1e-9 at mid-nodes (got " + std::to_string(worst) + ")");
    check(coarse.derivative_at(coarse.t0() - 1.0).lambda == 0.0, "derivative is zero outside the integration range");

    // Very heavy top: y_t becomes non-perturbative below M_Pl.
    Inputs heavy;
    heavy.Mt = 240.0;
    heavy.matching = MatchingScheme::Extended;
    const RGSolution rh(heavy);
    check(!rh.perturbative(), "M_t = 240 GeV is flagged non-perturbative");
    check(rh.t_end() < rh.t_planck(), "non-perturbative run stops before M_Pl");
    return smvac_test::finish();
}
