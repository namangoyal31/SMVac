// Tests of the numerical bounce solver (Chebyshev collocation).
// Reference values are the consensus of the independent Python implementations M1 (shooting, physical variables),
// M2 (shooting, scale-covariant variables), M3 (adaptive collocation) and M4 (spectral collocation), which agree
// to <= 1e-11 relative (validation/results/method_comparison.csv). They are NOT outputs of this C++ code.
#include "smvac/analytic.hpp"
#include "smvac/bounce.hpp"
#include "smvac/constants.hpp"
#include "smvac/rge.hpp"
#include "test_util.hpp"

#include <cmath>

using namespace smvac;
using smvac_test::check;
using smvac_test::close;

namespace {
void run(const Potential& p, double lo, double hi, double ref, double tol, const char* name) {
    const AnalyticEstimates a = analytic_estimates(p, lo, hi);
    const BounceResult b = solve_bounce(p, a);
    check(b.ok, std::string(name) + ": solver converged (" + b.message + ")");
    close(b.S, ref, tol, 0, std::string(name) + ": S vs independent multi-method reference");
    check(std::abs(b.derrick) < 1e-9 && std::abs(b.dilatation) < 1e-9, std::string(name) + ": identities < 1e-9");
    check(a.S_FL <= b.S * (1 + 1e-12), std::string(name) + ": rigorous lower bound S_FL <= S");
    check(b.dS < 1e-6 * b.S, std::string(name) + ": self-convergence estimate dS/S < 1e-6");
}
}  // namespace

int main() {
    // Constant lambda: every Fubini-Lipatov profile is a solution, S = 8 pi^2 / (3 |lambda|).
    {
        const ConstantPotential c(-0.0125);
        const AnalyticEstimates a = analytic_estimates(c, 25, 55);
        const BounceResult b = solve_bounce(c, a);
        close(b.S, kEightPi2Over3 / 0.0125, 1e-10, 0, "constant lambda: Fubini-Lipatov action reproduced");
    }
    run(ParabolaPotential(-0.0125, 5e-5, 39.0), 25, 55, 2107.0357973921, 1e-10, "parabola q=250 (thick wall)");
    run(ParabolaPotential(-2.5 * 5e-5, 5e-5, 39.0), 25, 55, 225400.864807, 1e-9, "parabola q=2.5 (thin-wall onset)");
    {
        const RGSolution rg(Inputs{});
        const auto pot = make_sm_potential(rg, Prescription::CW1);
        const AnalyticEstimates a = analytic_estimates(*pot, pot->x_lo(), pot->x_hi());
        const BounceResult b = solve_bounce(*pot, a);
        check(b.ok, "SM central point: solver converged");
        close(b.S, 1856.01513783, 2e-10, 0, "SM central point (cw1): S vs independent multi-method reference");
        check(b.S <= a.S_CGM, "SM central point: S below the variational upper estimate");
    }
    return smvac_test::finish();
}
