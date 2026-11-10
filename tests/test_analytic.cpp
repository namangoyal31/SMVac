// Tests of the analytical approximations and bounds.
//  1. Constant lambda: S_FL = S_CGM = 8 pi^2/(3|lambda|) (exact Fubini-Lipatov result).
//  2. Parabolic lambda: S_curv = S_CGM = 8 pi^2 / (3(|l0| - 13 l2/72)) exactly (the CGM estimate of a quadratic
//     lambda with the FL trial profile is the closed-form curvature correction).
//  3. Ordering S_FL <= S_CGM.
//  4. SM central point: agreement with the independent Python reference implementation.
#include "smvac/analytic.hpp"
#include "smvac/constants.hpp"
#include "smvac/rge.hpp"
#include "test_util.hpp"

#include <cmath>

using namespace smvac;
using smvac_test::check;
using smvac_test::close;

int main() {
    const ConstantPotential c(-0.0125);
    const AnalyticEstimates a = analytic_estimates(c, 25.0, 55.0);
    close(a.S_FL, kEightPi2Over3 / 0.0125, 1e-14, 0, "constant lambda: S_FL exact");
    close(a.S_CGM, kEightPi2Over3 / 0.0125, 1e-12, 0, "constant lambda: S_CGM exact");

    const double l0 = -0.0125, l2 = 5e-5;
    const ParabolaPotential p(l0, l2, 39.0);
    const AnalyticEstimates b = analytic_estimates(p, 25.0, 55.0);
    const double closed = kEightPi2Over3 / (std::abs(l0) - kCurvatureCoefficient * l2);
    close(b.S_curv, closed, 1e-13, 0, "parabola: S_curv closed form");
    close(b.S_CGM, closed, 1e-11, 0, "parabola: S_CGM equals the closed-form curvature correction");
    // Locating a minimum is limited to |dx| ~ sqrt(2 eps |lambda| / lambda'') ~ 3e-7 in double precision.
    close(b.x_min, 39.0, 0, 1e-6, "parabola: location of the minimum");
    close(b.lnh_CGM, 39.0 + 5.0 / 6.0, 0, 1e-6, "parabola: optimal FL amplitude ln h = x0 + 5/6");
    check(b.S_FL <= b.S_CGM, "S_FL <= S_CGM");

    // SM central point, prescription cw1; reference values from python/smvac/reference (independent code).
    const RGSolution rg(Inputs{});
    const auto pot = make_sm_potential(rg, Prescription::CW1);
    const AnalyticEstimates s = analytic_estimates(*pot, pot->x_lo(), pot->x_hi());
    close(s.lambda_min, -0.0141904768, 0, 2e-10, "SM cw1 lambda_min vs Python reference");
    close(s.S_FL, 1854.6906782417, 1e-10, 0, "SM cw1 S_FL vs Python reference");
    close(s.S_CGM, 1856.0164919353, 1e-10, 0, "SM cw1 S_CGM vs Python reference");
    check(!s.min_at_boundary && s.CGM_interior, "SM central point: interior minimum and interior CGM optimum");
    return smvac_test::finish();
}
