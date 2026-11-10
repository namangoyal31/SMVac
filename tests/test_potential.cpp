// Tests of the effective potential.
//  1. lambda_cw1 equals the Coleman-Weinberg sum sum_i n_i M_i^4 (ln M_i^2/mu^2 - c_i)/(64 pi^2) divided by h^4/4 at
//     mu = h, with M_t^2 = y_t^2 h^2/2, M_W^2 = g2^2 h^2/4, M_Z^2 = (g2^2 + gY^2) h^2/4 and gY^2 = 3/5 g1^2.
//  2. The Chebyshev representation reproduces the potential it was built from; frozen outside the range.
//  3. The phi^6 deformation is opt-in and changes lambda by exactly (2/3) c6 phi^2/Lambda^2.
#include "smvac/constants.hpp"
#include "smvac/potential.hpp"
#include "test_util.hpp"

#include <cmath>

using namespace smvac;
using smvac_test::check;
using smvac_test::close;

int main() {
    const Couplings c{0.58, 0.55, 0.62, 0.52, 0.01, 0.008, -0.011};
    const double h = 3.0e15, mu = h;
    const double gy2 = 0.6 * c.g1 * c.g1;
    const double mt2 = c.yt * c.yt * h * h / 2, mw2 = c.g2 * c.g2 * h * h / 4, mz2 = (c.g2 * c.g2 + gy2) * h * h / 4;
    const double V1 = (-12 * mt2 * mt2 * (std::log(mt2 / (mu * mu)) - 1.5) + 6 * mw2 * mw2 * (std::log(mw2 / (mu * mu)) - 5.0 / 6.0) +
                       3 * mz2 * mz2 * (std::log(mz2 / (mu * mu)) - 5.0 / 6.0)) / (64 * kPi2);
    close(lambda_cw1(c), c.lambda + V1 / (h * h * h * h / 4), 1e-12, 0, "lambda_cw1 equals the CW sum at mu = h");

    // Legacy defect check: using g1^2 (GUT) instead of gY^2 in M_Z would give a different value.
    const double mz2_wrong = (c.g2 * c.g2 + c.g1 * c.g1) * h * h / 4;
    const double V1w = (-12 * mt2 * mt2 * (std::log(mt2 / (mu * mu)) - 1.5) + 6 * mw2 * mw2 * (std::log(mw2 / (mu * mu)) - 5.0 / 6.0) +
                        3 * mz2_wrong * mz2_wrong * (std::log(mz2_wrong / (mu * mu)) - 5.0 / 6.0)) / (64 * kPi2);
    check(std::abs(lambda_cw1(c) - (c.lambda + V1w / (h * h * h * h / 4))) > 1e-5, "Z mass uses gY^2, not GUT g1^2");

    const FrozenChebyshevPotential p([](double x) { return -0.01 + 1e-4 * std::sin(0.3 * x) + 2e-5 * (x - 30) * (x - 30) / 50; },
                                     10.0, 44.0);
    close(p.lam(27.3), -0.01 + 1e-4 * std::sin(0.3 * 27.3) + 2e-5 * (27.3 - 30) * (27.3 - 30) / 50, 0, 1e-15, "Chebyshev value");
    close(p.dlam(27.3), 0.3e-4 * std::cos(0.3 * 27.3) + 4e-5 * (27.3 - 30) / 50, 0, 1e-14, "Chebyshev first derivative");
    close(p.d2lam(27.3), -0.09e-4 * std::sin(0.3 * 27.3) + 4e-5 / 50, 0, 1e-13, "Chebyshev second derivative");
    close(p.lam(50.0), p.lam(44.0), 0, 0, "frozen above x_hi");
    check(p.dlam(50.0) == 0.0 && p.d2lam(3.0) == 0.0, "derivatives vanish in frozen regions");

    const ConstantPotential base(-0.01);
    const Phi6Deformation d(base, 1.0, kPlanckMass);
    const double x = std::log(1e18);
    close(d.lam(x) - base.lam(x), (2.0 / 3.0) * std::exp(2 * x) / (kPlanckMass * kPlanckMass), 1e-13, 0,
          "phi^6 deformation shifts lambda by (2/3) phi^2 / Lambda^2");
    return smvac_test::finish();
}
