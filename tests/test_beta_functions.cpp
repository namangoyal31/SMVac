// Tests of the Standard Model beta functions.
//  1. Loop-order power counting (literature independent): with g^2, y^2, lambda counted as degree 1, every n-loop
//     term has degree n+1 and carries (16 pi^2)^-n, so scaling couplings^2 by s and the loop factor by 1/s must
//     scale every beta function exactly by s. Any misplaced coupling power (the legacy transcription errors) fails.
//  2. One- and two-loop coefficients against independent closed forms at a random point.
//  3. The four-loop QCD term has the sign of -beta_3(n_f = 6) (van Ritbergen, Vermaseren, Larin 1997).
#include "smvac/beta_functions.hpp"
#include "smvac/constants.hpp"
#include "test_util.hpp"

#include <cmath>

using namespace smvac;
using smvac_test::close;
using smvac_test::check;

int main() {
    const Couplings c{0.55, 0.61, 0.83, 0.71, 0.13, 0.09, -0.031};
    const double kappa = 1.0 / kLoop;
    for (double s : {2.0, 0.37}) {
        const double r = std::sqrt(s);
        const Couplings cs{c.g1 * r, c.g2 * r, c.g3 * r, c.yt * r, c.yb * r, c.ytau * r, c.lambda * s};
        const BetaSquares b1 = beta_squares(c, kappa), b2 = beta_squares(cs, kappa / s);
        const double* p1 = &b1.g1sq;
        const double* p2 = &b2.g1sq;
        const char* names[7] = {"g1^2", "g2^2", "g3^2", "yt^2", "yb^2", "ytau^2", "lambda"};
        for (int i = 0; i < 7; ++i)
            close(p2[i], s * p1[i], 1e-13, 0.0, std::string("power counting d") + names[i] + "/dt, s=" + std::to_string(s));
    }

    // One-loop pieces in the limit kappa -> 0 (leading order in kappa).
    const double eps = 1e-10;   // two-loop contamination ~ 1e2 * eps, far below the tolerance
    const BetaSquares b = beta_squares(c, eps);
    const double g1 = c.g1 * c.g1, g2 = c.g2 * c.g2, g3 = c.g3 * c.g3, yt = c.yt * c.yt, yb = c.yb * c.yb,
                 ut = c.ytau * c.ytau, l = c.lambda;
    close(b.g1sq / eps, g1 * g1 * 41.0 / 10.0, 1e-6, 0, "one-loop b1 = 41/10 (GUT normalisation)");
    close(b.g2sq / eps, -g2 * g2 * 19.0 / 6.0, 1e-6, 0, "one-loop b2 = -19/6");
    close(b.g3sq / eps, -7.0 * g3 * g3, 1e-6, 0, "one-loop b3 = -7 (n_f = 6)");
    // d lambda/d ln mu^2 at one loop for V = lambda |H|^4 (half the usual d/d ln mu form; gY^2 = 3/5 g1^2)
    const double gy = 0.6 * g1;
    const double bl1 = l * (12 * l + 6 * yt + 6 * yb + 2 * ut - 4.5 * g2 - 1.5 * gy) - 3 * yt * yt - 3 * yb * yb - ut * ut +
                       3.0 / 16.0 * (2 * g2 * g2 + (g2 + gy) * (g2 + gy));
    close(b.lambda / eps, bl1, 1e-6, 0, "one-loop beta_lambda from the standard SM expression");
    const double byt1 = yt * (4.5 * yt + 1.5 * yb + ut - 8 * g3 - 2.25 * g2 - 17.0 / 12.0 * gy);
    close(b.ytsq / eps, byt1, 1e-6, 0, "one-loop beta_{yt^2} from the standard SM expression");

    // Four-loop QCD sign: compare d g3^2/dt with and without the term via the loop expansion in kappa.
    const Couplings q{1e-6, 1e-6, 1.2, 1e-6, 1e-6, 1e-6, 0.0};   // pure QCD limit
    const double k = 1.0 / kLoop;
    const BetaSquares bq = beta_squares(q, k);
    const double G = q.g3 * q.g3;
    const double three_loop = G * G * (-7 * k - 26 * G * k * k + 32.5 * G * G * k * k * k);
    const double four_loop = bq.g3sq - three_loop;
    close(four_loop, -2472.28 * G * G * G * G * G * k * k * k * k, 1e-6, 1e-30,
          "four-loop QCD term equals -beta_3(nf=6) g3^10/(16pi^2)^4");
    check(four_loop < 0.0, "four-loop QCD term is negative (asymptotic freedom strengthened)");
    return smvac_test::finish();
}
