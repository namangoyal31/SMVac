// SMVac — Standard Model beta functions (Buttazzo et al. 2013, App. B), t = ln(mu^2).
#pragma once

#include "smvac/couplings.hpp"

namespace smvac {

// d(x^2)/dt for the gauge and Yukawa couplings and d(lambda)/dt.
struct BetaSquares {
    double g1sq = 0, g2sq = 0, g3sq = 0, ytsq = 0, ybsq = 0, ytausq = 0, lambda = 0;
};

BetaSquares beta_squares(const Couplings& c);
// Same with an explicit loop-counting factor (1/(16 pi^2) physically); used by the power-counting test.
BetaSquares beta_squares(const Couplings& c, double loop_factor);

// dc/dt for the coupling vector itself (dx/dt = (dx^2/dt) / (2x)).
Couplings coupling_derivatives(const Couplings& c);

}  // namespace smvac
