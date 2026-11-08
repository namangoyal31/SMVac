// SMVac — numerical reference solution of the O(4) bounce equation for V = lambda(phi) phi^4 / 4.
//
// Variables tau = ln r, f = r phi, x = ln phi = ln f - tau. Bounce equation
//     f'' = f + f^3 G(x),   G = lambda + (1/4) d lambda/dx,
// with f ~ e^{+tau} as r -> 0 and f ~ e^{-tau} as r -> infinity; action S = T + V,
//     T = 2 pi^2 int (f' - f)^2 / 2 dtau,   V = 2 pi^2 int lambda f^4 / 4 dtau.
// A solution satisfies the Derrick identity T + 2V = 0 and the dilatation identity int beta f^4 dtau = 0.
// The result is a NUMERICAL solution with an estimated uncertainty, never an "exact" value.
// The method was selected by the comparison in docs/NUMERICAL_METHOD_COMPARISON.md.
#pragma once

#include "smvac/analytic.hpp"
#include "smvac/potential.hpp"

#include <string>

namespace smvac {

struct BounceOptions {
    int nodes = 256;          // initial number of Chebyshev collocation nodes
    int nodes_check = 384;    // second resolution used for the self-convergence estimate (scaled with nodes)
    int max_nodes = 1024;     // resolution is doubled (pairs 256/384, 512/768) until the identities hold
    double identity_tol = 1e-9;
    // Reduced-accuracy tier: if identity_tol cannot be met (e.g. the bounce crosses the frozen-coupling kink at
    // M_Pl, where spectral convergence becomes algebraic), accept identities below this value and flag the result.
    double identity_tol_reduced = 1e-5;
    double half_width = 15.0; // half-width L of the tau interval (Robin conditions are exact up to O(f^3) there)
    double newton_tol = 1e-14;
    int max_newton = 60;
};

struct BounceResult {
    bool ok = false;
    std::string method = "chebyshev-collocation";
    std::string message;
    double S = 0, T = 0, V = 0;
    double derrick = 0;       // (T + 2V) / T
    double dilatation = 0;    // int beta f^4 / int |lambda| f^4
    double lnphi0 = 0;        // ln phi(r = 0)
    double lnR = 0;           // tau at the maximum of f = r phi (the radius R for a Fubini-Lipatov profile)
    double dS = 0;            // estimated numerical uncertainty of S
    bool reduced_accuracy = false;   // accepted in the reduced-accuracy tier (see BounceOptions)
    int newton_iterations = 0;
};

BounceResult solve_bounce(const Potential& pot, const AnalyticEstimates& guess, const BounceOptions& opt = BounceOptions());

}  // namespace smvac
