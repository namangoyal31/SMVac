// SMVac — analytical approximations and bounds for the tree-level O(4) bounce action of V = lambda(phi) phi^4/4.
//
//   S_FL   = 8 pi^2 / (3 |lambda_min|)                          literature estimate; RIGOROUS LOWER BOUND on every
//                                                               finite-action bounce (Pohozaev + sharp Sobolev)
//   S_curv = 8 pi^2 / (3 (|lambda_min| - (13/72) lambda''))     analytical approximation (first-order perturbation
//                                                               theory about the FL instanton, quadratic lambda)
//   S_CGM  = 8 pi^2 / (3 |min_h lambdabar(h)|)                  variational UPPER ESTIMATE (a rigorous bound only if
//                                                               the Coleman-Glaser-Martin infimum is attained)
//   lambdabar(h) = 6 int_0^1 w (1 - w) lambda(h w) dw,  x = ln(phi/GeV).
// See docs/SCIENTIFIC_METHODOLOGY.md section 6 for the derivations.
#pragma once

#include "smvac/potential.hpp"

#include <string>
#include <vector>

namespace smvac {

constexpr double kCurvatureCoefficient = 13.0 / 72.0;   // (1/2) Var(ln w) under the FL weight 6 w (1 - w)

struct AnalyticEstimates {
    double lambda_min = 0, x_min = 0;
    bool min_at_boundary = false;      // minimum on a frozen boundary (Planck-sensitive / degenerate)
    double lambda2 = 0;                // d^2 lambda/dx^2 at the minimum
    double S_FL = 0, S_curv = 0, S_CGM = 0;
    double lnh_CGM = 0;                // ln phi(0) of the optimal FL trial profile
    bool CGM_interior = false;
    double lnLambdaB_literature = 0;   // Buttazzo identification: Lambda_B = argmin lambda
    double lnLambdaB_FL = 0;           // 1/R of the optimal FL trial instanton
    std::vector<std::string> notes;
};

// lambdabar(h) for a given potential (adaptive Gauss-Kronrod in y = ln w, split at the potential's kinks).
double lambda_bar(const Potential& pot, double lnh);

// Global minimum of lambda over [x_lo, x_hi] (dense scan + Brent); frozen potentials are constant outside.
void global_minimum(const Potential& pot, double x_lo, double x_hi, double& lam_min, double& x_min, bool& at_boundary);

AnalyticEstimates analytic_estimates(const Potential& pot, double x_lo, double x_hi);

}  // namespace smvac
