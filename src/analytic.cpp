// SMVac — analytical approximations and bounds.
#include "smvac/analytic.hpp"
#include "smvac/constants.hpp"
#include "smvac/numerics/brent.hpp"
#include "smvac/numerics/quadrature.hpp"

#include <cmath>
#include <limits>

namespace smvac {

double lambda_bar(const Potential& pot, double lnh) {
    auto f = [&](double y) { return 6.0 * std::exp(2.0 * y) * (1.0 - std::exp(y)) * pot.lam(lnh + y); };
    std::vector<double> breaks;
    for (double k : pot.kinks()) breaks.push_back(k - lnh);
    return numerics::integrate_split(f, -45.0, 0.0, breaks, 1e-16, 1e-13);
}

void global_minimum(const Potential& pot, double x_lo, double x_hi, double& lam_min, double& x_min, bool& at_boundary) {
    const int n = 4000;
    int ibest = 0;
    double best = std::numeric_limits<double>::infinity();
    for (int i = 0; i <= n; ++i) {
        const double v = pot.lam(x_lo + (x_hi - x_lo) * i / n);
        if (v < best) { best = v; ibest = i; }
    }
    if (ibest == 0 || ibest == n) {
        lam_min = best;
        x_min = (ibest == 0) ? x_lo : x_hi;
        at_boundary = true;
        return;
    }
    const double h = (x_hi - x_lo) / n;
    const auto r = numerics::brent_minimize([&](double x) { return pot.lam(x); }, x_lo + (ibest - 1) * h,
                                            x_lo + (ibest + 1) * h, 1e-12);
    lam_min = r.fx;
    x_min = r.x;
    at_boundary = false;
}

AnalyticEstimates analytic_estimates(const Potential& pot, double x_lo, double x_hi) {
    AnalyticEstimates a;
    global_minimum(pot, x_lo, x_hi, a.lambda_min, a.x_min, a.min_at_boundary);
    const double inf = std::numeric_limits<double>::infinity();
    a.S_FL = a.S_curv = a.S_CGM = inf;
    if (a.lambda_min >= 0.0) {
        a.notes.push_back("absolutely stable: lambda_min >= 0");
        return a;
    }
    a.S_FL = kEightPi2Over3 / std::abs(a.lambda_min);
    a.lnLambdaB_literature = a.x_min;
    if (!a.min_at_boundary) {
        a.lambda2 = pot.d2lam(a.x_min);
        const double den = std::abs(a.lambda_min) - kCurvatureCoefficient * a.lambda2;
        a.S_curv = den > 0 ? kEightPi2Over3 / den : inf;
        if (den <= 0) a.notes.push_back("curvature correction exceeds |lambda_min|: first-order expansion invalid");
    } else {
        a.S_curv = std::numeric_limits<double>::quiet_NaN();
        a.notes.push_back("lambda minimum on a frozen boundary (Planck-sensitive / degenerate)");
    }
    // CGM variational estimate: scan ln h, then Brent around the best grid point.
    const double lo = a.x_min - 6.0, hi = a.x_min + 10.0;
    const int n = 160;
    int ib = 0;
    double best = inf;
    for (int i = 0; i <= n; ++i) {
        const double v = lambda_bar(pot, lo + (hi - lo) * i / n);
        if (v < best) { best = v; ib = i; }
    }
    double lnh = lo + (hi - lo) * ib / n, lb = best;
    if (ib > 0 && ib < n) {
        const double h = (hi - lo) / n;
        const auto r = numerics::brent_minimize([&](double g) { return lambda_bar(pot, g); }, lnh - h, lnh + h, 1e-9);
        if (r.fx <= lb) { lb = r.fx; lnh = r.x; }
        a.CGM_interior = true;
    } else {
        a.notes.push_back("CGM optimum at the edge of the search range (no interior optimum; degenerate)");
    }
    a.lnh_CGM = lnh;
    if (lb < 0) {
        a.S_CGM = kEightPi2Over3 / std::abs(lb);
        a.lnLambdaB_FL = lnh + 0.5 * std::log(std::abs(lb) / 8.0);
    }
    return a;
}

}  // namespace smvac
