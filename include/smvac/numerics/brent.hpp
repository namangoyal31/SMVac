// SMVac — Brent's methods for 1-D root finding and minimisation (R.P. Brent, "Algorithms for Minimization
// without Derivatives", 1973). Header-only.
#pragma once

#include <cmath>
#include <stdexcept>

namespace smvac {
namespace numerics {

struct MinResult {
    double x, fx;
    int iterations;
};

template <class F>
double brent_root(F f, double a, double b, double xtol = 1e-14, int max_iter = 200) {
    double fa = f(a), fb = f(b);
    if (fa == 0.0) return a;
    if (fb == 0.0) return b;
    if ((fa > 0) == (fb > 0)) throw std::runtime_error("brent_root: root not bracketed");
    double c = a, fc = fa, d = b - a, e = d;
    for (int it = 0; it < max_iter; ++it) {
        if ((fb > 0) == (fc > 0)) { c = a; fc = fa; d = e = b - a; }
        if (std::abs(fc) < std::abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
        const double tol = 2.0 * 2.2e-16 * std::abs(b) + 0.5 * xtol, m = 0.5 * (c - b);
        if (std::abs(m) <= tol || fb == 0.0) return b;
        if (std::abs(e) >= tol && std::abs(fa) > std::abs(fb)) {
            double p, q, r, s = fb / fa;
            if (a == c) { p = 2.0 * m * s; q = 1.0 - s; }
            else {
                q = fa / fc; r = fb / fc;
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0) q = -q; else p = -p;
            if (2.0 * p < std::fmin(3.0 * m * q - std::abs(tol * q), std::abs(e * q))) { e = d; d = p / q; }
            else { d = m; e = m; }
        } else { d = m; e = m; }
        a = b; fa = fb;
        b += (std::abs(d) > tol) ? d : (m > 0 ? tol : -tol);
        fb = f(b);
    }
    throw std::runtime_error("brent_root: no convergence");
}

template <class F>
MinResult brent_minimize(F f, double a, double b, double xtol = 1e-10, int max_iter = 500) {
    const double cg = 0.3819660112501051;
    double x = a + cg * (b - a), w = x, v = x, fx = f(x), fw = fx, fv = fx, d = 0.0, e = 0.0;
    int it = 0;
    for (; it < max_iter; ++it) {
        const double xm = 0.5 * (a + b), tol1 = xtol + 1e-15 * std::abs(x), tol2 = 2.0 * tol1;
        if (std::abs(x - xm) <= tol2 - 0.5 * (b - a)) break;
        bool golden = true;
        if (std::abs(e) > tol1) {
            const double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw), p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0) p = -p; else q = -q;
            const double etemp = e;
            e = d;
            if (std::abs(p) < std::abs(0.5 * q * etemp) && p > q * (a - x) && p < q * (b - x)) {
                d = p / q;
                const double u = x + d;
                if (u - a < tol2 || b - u < tol2) d = (xm - x >= 0) ? tol1 : -tol1;
                golden = false;
            }
        }
        if (golden) { e = (x >= xm) ? a - x : b - x; d = cg * e; }
        const double u = (std::abs(d) >= tol1) ? x + d : x + (d > 0 ? tol1 : -tol1);
        const double fu = f(u);
        if (fu <= fx) {
            if (u >= x) a = x; else b = x;
            v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
            else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
        }
    }
    return MinResult{x, fx, it};
}

}  // namespace numerics
}  // namespace smvac
