// SMVac — adaptive Gauss-Kronrod (7-15) quadrature with absolute/relative tolerance. Header-only.
#pragma once

#include <cmath>
#include <algorithm>
#include <vector>

namespace smvac {
namespace numerics {

namespace detail {
constexpr double kXgk[8] = {0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
                            0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
                            0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
                            0.207784955007898467600689403773245, 0.000000000000000000000000000000000};
constexpr double kWgk[8] = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
                            0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
                            0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
                            0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
constexpr double kWg[4] = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
                           0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

template <class F>
void gk15(F& f, double a, double b, double& res, double& err) {
    const double c = 0.5 * (a + b), h = 0.5 * (b - a);
    const double fc = f(c);
    double rk = fc * kWgk[7], rg = fc * kWg[3];
    for (int j = 0; j < 7; ++j) {
        const double x = h * kXgk[j];
        const double s = f(c - x) + f(c + x);
        rk += kWgk[j] * s;
        if (j % 2 == 1) rg += kWg[j / 2] * s;
    }
    res = rk * h;
    err = std::abs((rk - rg) * h);
}

template <class F>
double adapt(F& f, double a, double b, double whole, double err, double eps_abs, double eps_rel, int depth) {
    if (err <= std::fmax(eps_abs, eps_rel * std::abs(whole)) || depth > 50) return whole;
    const double m = 0.5 * (a + b);
    double r1, e1, r2, e2;
    gk15(f, a, m, r1, e1);
    gk15(f, m, b, r2, e2);
    return adapt(f, a, m, r1, e1, 0.5 * eps_abs, eps_rel, depth + 1) +
           adapt(f, m, b, r2, e2, 0.5 * eps_abs, eps_rel, depth + 1);
}
}  // namespace detail

template <class F>
double integrate(F f, double a, double b, double eps_abs = 1e-15, double eps_rel = 1e-13) {
    double r, e;
    detail::gk15(f, a, b, r, e);
    return detail::adapt(f, a, b, r, e, eps_abs, eps_rel, 0);
}

// Integrate over [a, b] splitting at the given interior break points (non-smooth points of the integrand).
template <class F>
double integrate_split(F f, double a, double b, std::vector<double> breaks, double eps_abs = 1e-15,
                       double eps_rel = 1e-13) {
    std::sort(breaks.begin(), breaks.end());
    std::vector<double> e{a};
    for (double x : breaks)
        if (x > a && x < b) e.push_back(x);
    e.push_back(b);
    double s = 0.0;
    for (std::size_t i = 0; i + 1 < e.size(); ++i) s += integrate(f, e[i], e[i + 1], eps_abs, eps_rel);
    return s;
}

}  // namespace numerics
}  // namespace smvac
