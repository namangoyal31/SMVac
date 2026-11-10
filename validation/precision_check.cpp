// Floating-point precision check of the Chebyshev-collocation bounce solver.
//
// The production solver (src/bounce.cpp) works in IEEE double. This program re-implements the same discretisation
// as a template on the scalar type and solves the same problems in double and in long double (64-bit mantissa on
// x86 GCC, eps ~ 1.1e-19). With the discretisation held fixed (same N, same domain), the difference
// S(long double) - S(double) measures the round-off error of the double-precision solution. The potentials are
// evaluated natively in the wider type: closed form for the parabolic model; Clenshaw recurrence in long double on
// the (double) Chebyshev coefficients for the SM potential, which DEFINES the potential for this test.
//
// Build: cmake --build build --target precision_check
//        (or: g++ -std=c++14 -O2 -Iinclude validation/precision_check.cpp src/*.cpp -o precision_check)
// Output: one line per case and resolution (printed through double, which suffices for the differences).
#include "smvac/analytic.hpp"
#include "smvac/potential.hpp"
#include "smvac/rge.hpp"

#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

using namespace smvac;

template <class R>
struct Pot {   // lambda(x), lambda'(x), lambda''(x) in type R
    std::function<void(R, R&, R&, R&)> eval;
};

template <class R>
Pot<R> parabola(R l0, R l2, R x0) {
    return Pot<R>{[=](R x, R& l, R& d, R& dd) { l = l0 + l2 * (x - x0) * (x - x0) / 2; d = l2 * (x - x0); dd = l2; }};
}

// Clenshaw for a Chebyshev series with derivative series computed in R.
template <class R>
struct ChebR {
    R a, b;
    std::vector<R> c, c1, c2;
    static std::vector<R> deriv(const std::vector<R>& cc, R a, R b) {
        const int n = static_cast<int>(cc.size()) - 1;
        std::vector<R> d(n + 1, 0);
        for (int j = n - 1; j >= 0; --j) d[j] = (j + 2 <= n ? d[j + 2] : R(0)) + 2 * R(j + 1) * cc[j + 1];
        d[0] /= 2;
        for (R& v : d) v *= 2 / (b - a);
        d.resize(n);
        return d;
    }
    static R clenshaw(const std::vector<R>& cc, R u) {
        R b1 = 0, b2 = 0;
        for (int j = static_cast<int>(cc.size()) - 1; j >= 1; --j) {
            const R t = 2 * u * b1 - b2 + cc[j];
            b2 = b1;
            b1 = t;
        }
        return u * b1 - b2 + cc[0];
    }
};

template <class R>
Pot<R> sm_pot(const FrozenChebyshevPotential& p) {
    auto s = std::make_shared<ChebR<R>>();
    s->a = p.x_lo();
    s->b = p.x_hi();
    for (double v : p.series().coefficients()) s->c.push_back(R(v));
    s->c1 = ChebR<R>::deriv(s->c, s->a, s->b);
    s->c2 = ChebR<R>::deriv(s->c1, s->a, s->b);
    return Pot<R>{[s](R x, R& l, R& d, R& dd) {
        const R xc = x < s->a ? s->a : (x > s->b ? s->b : x);
        const R u = (2 * xc - s->a - s->b) / (s->b - s->a);
        l = ChebR<R>::clenshaw(s->c, u);
        const bool inside = x > s->a && x < s->b;
        d = inside ? ChebR<R>::clenshaw(s->c1, u) : R(0);
        dd = inside ? ChebR<R>::clenshaw(s->c2, u) : R(0);
    }};
}

template <class R>
bool solve_lu(std::vector<R> A, std::vector<R>& b, int n) {
    for (int k = 0; k < n; ++k) {
        int p = k;
        for (int i = k + 1; i < n; ++i)
            if (std::fabs(A[i * n + k]) > std::fabs(A[p * n + k])) p = i;
        if (A[p * n + k] == 0) return false;
        if (p != k) {
            for (int j = 0; j < n; ++j) std::swap(A[k * n + j], A[p * n + j]);
            std::swap(b[k], b[p]);
        }
        for (int i = k + 1; i < n; ++i) {
            const R m = A[i * n + k] / A[k * n + k];
            for (int j = k; j < n; ++j) A[i * n + j] -= m * A[k * n + j];
            b[i] -= m * b[k];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        R s = b[i];
        for (int j = i + 1; j < n; ++j) s -= A[i * n + j] * b[j];
        b[i] = s / A[i * n + i];
    }
    return true;
}

struct Out { double S, derrick, dilatation; int it; bool ok; };

template <class R>
Out collocate(const Pot<R>& pot, R lam_min, R lnh, int N, R L) {
    const R pi = std::acos(R(-1));
    const int n = N + 1;
    std::vector<R> x(n), c(n, 1), D(n * n, 0), D1, D2(n * n, 0);
    for (int j = 0; j < n; ++j) x[j] = std::cos(pi * j / N);
    c[0] = c[N] = 2;
    for (int j = 0; j < n; ++j) if (j % 2) c[j] = -c[j];
    for (int i = 0; i < n; ++i) {
        R s = 0;
        for (int j = 0; j < n; ++j) if (i != j) { D[i * n + j] = (c[i] / c[j]) / (x[i] - x[j]); s += D[i * n + j]; }
        D[i * n + i] = -s;
    }
    D1 = D;
    for (R& v : D1) v /= L;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < n; ++j) D2[i * n + j] += D1[i * n + k] * D1[k * n + j];
    const R A = std::sqrt(2 / std::fabs(lam_min)), tc = std::log(2 * A) - lnh;
    std::vector<R> tau(n), f(n);
    for (int j = 0; j < n; ++j) { tau[j] = tc + L * x[j]; f[j] = A / std::cosh(tau[j] - tc); }
    auto resid = [&](const std::vector<R>& g, std::vector<R>& Rv, std::vector<R>& Gv, std::vector<R>& Gx) {
        Rv.assign(n, 0); Gv.resize(n); Gx.resize(n);
        for (int i = 0; i < n; ++i) {
            R l, d, dd;
            pot.eval(std::log(g[i] > 0 ? g[i] : R(1e-300)) - tau[i], l, d, dd);
            Gv[i] = l + d / 4; Gx[i] = d + dd / 4;
            R s = 0;
            for (int j = 0; j < n; ++j) s += D2[i * n + j] * g[j];
            Rv[i] = s - g[i] - g[i] * g[i] * g[i] * Gv[i];
        }
        R d0 = 0, dN = 0;
        for (int j = 0; j < n; ++j) { d0 += D1[j] * g[j]; dN += D1[N * n + j] * g[j]; }
        Rv[0] = d0 + g[0];
        Rv[N] = dN - g[N];
    };
    std::vector<R> Rv, Gv, Gx, Rn, Gn, Gxn;
    resid(f, Rv, Gv, Gx);
    int it = 0;
    for (it = 1; it <= 60; ++it) {
        std::vector<R> J(D2);
        for (int i = 1; i < N; ++i) J[i * n + i] -= 1 + 3 * f[i] * f[i] * Gv[i] + f[i] * f[i] * Gx[i];
        for (int j = 0; j < n; ++j) { J[j] = D1[j] + (j == 0); J[N * n + j] = D1[N * n + j] - (j == N); }
        std::vector<R> dlt(n);
        for (int j = 0; j < n; ++j) dlt[j] = -Rv[j];
        if (!solve_lu(J, dlt, n)) return Out{0, 0, 0, it, false};
        R dm = 0, fm = 0;
        for (int j = 0; j < n; ++j) { f[j] += dlt[j]; dm = std::fmax(dm, std::fabs(dlt[j])); fm = std::fmax(fm, std::fabs(f[j])); }
        resid(f, Rv, Gv, Gx);
        if (dm < R(1e-17) * fm) break;   // near the round-off level of the wider type (and below that of double)
    }
    // Clenshaw-Curtis weights
    std::vector<R> w(n, 0), v(N - 1, 1);
    w[0] = w[N] = 1 / (R(N) * N - 1);
    for (int k = 1; k < N / 2; ++k)
        for (int j = 1; j < N; ++j) v[j - 1] -= 2 * std::cos(2 * k * pi * j / N) / (4 * R(k) * k - 1);
    for (int j = 1; j < N; ++j) { v[j - 1] -= std::cos(N * pi * j / N) / (R(N) * N - 1); w[j] = 2 * v[j - 1] / N; }
    R T = 0, V = 0, Dl = 0, W = 0;
    for (int j = 0; j < n; ++j) {
        R fp = 0;
        for (int k = 0; k < n; ++k) fp += D1[j * n + k] * f[k];
        R l, d, dd;
        pot.eval(std::log(f[j]) - tau[j], l, d, dd);
        const R f4 = f[j] * f[j] * f[j] * f[j];
        T += w[j] * L * (fp - f[j]) * (fp - f[j]) / 2;
        V += w[j] * L * l * f4 / 4;
        Dl += w[j] * L * d * f4;
        W += w[j] * L * std::fabs(l) * f4;
    }
    const R tp2 = 2 * pi * pi;
    T = tp2 * (T + f[0] * f[0]);
    V *= tp2;
    return Out{static_cast<double>(T + V), static_cast<double>((T + 2 * V) / T), static_cast<double>(Dl / W), it, true};
}

int main() {
    struct Case { const char* name; double lam_min, lnh; std::function<Out(int, bool)> run; };
    const RGSolution rg(Inputs{});
    const auto smp = make_sm_potential(rg, Prescription::CW1);
    const AnalyticEstimates a_sm = analytic_estimates(*smp, smp->x_lo(), smp->x_hi());
    const ParabolaPotential pthick(-0.0125, 5e-5, 39.0), pthin(-2.5 * 5e-5, 5e-5, 39.0);
    const AnalyticEstimates a_thick = analytic_estimates(pthick, 25, 55), a_thin = analytic_estimates(pthin, 25, 55);
    std::vector<Case> cases = {
        {"parabola q=250", a_thick.lambda_min, a_thick.lnh_CGM, [&](int N, bool ld) {
             return ld ? collocate<long double>(parabola<long double>(-0.0125L, 5e-5L, 39.0L), -0.0125L, a_thick.lnh_CGM, N, 15.0L)
                       : collocate<double>(parabola<double>(-0.0125, 5e-5, 39.0), -0.0125, a_thick.lnh_CGM, N, 15.0); }},
        {"parabola q=2.5", a_thin.lambda_min, a_thin.lnh_CGM, [&](int N, bool ld) {
             return ld ? collocate<long double>(parabola<long double>(-1.25e-4L, 5e-5L, 39.0L), -1.25e-4L, a_thin.lnh_CGM, N, 15.0L)
                       : collocate<double>(parabola<double>(-1.25e-4, 5e-5, 39.0), -1.25e-4, a_thin.lnh_CGM, N, 15.0); }},
        {"SM cw1 central", a_sm.lambda_min, a_sm.lnh_CGM, [&](int N, bool ld) {
             return ld ? collocate<long double>(sm_pot<long double>(*smp), a_sm.lambda_min, a_sm.lnh_CGM, N, 15.0L)
                       : collocate<double>(sm_pot<double>(*smp), a_sm.lambda_min, a_sm.lnh_CGM, N, 15.0); }},
    };
    std::printf("case, N, S_double, S_longdouble - S_double, rel_diff, derrick_double, derrick_longdouble\n");
    for (auto& c : cases)
        for (int N : {192, 256, 384}) {
            const Out d = c.run(N, false), l = c.run(N, true);
            std::printf("%s, %d, %.15g, %.3e, %.3e, %.2e, %.2e\n", c.name, N, d.S, l.S - d.S, (l.S - d.S) / d.S,
                        d.derrick, l.derrick);
        }
    return 0;
}
