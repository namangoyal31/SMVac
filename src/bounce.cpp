// SMVac — numerical reference bounce: Chebyshev spectral collocation with damped Newton iteration.
//
// Domain tau in [tau_c - L, tau_c + L] (Chebyshev-Gauss-Lobatto nodes); Robin conditions f' = +f (left) and
// f' = -f (right) are the exact linear asymptotics, so the truncation error is O(f^3) at the ends. The initial guess
// is the Fubini-Lipatov profile placed by the CGM variational estimate. The method choice and its validation are
// documented in docs/NUMERICAL_METHOD_COMPARISON.md.
#include "smvac/bounce.hpp"
#include "smvac/constants.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace smvac {
namespace {

using Mat = std::vector<double>;   // dense row-major (n x n)

// Chebyshev differentiation matrix on Lobatto points x_j = cos(pi j / N) (Trefethen, Spectral Methods in MATLAB).
void cheb_matrix(int N, Mat& D, std::vector<double>& x) {
    const int n = N + 1;
    x.resize(n);
    for (int j = 0; j < n; ++j) x[j] = std::cos(kPi * j / N);
    std::vector<double> c(n, 1.0);
    c[0] = c[N] = 2.0;
    for (int j = 0; j < n; ++j)
        if (j % 2) c[j] = -c[j];
    D.assign(static_cast<std::size_t>(n) * n, 0.0);
    for (int i = 0; i < n; ++i) {
        double rowsum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            const double v = (c[i] / c[j]) / (x[i] - x[j]);
            D[static_cast<std::size_t>(i) * n + j] = v;
            rowsum += v;
        }
        D[static_cast<std::size_t>(i) * n + i] = -rowsum;   // negative-sum trick for accuracy
    }
}

// Clenshaw-Curtis weights on the same Lobatto points (Trefethen's clencurt), for [-1, 1].
std::vector<double> clenshaw_curtis(int N) {
    std::vector<double> w(N + 1, 0.0), v(N - 1, 1.0);
    std::vector<double> th(N + 1);
    for (int j = 0; j <= N; ++j) th[j] = kPi * j / N;
    if (N % 2 == 0) {
        w[0] = w[N] = 1.0 / (static_cast<double>(N) * N - 1.0);
        for (int k = 1; k < N / 2; ++k)
            for (int j = 1; j < N; ++j) v[j - 1] -= 2.0 * std::cos(2.0 * k * th[j]) / (4.0 * k * k - 1.0);
        for (int j = 1; j < N; ++j) v[j - 1] -= std::cos(N * th[j]) / (static_cast<double>(N) * N - 1.0);
    } else {
        w[0] = w[N] = 1.0 / (static_cast<double>(N) * N);
        for (int k = 1; k <= (N - 1) / 2; ++k)
            for (int j = 1; j < N; ++j) v[j - 1] -= 2.0 * std::cos(2.0 * k * th[j]) / (4.0 * k * k - 1.0);
    }
    for (int j = 1; j < N; ++j) w[j] = 2.0 * v[j - 1] / N;
    return w;
}

// Solve A x = b in place by Gaussian elimination with partial pivoting. Returns false if singular.
bool lu_solve(Mat A, std::vector<double>& b, int n) {
    for (int k = 0; k < n; ++k) {
        int p = k;
        double amax = std::abs(A[static_cast<std::size_t>(k) * n + k]);
        for (int i = k + 1; i < n; ++i) {
            const double a = std::abs(A[static_cast<std::size_t>(i) * n + k]);
            if (a > amax) { amax = a; p = i; }
        }
        if (amax == 0.0) return false;
        if (p != k) {
            for (int j = 0; j < n; ++j) std::swap(A[static_cast<std::size_t>(k) * n + j], A[static_cast<std::size_t>(p) * n + j]);
            std::swap(b[k], b[p]);
        }
        const double piv = A[static_cast<std::size_t>(k) * n + k];
        for (int i = k + 1; i < n; ++i) {
            const double m = A[static_cast<std::size_t>(i) * n + k] / piv;
            if (m == 0.0) continue;
            for (int j = k; j < n; ++j) A[static_cast<std::size_t>(i) * n + j] -= m * A[static_cast<std::size_t>(k) * n + j];
            b[i] -= m * b[k];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < n; ++j) s -= A[static_cast<std::size_t>(i) * n + j] * b[j];
        b[i] = s / A[static_cast<std::size_t>(i) * n + i];
    }
    return true;
}

struct Solution {
    bool ok = false;
    std::string message;
    std::vector<double> tau, f, fp, w;
    int iterations = 0;
};

Solution collocate(const Potential& pot, double A, double tc, int N, double L, double tol, int max_it) {
    Solution s;
    Mat D;
    std::vector<double> xi;
    cheb_matrix(N, D, xi);
    const int n = N + 1;
    Mat D1(D), D2(static_cast<std::size_t>(n) * n, 0.0);
    for (double& v : D1) v /= L;
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k) {
            const double a = D1[static_cast<std::size_t>(i) * n + k];
            if (a == 0.0) continue;
            for (int j = 0; j < n; ++j) D2[static_cast<std::size_t>(i) * n + j] += a * D1[static_cast<std::size_t>(k) * n + j];
        }
    s.tau.resize(n);
    std::vector<double> f(n);
    for (int j = 0; j < n; ++j) {
        s.tau[j] = tc + L * xi[j];
        f[j] = A / std::cosh(s.tau[j] - tc);
    }
    auto residual = [&](const std::vector<double>& g, std::vector<double>& R, std::vector<double>& G) {
        R.assign(n, 0.0);
        G.resize(n);
        for (int i = 0; i < n; ++i) {
            const double x = std::log(std::fmax(g[i], 1e-300)) - s.tau[i];
            G[i] = pot.G(x);
            double d2 = 0.0;
            for (int j = 0; j < n; ++j) d2 += D2[static_cast<std::size_t>(i) * n + j] * g[j];
            R[i] = d2 - g[i] - g[i] * g[i] * g[i] * G[i];
        }
        double d10 = 0.0, d1N = 0.0;
        for (int j = 0; j < n; ++j) {
            d10 += D1[j] * g[j];
            d1N += D1[static_cast<std::size_t>(N) * n + j] * g[j];
        }
        R[0] = d10 + g[0];      // right end (xi = +1): f' = -f
        R[N] = d1N - g[N];      // left end (xi = -1):  f' = +f
    };
    auto norm = [](const std::vector<double>& v) {
        double s2 = 0.0;
        for (double a : v) s2 += a * a;
        return std::sqrt(s2);
    };
    std::vector<double> R, G, Rn, Gn;
    residual(f, R, G);
    auto rel_residual = [&]() {
        double rmax = 0.0, fm = 0.0;
        for (int j = 0; j < n; ++j) { rmax = std::fmax(rmax, std::abs(R[j])); fm = std::fmax(fm, std::abs(f[j])); }
        return rmax / fm;
    };
    // Stopping rules: (i) residual below 1e-12 (degenerate zero-mode directions need no Newton step, and S is stationary
    // along them); (ii) Newton step below tol; (iii) stagnation of the residual at the discretisation floor (< 1e-9).
    double prev = rel_residual(), prev2 = 1e300;
    for (int it = 1; it <= max_it; ++it) {
        if (prev < 1e-12) break;
        Mat J(D2);
        for (int i = 1; i < N; ++i) {
            const double x = std::log(std::fmax(f[i], 1e-300)) - s.tau[i];
            const double Gx = pot.dlam(x) + 0.25 * pot.d2lam(x);
            J[static_cast<std::size_t>(i) * n + i] -= 1.0 + 3.0 * f[i] * f[i] * G[i] + f[i] * f[i] * Gx;
        }
        for (int j = 0; j < n; ++j) {
            J[j] = D1[j] + (j == 0 ? 1.0 : 0.0);
            J[static_cast<std::size_t>(N) * n + j] = D1[static_cast<std::size_t>(N) * n + j] - (j == N ? 1.0 : 0.0);
        }
        std::vector<double> delta(R);
        for (double& v : delta) v = -v;
        if (!lu_solve(J, delta, n)) {
            s.message = "singular Newton matrix";
            return s;
        }
        const double n0 = norm(R);
        double step = 1.0, dcap = 0.0, fcap = 0.0;
        for (int j = 0; j < n; ++j) { dcap = std::fmax(dcap, std::abs(delta[j])); fcap = std::fmax(fcap, std::abs(f[j])); }
        if (dcap > 0.3 * fcap) step = 0.3 * fcap / dcap;   // trust-region-like cap on the Newton step
        std::vector<double> fn(n);
        for (;;) {
            for (int j = 0; j < n; ++j) fn[j] = f[j] + step * delta[j];
            residual(fn, Rn, Gn);
            if (norm(Rn) < n0 || step < 1e-4) break;
            step *= 0.5;
        }
        double dmax = 0.0, fmax = 0.0;
        for (int j = 0; j < n; ++j) {
            dmax = std::fmax(dmax, std::abs(step * delta[j]));
            fmax = std::fmax(fmax, std::abs(fn[j]));
        }
        f.swap(fn);
        R.swap(Rn);
        G.swap(Gn);
        s.iterations = it;
        if (dmax < tol * fmax) break;
        const double rr = rel_residual();
        if (rr < 1e-9 && rr > 0.5 * prev2) break;   // stagnation at the discretisation floor
        prev2 = prev;
        prev = rr;
        if (it == max_it) {
            s.message = "Newton iteration did not converge";
            return s;
        }
    }
    double fmax = 0.0, fmin = 1e300;
    for (double v : f) { fmax = std::fmax(fmax, v); fmin = std::fmin(fmin, v); }
    if (!(fmin > 0.0) || fmax < 0.1 * A) {
        s.message = "converged to a non-bounce solution (trivial or sign-changing)";
        return s;
    }
    s.fp.assign(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) s.fp[i] += D1[static_cast<std::size_t>(i) * n + j] * f[j];
    s.f = f;
    s.w = clenshaw_curtis(N);
    for (double& v : s.w) v *= L;
    s.ok = true;
    return s;
}

void functionals(const Potential& pot, const Solution& s, BounceResult& r) {
    const double tp2 = 2.0 * kPi2;
    double T = 0, V = 0, D = 0, W = 0;
    std::size_t imax = 0;
    for (std::size_t j = 0; j < s.f.size(); ++j) {
        const double f = s.f[j], x = std::log(f) - s.tau[j];
        const double lam = pot.lam(x), beta = pot.dlam(x), f4 = f * f * f * f;
        T += s.w[j] * 0.5 * (s.fp[j] - f) * (s.fp[j] - f);
        V += s.w[j] * 0.25 * lam * f4;
        D += s.w[j] * beta * f4;
        W += s.w[j] * std::abs(lam) * f4;
        if (f > s.f[imax]) imax = j;
    }
    T = tp2 * (T + s.f[0] * s.f[0]);   // + right tail beyond tau_c + L (f ~ e^{-tau}); node 0 is the right end
    V *= tp2;
    r.T = T;
    r.V = V;
    r.S = T + V;
    r.derrick = (T + 2.0 * V) / T;
    r.dilatation = D / W;
    r.lnR = s.tau[imax];
    r.lnphi0 = std::log(s.f.back()) - s.tau.back();   // node N is the left end, where f ~ phi(0) e^{tau}
}

}  // namespace

BounceResult solve_bounce(const Potential& pot, const AnalyticEstimates& guess, const BounceOptions& opt) {
    BounceResult r;
    if (!(guess.lambda_min < 0)) {
        r.message = "lambda_min >= 0: no bounce";
        return r;
    }
    const double A = std::sqrt(2.0 / std::abs(guess.lambda_min));
    const double lnh = guess.CGM_interior ? guess.lnh_CGM : guess.x_min + 5.0 / 6.0;
    const double tc = std::log(2.0 * A) - lnh;   // FL: f = A sech(tau - tc), phi(0) = 2 A e^{-tc}
    // Increase the resolution until both consistency identities hold (thin walls need more nodes). Keep the best
    // converged pair in case only the reduced-accuracy tier can be reached.
    int n1 = opt.nodes, n2 = opt.nodes_check;
    BounceResult best;
    double best_res = 1e300;
    for (;;) {
        const Solution s1 = collocate(pot, A, tc, n1, opt.half_width, opt.newton_tol, opt.max_newton);
        const Solution s2 = s1.ok ? collocate(pot, A, tc, n2, opt.half_width, opt.newton_tol, opt.max_newton) : s1;
        if (s1.ok && s2.ok) {
            BounceResult r1;
            functionals(pot, s1, r1);
            functionals(pot, s2, r);
            r.newton_iterations = s2.iterations;
            r.dS = std::abs(r.S - r1.S) + (std::abs(r.derrick) + std::abs(r.dilatation)) * std::abs(r.S);
            const double res = std::fmax(std::abs(r.derrick), std::abs(r.dilatation));
            if (res < opt.identity_tol) {
                r.ok = std::isfinite(r.S);
                r.message = "N=" + std::to_string(n2);
                return r;
            }
            if (res < best_res) {
                best = r;
                best_res = res;
                best.message = "reduced accuracy: identities " + std::to_string(res) + " at N=" + std::to_string(n2);
            }
            r.message = "identities not satisfied at N=" + std::to_string(n2);
        } else {
            r.message = "N=" + std::to_string(s1.ok ? n2 : n1) + ": " + (s1.ok ? s2.message : s1.message);
        }
        if (2 * n2 > opt.max_nodes) break;
        n1 *= 2;
        n2 *= 2;
    }
    if (best_res < opt.identity_tol_reduced && std::isfinite(best.S)) {
        best.ok = true;
        best.reduced_accuracy = true;
        return best;
    }
    r.ok = false;
    return r;
}

}  // namespace smvac
