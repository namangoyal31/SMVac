// SMVac — Chebyshev series on [a, b]: interpolation at Chebyshev-Gauss-Lobatto points, Clenshaw evaluation,
// exact differentiation of the series.
#pragma once

#include <cmath>
#include <functional>
#include <vector>

namespace smvac {
namespace numerics {

class ChebyshevSeries {
public:
    ChebyshevSeries() = default;
    ChebyshevSeries(const std::function<double(double)>& f, double a, double b, int n) : a_(a), b_(b) {
        constexpr double pi = 3.14159265358979323846;
        std::vector<double> fv(n + 1);
        for (int k = 0; k <= n; ++k) fv[k] = f(map(std::cos(pi * k / n)));
        c_.assign(n + 1, 0.0);
        for (int j = 0; j <= n; ++j) {
            double s = 0.0;
            for (int k = 0; k <= n; ++k) {
                const double w = (k == 0 || k == n) ? 0.5 : 1.0;
                s += w * fv[k] * std::cos(pi * j * k / n);
            }
            c_[j] = 2.0 * s / n;
        }
        c_[0] *= 0.5;
        c_[n] *= 0.5;
    }

    double operator()(double x) const {
        const double u = (2.0 * x - a_ - b_) / (b_ - a_);
        double b1 = 0.0, b2 = 0.0;
        for (int j = static_cast<int>(c_.size()) - 1; j >= 1; --j) {
            const double t = 2.0 * u * b1 - b2 + c_[j];
            b2 = b1;
            b1 = t;
        }
        return u * b1 - b2 + c_[0];
    }

    ChebyshevSeries derivative() const {
        ChebyshevSeries d;
        d.a_ = a_;
        d.b_ = b_;
        const int n = static_cast<int>(c_.size()) - 1;
        if (n < 1) {
            d.c_ = {0.0};
            return d;
        }
        std::vector<double> cd(n + 1, 0.0);   // cd has degree n-1; extra slot stays zero
        for (int j = n - 1; j >= 0; --j) cd[j] = (j + 2 <= n ? cd[j + 2] : 0.0) + 2.0 * (j + 1) * c_[j + 1];
        cd[0] *= 0.5;
        const double scale = 2.0 / (b_ - a_);
        for (double& v : cd) v *= scale;
        cd.resize(n);
        d.c_ = cd;
        return d;
    }

    // max |last 6 coefficients| / max |coefficients|: a convergence diagnostic for the interpolation degree.
    double tail_ratio() const {
        double mx = 0.0, tail = 0.0;
        for (std::size_t j = 0; j < c_.size(); ++j) {
            mx = std::fmax(mx, std::abs(c_[j]));
            if (j + 6 >= c_.size()) tail = std::fmax(tail, std::abs(c_[j]));
        }
        return mx > 0 ? tail / mx : 0.0;
    }
    int degree() const { return static_cast<int>(c_.size()) - 1; }
    const std::vector<double>& coefficients() const { return c_; }
    double a() const { return a_; }
    double b() const { return b_; }

private:
    double map(double u) const { return 0.5 * (a_ + b_) + 0.5 * (b_ - a_) * u; }
    double a_ = -1.0, b_ = 1.0;
    std::vector<double> c_;
};

}  // namespace numerics
}  // namespace smvac
