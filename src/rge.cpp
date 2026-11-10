// SMVac — RG evolution.
#include "smvac/rge.hpp"
#include "smvac/beta_functions.hpp"
#include "smvac/constants.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace smvac {
namespace {

using Arr = std::array<double, Couplings::kSize>;

Arr deriv(const Arr& y) { return coupling_derivatives(Couplings::from_array(y)).to_array(); }

Arr axpy(const Arr& y, double h, const Arr& k) {
    Arr r;
    for (int i = 0; i < Couplings::kSize; ++i) r[i] = y[i] + h * k[i];
    return r;
}

bool sane(const Arr& y) {
    for (double v : y)
        if (!std::isfinite(v)) return false;
    return std::abs(y[6]) <= 4.0 * kPi && y[3] <= 4.0 * kPi && y[3] > 0.0;
}

}  // namespace

RGSolution::RGSolution(const Inputs& in, double dt) : in_(in), dt_(dt), t_planck_(2.0 * std::log(kPlanckMass)) {
    Arr y = match_at_top(in).to_array();
    double t = 2.0 * std::log(in.Mt);
    t_.push_back(t);
    y_.push_back(Couplings::from_array(y));
    dy_.push_back(Couplings::from_array(deriv(y)));
    const int nsteps = static_cast<int>(std::ceil((t_planck_ - t) / dt));
    const double h = (t_planck_ - t) / nsteps;   // uniform steps landing exactly on M_Pl
    for (int n = 0; n < nsteps; ++n) {
        const Arr k1 = deriv(y);
        const Arr k2 = deriv(axpy(y, 0.5 * h, k1));
        const Arr k3 = deriv(axpy(y, 0.5 * h, k2));
        const Arr k4 = deriv(axpy(y, h, k3));
        Arr yn;
        for (int i = 0; i < Couplings::kSize; ++i) yn[i] = y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        if (!sane(yn)) {
            perturbative_ = false;
            break;
        }
        y = yn;
        t = (n + 1 == nsteps) ? t_planck_ : t + h;
        t_.push_back(t);
        y_.push_back(Couplings::from_array(y));
        dy_.push_back(Couplings::from_array(deriv(y)));
    }
}

Couplings RGSolution::at(double t) const {
    if (t <= t_.front()) return y_.front();
    if (t >= t_.back()) return y_.back();
    const std::size_t i = std::min<std::size_t>(std::upper_bound(t_.begin(), t_.end(), t) - t_.begin() - 1, t_.size() - 2);
    const double h = t_[i + 1] - t_[i], s = (t - t_[i]) / h;
    const double h00 = (1 + 2 * s) * (1 - s) * (1 - s), h10 = s * (1 - s) * (1 - s);
    const double h01 = s * s * (3 - 2 * s), h11 = s * s * (s - 1);
    const Arr a = y_[i].to_array(), b = y_[i + 1].to_array(), da = dy_[i].to_array(), db = dy_[i + 1].to_array();
    Arr r;
    for (int k = 0; k < Couplings::kSize; ++k) r[k] = h00 * a[k] + h10 * h * da[k] + h01 * b[k] + h11 * h * db[k];
    return Couplings::from_array(r);
}

Couplings RGSolution::derivative_at(double t) const {
    if (t <= t_.front() || t >= t_.back()) return Couplings{};
    const std::size_t i = std::min<std::size_t>(std::upper_bound(t_.begin(), t_.end(), t) - t_.begin() - 1, t_.size() - 2);
    const double h = t_[i + 1] - t_[i], s = (t - t_[i]) / h;
    const double d00 = 6 * s * (s - 1) / h, d10 = (1 - s) * (1 - 3 * s), d01 = -6 * s * (s - 1) / h, d11 = s * (3 * s - 2);
    const Arr a = y_[i].to_array(), b = y_[i + 1].to_array(), da = dy_[i].to_array(), db = dy_[i + 1].to_array();
    Arr r;
    for (int k = 0; k < Couplings::kSize; ++k) r[k] = d00 * a[k] + d10 * da[k] + d01 * b[k] + d11 * db[k];
    return Couplings::from_array(r);
}

}  // namespace smvac
