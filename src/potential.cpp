// SMVac — RG-improved potentials.
#include "smvac/potential.hpp"
#include "smvac/constants.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace smvac {

double lambda_cw1(const Couplings& c) {
    const double xt = c.yt * c.yt, xw = c.g2 * c.g2, xz = 0.6 * c.g1 * c.g1 + xw;
    return c.lambda + (-3.0 * xt * xt * (std::log(0.5 * xt) - 1.5) + 0.375 * xw * xw * (std::log(0.25 * xw) - 5.0 / 6.0)
                       + 0.1875 * xz * xz * (std::log(0.25 * xz) - 5.0 / 6.0)) / kLoop;
}

Prescription parse_prescription(const std::string& name) {
    if (name == "tree") return Prescription::Tree;
    if (name == "cw1") return Prescription::CW1;
    throw std::invalid_argument("unknown prescription: " + name);
}

std::string to_string(Prescription p) { return p == Prescription::Tree ? "tree" : "cw1"; }

FrozenChebyshevPotential::FrozenChebyshevPotential(const std::function<double(double)>& f, double x_lo, double x_hi,
                                                   int max_degree)
    : x_lo_(x_lo), x_hi_(x_hi) {
    if (!(x_hi > x_lo)) throw std::invalid_argument("FrozenChebyshevPotential: empty interval");
    for (int n = 32; n <= max_degree; n *= 2) {   // raise the degree until the tail is at round-off level
        s_ = numerics::ChebyshevSeries(f, x_lo, x_hi, n);
        if (s_.tail_ratio() < 1e-14) break;
    }
    d1_ = s_.derivative();
    d2_ = d1_.derivative();
    lam_lo_ = s_(x_lo);
    lam_hi_ = s_(x_hi);
}

double FrozenChebyshevPotential::lam(double x) const {
    if (x <= x_lo_) return lam_lo_;
    if (x >= x_hi_) return lam_hi_;
    return s_(x);
}
double FrozenChebyshevPotential::dlam(double x) const { return (x <= x_lo_ || x >= x_hi_) ? 0.0 : d1_(x); }
double FrozenChebyshevPotential::d2lam(double x) const { return (x <= x_lo_ || x >= x_hi_) ? 0.0 : d2_(x); }

std::unique_ptr<FrozenChebyshevPotential> make_sm_potential(const RGSolution& rg, Prescription p) {
    const double x_lo = 0.5 * std::max(rg.t0(), 2.0 * std::log(kVev));
    const double x_hi = 0.5 * rg.t_end();
    std::function<double(double)> f;
    if (p == Prescription::Tree)
        f = [&rg](double x) { return rg.at(2.0 * x).lambda; };
    else
        f = [&rg](double x) { return lambda_cw1(rg.at(2.0 * x)); };
    return std::unique_ptr<FrozenChebyshevPotential>(new FrozenChebyshevPotential(f, x_lo, x_hi));
}

double Phi6Deformation::lam(double x) const {
    return base_.lam(x) + (2.0 / 3.0) * c6_ * std::exp(2.0 * (x - std::log(scale_)));
}
double Phi6Deformation::dlam(double x) const {
    return base_.dlam(x) + (4.0 / 3.0) * c6_ * std::exp(2.0 * (x - std::log(scale_)));
}
double Phi6Deformation::d2lam(double x) const {
    return base_.d2lam(x) + (8.0 / 3.0) * c6_ * std::exp(2.0 * (x - std::log(scale_)));
}

}  // namespace smvac
