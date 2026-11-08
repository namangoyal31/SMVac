// SMVac — RG-improved quartic lambda(x), x = ln(phi/GeV), for V(phi) = lambda(phi) phi^4 / 4.
#pragma once

#include "smvac/couplings.hpp"
#include "smvac/numerics/chebyshev.hpp"
#include "smvac/rge.hpp"

#include <memory>
#include <string>
#include <vector>

namespace smvac {

class Potential {
public:
    virtual ~Potential() = default;
    virtual double lam(double x) const = 0;     // lambda
    virtual double dlam(double x) const = 0;    // d lambda / dx
    virtual double d2lam(double x) const = 0;   // d^2 lambda / dx^2
    // Points where the representation is not smooth (frozen-coupling boundaries), for quadrature splitting.
    virtual std::vector<double> kinks() const { return {}; }
    // V'(phi) = phi^3 G(x), G = lambda + beta/4.
    double G(double x) const { return lam(x) + 0.25 * dlam(x); }
};

enum class Prescription {
    Tree,   // MS-bar lambda at mu = phi
    CW1     // one-loop Coleman-Weinberg lambda_eff at mu = phi (Buttazzo et al. eq. 106 with Gamma = 0, no Goldstones)
};
Prescription parse_prescription(const std::string& name);
std::string to_string(Prescription p);

// One-loop large-field effective quartic (Landau gauge, MS-bar), with the Z coupling g2^2 + gY^2 = g2^2 + 3/5 g1^2.
double lambda_cw1(const Couplings& c);

// lambda(x) as a Chebyshev series on [x_lo, x_hi]; constant outside (frozen couplings).
class FrozenChebyshevPotential : public Potential {
public:
    FrozenChebyshevPotential(const std::function<double(double)>& f, double x_lo, double x_hi, int max_degree = 512);
    double lam(double x) const override;
    double dlam(double x) const override;
    double d2lam(double x) const override;
    std::vector<double> kinks() const override { return {x_lo_, x_hi_}; }
    double x_lo() const { return x_lo_; }
    double x_hi() const { return x_hi_; }
    int degree() const { return s_.degree(); }
    double tail_ratio() const { return s_.tail_ratio(); }

private:
    double x_lo_, x_hi_, lam_lo_, lam_hi_;
    numerics::ChebyshevSeries s_, d1_, d2_;
};

// The RG-improved SM potential of problem P (docs/STAGE_C_PLAN.md section 2).
std::unique_ptr<FrozenChebyshevPotential> make_sm_potential(const RGSolution& rg, Prescription p);

// Explicit beyond-the-SM deformation: V -> V + c6 phi^6 / (6 Lambda^2), i.e. lambda -> lambda + (2/3) c6 phi^2/Lambda^2.
// Never used unless requested; outputs record it in their metadata.
class Phi6Deformation : public Potential {
public:
    Phi6Deformation(const Potential& base, double c6, double scale) : base_(base), c6_(c6), scale_(scale) {}
    double lam(double x) const override;
    double dlam(double x) const override;
    double d2lam(double x) const override;
    std::vector<double> kinks() const override { return base_.kinks(); }

private:
    const Potential& base_;
    double c6_, scale_;
};

// Analytic test potentials.
class ConstantPotential : public Potential {
public:
    explicit ConstantPotential(double l0) : l0_(l0) {}
    double lam(double) const override { return l0_; }
    double dlam(double) const override { return 0.0; }
    double d2lam(double) const override { return 0.0; }

private:
    double l0_;
};

class ParabolaPotential : public Potential {   // lambda = l0 + (l2/2)(x - x0)^2
public:
    ParabolaPotential(double l0, double l2, double x0) : l0_(l0), l2_(l2), x0_(x0) {}
    double lam(double x) const override { return l0_ + 0.5 * l2_ * (x - x0_) * (x - x0_); }
    double dlam(double x) const override { return l2_ * (x - x0_); }
    double d2lam(double) const override { return l2_; }

private:
    double l0_, l2_, x0_;
};

}  // namespace smvac
