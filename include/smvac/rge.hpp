// SMVac — RG evolution from mu = M_t to M_Pl with C^1 dense output.
#pragma once

#include "smvac/couplings.hpp"

#include <vector>

namespace smvac {

// Fixed-step classical RK4 in t = ln(mu^2) (fourth-order convergence verified in
// docs/NUMERICAL_METHODS_REVIEW.md section 2), storing values AND exact derivatives at every node so that all
// couplings can be interpolated with cubic Hermite polynomials (C^1, O(dt^4)).
class RGSolution {
public:
    explicit RGSolution(const Inputs& in, double dt = 0.025);

    // Couplings at t (clamped to [t0, t_end]: couplings are frozen outside the integration range).
    Couplings at(double t) const;
    // dc/dt at t (zero outside the integration range).
    Couplings derivative_at(double t) const;

    double t0() const { return t_.front(); }
    double t_end() const { return t_.back(); }
    double t_planck() const { return t_planck_; }
    bool perturbative() const { return perturbative_; }   // reached M_Pl without |lambda|, y_t > 4 pi
    const Inputs& inputs() const { return in_; }
    std::size_t nodes() const { return t_.size(); }

private:
    Inputs in_;
    double dt_, t_planck_;
    bool perturbative_ = true;
    std::vector<double> t_;
    std::vector<Couplings> y_, dy_;
};

}  // namespace smvac
