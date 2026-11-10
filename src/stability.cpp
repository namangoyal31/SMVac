// SMVac — classification pipeline.
#include "smvac/stability.hpp"
#include "smvac/constants.hpp"
#include "smvac/numerics/brent.hpp"
#include "smvac/rge.hpp"

#include <cmath>
#include <memory>
#include <stdexcept>

namespace smvac {

std::string to_string(Status s) {
    switch (s) {
        case Status::Stable: return "stable";
        case Status::Metastable: return "metastable";
        case Status::Unstable: return "unstable";
        case Status::NonPerturbative: return "non-perturbative";
        case Status::NoBounce: return "no-bounce";
        case Status::SolverFailure: return "solver-failure";
        case Status::InvalidInput: return "invalid-input";
    }
    return "?";
}

namespace {

// First zero crossing of lambda(x) above x_lo (log10 of the scale in GeV), or NaN if lambda stays positive.
double instability_scale(const Potential& pot, double x_lo, double x_hi) {
    const int n = 2000;
    double xa = x_lo, la = pot.lam(x_lo);
    if (la <= 0) return x_lo / std::log(10.0);
    for (int i = 1; i <= n; ++i) {
        const double xb = x_lo + (x_hi - x_lo) * i / n, lb = pot.lam(xb);
        if (lb < 0) {
            const double r = numerics::brent_root([&](double x) { return pot.lam(x); }, xa, xb, 1e-12);
            return r / std::log(10.0);
        }
        xa = xb;
        la = lb;
    }
    return std::nan("");
}

}  // namespace

PointResult analyze_point(const Inputs& in, const PointOptions& opt) {
    PointResult r;
    r.in = in;
    std::unique_ptr<RGSolution> rg;
    try {
        rg.reset(new RGSolution(in, opt.rge_dt));
    } catch (const std::invalid_argument& e) {
        r.status = Status::InvalidInput;
        r.message = e.what();
        return r;
    }
    if (!rg->perturbative()) {
        r.status = Status::NonPerturbative;
        r.message = "|lambda| or y_t exceeded 4 pi below M_Pl";
        return r;
    }
    const auto sm = make_sm_potential(*rg, opt.prescription);
    std::unique_ptr<Phi6Deformation> bsm;
    const Potential* pot = sm.get();
    if (opt.bsm_phi6) {
        bsm.reset(new Phi6Deformation(*sm, opt.bsm_c6, opt.bsm_scale));
        pot = bsm.get();
    }
    // The phi^6 deformation is not frozen at M_Pl, so extend the search range when it is active.
    const double x_lo = sm->x_lo(), x_hi = opt.bsm_phi6 ? sm->x_hi() + 6.0 : sm->x_hi();
    r.analytic = analytic_estimates(*pot, x_lo, x_hi);
    r.lambda_min = r.analytic.lambda_min;
    r.x_min = r.analytic.x_min;
    if (r.lambda_min >= 0) {
        r.status = Status::Stable;
        r.log10_mu_instability = std::nan("");
        return r;
    }
    r.log10_mu_instability = instability_scale(*pot, x_lo, x_hi);
    r.log10_p_FL_literature = log10_decay_probability(r.analytic.S_FL, r.analytic.lnLambdaB_literature, opt.decay);
    r.log10_p_FL_bounce_scale = log10_decay_probability(r.analytic.S_FL, r.analytic.lnLambdaB_FL, opt.decay);
    r.metastable_FL_literature = r.log10_p_FL_literature < 0;
    const double x_planck = std::log(kPlanckMass);
    if (r.analytic.min_at_boundary && std::abs(r.x_min - sm->x_hi()) < 1e-6) {
        r.planck_sensitive = true;
        r.status = Status::NoBounce;
        r.message = "lambda decreases up to M_Pl (frozen boundary): the infimum 8pi^2/(3|lambda(M_Pl)|) is not attained";
        return r;
    }
    if (!opt.numerical_bounce) {
        r.status = r.metastable_FL_literature ? Status::Metastable : Status::Unstable;
        r.message = "classification from the analytical FL estimate only (numerical bounce disabled)";
        return r;
    }
    r.bounce = solve_bounce(*pot, r.analytic, opt.bounce);
    if (!r.bounce.ok) {
        r.status = Status::SolverFailure;
        r.message = r.bounce.message;
        return r;
    }
    r.planck_sensitive = (r.bounce.lnphi0 > x_planck - opt.planck_margin_decades * std::log(10.0));
    r.log10_p_reference = log10_decay_probability(r.bounce.S, -r.bounce.lnR, opt.decay);
    r.status = r.log10_p_reference < 0 ? Status::Metastable : Status::Unstable;
    return r;
}

}  // namespace smvac
