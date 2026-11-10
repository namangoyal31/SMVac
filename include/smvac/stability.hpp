// SMVac — classification pipeline for problem P (docs/STAGE_C_PLAN.md section 2).
#pragma once

#include "smvac/analytic.hpp"
#include "smvac/bounce.hpp"
#include "smvac/couplings.hpp"
#include "smvac/decay.hpp"
#include "smvac/potential.hpp"

#include <string>

namespace smvac {

enum class Status {
    Stable,           // min lambda_P >= 0 up to M_Pl
    Metastable,       // numerical reference: p < 1
    Unstable,         // numerical reference: p >= 1
    NonPerturbative,  // |lambda| or y_t exceeded 4 pi below M_Pl
    NoBounce,         // lambda_min < 0 but no finite bounce (e.g. minimum on the frozen M_Pl boundary)
    SolverFailure,    // the numerical solver did not converge or failed its consistency checks
    InvalidInput
};
std::string to_string(Status s);

struct PointOptions {
    Prescription prescription = Prescription::CW1;
    double rge_dt = 0.025;
    bool numerical_bounce = true;
    BounceOptions bounce;
    DecayParameters decay;
    // Explicit beyond-the-SM deformation c6 phi^6/(6 Lambda^2): OFF by default, never part of SM results.
    bool bsm_phi6 = false;
    double bsm_c6 = 0.0, bsm_scale = 1.22e19;
    // A point is flagged Planck-sensitive if phi(0) of the bounce lies above M_Pl / 10^planck_margin_decades
    // (default 0: the bounce enters the region where the couplings are frozen) or lambda_min sits on that boundary.
    double planck_margin_decades = 0.0;
};

struct PointResult {
    Inputs in;
    Status status = Status::InvalidInput;
    std::string message;
    bool planck_sensitive = false;
    double lambda_min = 0, x_min = 0, log10_mu_instability = 0;   // first zero crossing of lambda_P
    AnalyticEstimates analytic;
    BounceResult bounce;
    // Decay probability (log10 p) from the numerical reference with Lambda_B = 1/R of the bounce,
    // and the literature analytical estimate (S_FL with Lambda_B = argmin lambda, Buttazzo et al. eq. 81).
    double log10_p_reference = 0, log10_p_FL_literature = 0, log10_p_FL_bounce_scale = 0;
    bool metastable_FL_literature = false;
};

PointResult analyze_point(const Inputs& in, const PointOptions& opt = PointOptions());

}  // namespace smvac
