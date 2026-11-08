// End-to-end tests of the classification pipeline, including failure modes.
#include "smvac/stability.hpp"
#include "test_util.hpp"

#include <cmath>

using namespace smvac;
using smvac_test::check;
using smvac_test::close;

int main() {
    Inputs central;   // Buttazzo et al. central values: metastable
    const PointResult r = analyze_point(central);
    check(r.status == Status::Metastable, "central point is metastable (got " + to_string(r.status) + ")");
    check(!r.planck_sensitive, "central point bounce stays below M_Pl");
    check(r.analytic.S_FL <= r.bounce.S && r.bounce.S <= r.analytic.S_CGM, "S_FL <= S_num <= S_CGM at the central point");
    // S ~ 1856 and S_thr ~ 530 give log10 p ~ -576.
    check(r.log10_p_reference < -500, "decay probability is negligible at the central point");
    check(r.log10_mu_instability > 9.5 && r.log10_mu_instability < 12, "instability scale ~1e10-1e11 GeV");

    Inputs light = central;
    light.Mt = 170.0;   // below the stability boundary (~171.4 GeV for cw1)
    check(analyze_point(light).status == Status::Stable, "M_t = 170 GeV is absolutely stable");

    Inputs heavy = central;
    heavy.Mt = 240.0;
    check(analyze_point(heavy).status == Status::NonPerturbative, "M_t = 240 GeV is non-perturbative");

    Inputs bad = central;
    bad.Mt = -1.0;
    check(analyze_point(bad).status == Status::InvalidInput, "negative M_t is rejected");

    // The phi^6 deformation is opt-in and changes the action; the SM default never includes it.
    PointOptions bsm;
    bsm.bsm_phi6 = true;
    bsm.bsm_c6 = 1.0;
    const PointResult rb = analyze_point(central, bsm);
    check(rb.bounce.ok && rb.bounce.S > r.bounce.S, "phi^6 deformation (c6 = 1) increases the action when requested");
    check(PointOptions().bsm_phi6 == false, "SM default excludes the phi^6 deformation");

    // Analytical-only mode reports the FL classification and says so.
    PointOptions fast;
    fast.numerical_bounce = false;
    const PointResult rf = analyze_point(central, fast);
    check(rf.status == Status::Metastable && rf.message.find("analytical") != std::string::npos,
          "analytical-only mode is labelled as such");
    return smvac_test::finish();
}
