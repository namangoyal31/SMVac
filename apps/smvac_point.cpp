// smvac_point — full analysis of one (M_h, M_t) point, printed as JSON.
// Usage: smvac_point <Mh> <Mt> [--prescription cw1|tree] [--matching buttazzo2013|extended|legacy] [--alpha3 x]
//                    [--nodes N] [--nodes-check N2] [--half-width L] [--h0 H0] [--no-bounce 1]
//                    [--bsm-phi6-c6 c6 [--bsm-phi6-scale Lambda]]
#include "smvac_common.hpp"

#include <cmath>
#include <cstdio>
#include <string>

namespace {
void num(const char* key, double v, bool comma = true) {
    if (std::isfinite(v)) std::printf("  \"%s\": %.15g%s\n", key, v, comma ? "," : "");
    else std::printf("  \"%s\": null%s\n", key, comma ? "," : "");
}
}  // namespace

int main(int argc, char** argv) {
    using namespace smvac;
    if (argc < 3) {
        std::fprintf(stderr, "usage: smvac_point <Mh> <Mt> [options]\n");
        return 2;
    }
    Inputs in;
    PointOptions opt;
    in.Mh = std::atof(argv[1]);
    in.Mt = std::atof(argv[2]);
    try {
        for (int i = 3; i + 1 < argc; i += 2)
            if (!smvac_app::parse_common(argv[i], argv[i + 1], in, opt))
                throw std::invalid_argument(std::string("unknown option ") + argv[i]);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
    const PointResult r = analyze_point(in, opt);
    std::printf("{\n");
    std::printf("  \"metadata\": {");
    const auto md = smvac_app::metadata(in, opt);
    for (std::size_t i = 0; i < md.size(); ++i)
        std::printf("%s\n    \"%s\": \"%s\"", i ? "," : "", md[i].first.c_str(), md[i].second.c_str());
    std::printf("\n  },\n");
    num("Mh", in.Mh);
    num("Mt", in.Mt);
    std::printf("  \"status\": \"%s\",\n  \"message\": \"%s\",\n  \"planck_sensitive\": %s,\n", to_string(r.status).c_str(),
                r.message.c_str(), r.planck_sensitive ? "true" : "false");
    num("lambda_min", r.lambda_min);
    num("log10_mu_lambda_min", r.x_min / std::log(10.0));
    num("log10_mu_instability", r.log10_mu_instability);
    num("lambda2_at_min", r.analytic.lambda2);
    num("S_FL_lower_bound", r.analytic.S_FL);
    num("S_curv", r.analytic.S_curv);
    num("S_CGM_upper_estimate", r.analytic.S_CGM);
    num("S_numerical", r.bounce.ok ? r.bounce.S : std::nan(""));
    num("dS_numerical", r.bounce.ok ? r.bounce.dS : std::nan(""));
    std::printf("  \"reduced_accuracy\": %s,\n", r.bounce.ok && r.bounce.reduced_accuracy ? "true" : "false");
    num("derrick_residual", r.bounce.ok ? r.bounce.derrick : std::nan(""));
    num("dilatation_residual", r.bounce.ok ? r.bounce.dilatation : std::nan(""));
    num("log10_phi0", r.bounce.ok ? r.bounce.lnphi0 / std::log(10.0) : std::nan(""));
    num("log10_LambdaB_bounce", r.bounce.ok ? -r.bounce.lnR / std::log(10.0) : std::nan(""));
    num("log10_p_numerical", r.bounce.ok ? r.log10_p_reference : std::nan(""));
    num("log10_p_FL_literature", r.log10_p_FL_literature);
    num("log10_p_FL_bounce_scale", r.log10_p_FL_bounce_scale, false);
    std::printf("}\n");
    return 0;
}
