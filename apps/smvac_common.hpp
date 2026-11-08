// SMVac — shared command-line handling for the applications.
#pragma once

#include "smvac/stability.hpp"

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace smvac_app {

constexpr const char* kVersion = "2.0.0";

// Parses --key value options common to the apps into Inputs/PointOptions. Unknown keys throw.
inline bool parse_common(const std::string& key, const char* val, smvac::Inputs& in, smvac::PointOptions& opt) {
    using namespace smvac;
    if (key == "--alpha3") in.alpha3 = std::atof(val);
    else if (key == "--mw") in.MW = std::atof(val);
    else if (key == "--matching") in.matching = parse_matching(val);
    else if (key == "--prescription") opt.prescription = parse_prescription(val);
    else if (key == "--rge-dt") opt.rge_dt = std::atof(val);
    else if (key == "--nodes") opt.bounce.nodes = std::atoi(val);
    else if (key == "--nodes-check") opt.bounce.nodes_check = std::atoi(val);
    else if (key == "--half-width") opt.bounce.half_width = std::atof(val);
    else if (key == "--h0") opt.decay.H0_km_s_Mpc = std::atof(val);
    else if (key == "--no-bounce") opt.numerical_bounce = (std::string(val) != "1");
    else if (key == "--bsm-phi6-c6") { opt.bsm_phi6 = true; opt.bsm_c6 = std::atof(val); }
    else if (key == "--bsm-phi6-scale") opt.bsm_scale = std::atof(val);
    else return false;
    return true;
}

inline std::vector<std::pair<std::string, std::string>> metadata(const smvac::Inputs& in, const smvac::PointOptions& opt) {
    using namespace smvac;
    auto d = [](double v) { char b[64]; std::snprintf(b, sizeof b, "%.10g", v); return std::string(b); };
    std::vector<std::pair<std::string, std::string>> m = {
        {"smvac_version", kVersion},
        {"problem", "tree-level O(4) bounce of V = lambda_P(phi) phi^4/4 (docs/STAGE_C_PLAN.md)"},
        {"matching", to_string(in.matching)}, {"alpha3", d(in.alpha3)}, {"MW", d(in.MW)},
        {"prescription", to_string(opt.prescription)}, {"rge_dt", d(opt.rge_dt)},
        {"bounce_method", "chebyshev-collocation"}, {"nodes", std::to_string(opt.bounce.nodes)},
        {"nodes_check", std::to_string(opt.bounce.nodes_check)}, {"half_width", d(opt.bounce.half_width)},
        {"numerical_bounce", opt.numerical_bounce ? "1" : "0"}, {"H0_km_s_Mpc", d(opt.decay.H0_km_s_Mpc)},
        {"light_cone", d(opt.decay.light_cone)}, {"bsm_phi6", opt.bsm_phi6 ? "1" : "0"}, {"bsm_c6", d(opt.bsm_c6)},
        {"bsm_scale", d(opt.bsm_scale)}};
    if (opt.bsm_phi6) m.push_back({"WARNING", "beyond-the-Standard-Model phi^6 deformation ACTIVE: not an SM result"});
    return m;
}

inline void print_metadata(std::FILE* out, const smvac::Inputs& in, const smvac::PointOptions& opt, const char* prefix) {
    for (const auto& kv : metadata(in, opt)) std::fprintf(out, "%s%s = %s\n", prefix, kv.first.c_str(), kv.second.c_str());
}

}  // namespace smvac_app
