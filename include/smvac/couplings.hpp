// SMVac — MS-bar coupling vector and inputs.
// Conventions (Buttazzo et al. 2013, arXiv:1307.3536): V = -m^2|H|^2/2 + lambda|H|^4, so M_h^2 = 2 lambda v^2;
// g1 is GUT-normalised (g1^2 = 5/3 gY^2); RG time t = ln(mu^2/GeV^2).
#pragma once

#include <array>
#include <string>

namespace smvac {

struct Couplings {
    double g1 = 0, g2 = 0, g3 = 0, yt = 0, yb = 0, ytau = 0, lambda = 0;

    static constexpr int kSize = 7;
    std::array<double, kSize> to_array() const { return {g1, g2, g3, yt, yb, ytau, lambda}; }
    static Couplings from_array(const std::array<double, kSize>& a) {
        return Couplings{a[0], a[1], a[2], a[3], a[4], a[5], a[6]};
    }
};

enum class MatchingScheme {
    Buttazzo2013,   // linear NNLO fits at mu = M_t, eqs. (55), (57)-(60): valid near (125.15, 173.34) GeV
    Extended,       // as Buttazzo2013 but lambda with tree-level M_h^2 dependence; for EXPLORATORY scans only
    Legacy          // original code (wrong g2, g1 slope, y_t slope); for traceability only
};

struct Inputs {
    double Mh = 125.15;      // Higgs pole mass [GeV]
    double Mt = 173.34;      // top pole mass [GeV]
    double alpha3 = 0.1184;  // alpha_s(M_Z)
    double MW = 80.384;      // W pole mass [GeV]
    MatchingScheme matching = MatchingScheme::Buttazzo2013;
};

// MS-bar couplings at mu = M_t. Throws std::invalid_argument for unphysical inputs (M_t <= 0, y_t <= 0, ...).
Couplings match_at_top(const Inputs& in);

MatchingScheme parse_matching(const std::string& name);
std::string to_string(MatchingScheme s);

}  // namespace smvac
