// SMVac — boundary conditions at mu = M_t.
#include "smvac/constants.hpp"
#include "smvac/couplings.hpp"

#include <cmath>
#include <stdexcept>

namespace smvac {

Couplings match_at_top(const Inputs& in) {
    if (!(in.Mt > 0.0) || !(in.Mh >= 0.0) || !std::isfinite(in.Mt) || !std::isfinite(in.Mh))
        throw std::invalid_argument("match_at_top: require M_t > 0 and M_h >= 0");
    const double das = (in.alpha3 - 0.1184) / 0.0007;
    const double dmw = (in.MW - 80.384) / 0.014;
    const double dt = in.Mt - 173.34, dh = in.Mh - 125.15;
    Couplings c;
    c.yb = std::sqrt(2.0) * kMbTree / kVev;
    c.ytau = std::sqrt(2.0) * kMtauTree / kVev;
    switch (in.matching) {
        case MatchingScheme::Buttazzo2013:
        case MatchingScheme::Extended: {
            // Buttazzo et al. (2013), eqs. (55), (57), (58), (59), (60).
            c.lambda = 0.12604 + 0.00206 * dh - 0.00004 * dt;
            if (in.matching == MatchingScheme::Extended)  // tree-level M_h^2/(2 v^2) shape, same central value
                c.lambda = in.Mh * in.Mh / (2.0 * kVev * kVev) + (0.12604 - 125.15 * 125.15 / (2.0 * kVev * kVev))
                           - 0.00004 * dt;
            c.yt = 0.93690 + 0.00556 * dt - 0.00042 * das;
            c.g3 = 1.1666 + 0.00314 * das - 0.00046 * dt;
            c.g2 = 0.64779 + 0.00004 * dt + 0.00011 * dmw;
            c.g1 = std::sqrt(5.0 / 3.0) * (0.35830 + 0.00011 * dt - 0.00020 * dmw);
            break;
        }
        case MatchingScheme::Legacy: {
            c.lambda = in.Mh * in.Mh / (2.0 * kVev * kVev) - 0.00313 - 0.00004 * dt;
            c.yt = std::sqrt(2.0) * in.Mt / kVev - 0.0587 - 0.00042 * das;
            c.g1 = 0.46266 + 0.00006 * dt;
            c.g2 = 0.65355 + 0.00002 * dt;
            c.g3 = 1.1666 + 0.00314 * das - 0.00046 * dt;
            break;
        }
    }
    if (!(c.yt > 0.0)) throw std::invalid_argument("match_at_top: y_t(M_t) <= 0 for these inputs");
    return c;
}

MatchingScheme parse_matching(const std::string& name) {
    if (name == "buttazzo2013") return MatchingScheme::Buttazzo2013;
    if (name == "extended") return MatchingScheme::Extended;
    if (name == "legacy") return MatchingScheme::Legacy;
    throw std::invalid_argument("unknown matching scheme: " + name);
}

std::string to_string(MatchingScheme s) {
    switch (s) {
        case MatchingScheme::Buttazzo2013: return "buttazzo2013";
        case MatchingScheme::Extended: return "extended";
        case MatchingScheme::Legacy: return "legacy";
    }
    return "?";
}

}  // namespace smvac
