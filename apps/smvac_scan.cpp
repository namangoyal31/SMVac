// smvac_scan — evaluate a rectangular (M_h, M_t) grid (or a contiguous chunk of it) and write CSV with a metadata
// header. Parallelism is by process: python/smvac/run_scan.py launches chunks and merges them.
// Usage: smvac_scan --mh MIN MAX STEP --mt MIN MAX STEP --out FILE [--chunk START COUNT] [common options]
#include "smvac_common.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    using namespace smvac;
    Inputs in;
    PointOptions opt;
    double mh[3] = {120, 130, 1}, mt[3] = {166, 182, 1};
    long start = 0, count = -1;
    std::string out;
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string k = argv[i];
            if (k == "--mh" && i + 3 < argc) { for (int j = 0; j < 3; ++j) mh[j] = std::atof(argv[++i]); }
            else if (k == "--mt" && i + 3 < argc) { for (int j = 0; j < 3; ++j) mt[j] = std::atof(argv[++i]); }
            else if (k == "--chunk" && i + 2 < argc) { start = std::atol(argv[++i]); count = std::atol(argv[++i]); }
            else if (k == "--out" && i + 1 < argc) out = argv[++i];
            else if (i + 1 < argc && smvac_app::parse_common(k, argv[i + 1], in, opt)) ++i;
            else throw std::invalid_argument("unknown or incomplete option " + k);
        }
        if (out.empty()) throw std::invalid_argument("--out is required");
        if (!(mh[2] > 0 && mt[2] > 0)) throw std::invalid_argument("grid steps must be positive");
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
    // Grid nodes computed from integer indices (no accumulation of round-off).
    const long nh = static_cast<long>(std::floor((mh[1] - mh[0]) / mh[2] + 1e-9)) + 1;
    const long nt = static_cast<long>(std::floor((mt[1] - mt[0]) / mt[2] + 1e-9)) + 1;
    const long total = nh * nt;
    if (count < 0 || start + count > total) count = total - start;
    std::FILE* f = std::fopen(out.c_str(), "w");
    if (!f) {
        std::fprintf(stderr, "error: cannot open %s\n", out.c_str());
        return 1;
    }
    smvac_app::print_metadata(f, in, opt, "# ");
    std::fprintf(f, "# grid_mh = %.10g %.10g %.10g\n# grid_mt = %.10g %.10g %.10g\n# chunk = %ld %ld of %ld\n", mh[0],
                 mh[1], mh[2], mt[0], mt[1], mt[2], start, count, total);
    std::fprintf(f, "Mh,Mt,status,planck_sensitive,lambda_min,log10_mu_min,log10_mu_inst,lambda2,S_FL,S_curv,S_CGM,"
                    "S_num,dS_num,derrick,dilatation,log10_phi0,log10_LambdaB,log10_p_num,log10_p_FL_lit,"
                    "log10_p_FL_bscale,reduced_accuracy,message\n");
    for (long idx = start; idx < start + count; ++idx) {
        Inputs p = in;
        p.Mt = mt[0] + (idx / nh) * mt[2];
        p.Mh = mh[0] + (idx % nh) * mh[2];
        const PointResult r = analyze_point(p, opt);
        const bool b = r.bounce.ok;
        const double nan = std::nan("");
        std::string msg = r.message;
        for (char& c : msg)
            if (c == ',' || c == '\n') c = ';';
        std::fprintf(f, "%.10g,%.10g,%s,%d,%.12g,%.10g,%.10g,%.10g,%.15g,%.15g,%.15g,%.15g,%.6g,%.3g,%.3g,%.10g,%.10g,%.10g,%.10g,%.10g,%d,%s\n",
                     p.Mh, p.Mt, to_string(r.status).c_str(), r.planck_sensitive ? 1 : 0, r.lambda_min,
                     r.x_min / std::log(10.0), r.log10_mu_instability, r.analytic.lambda2, r.analytic.S_FL,
                     r.analytic.S_curv, r.analytic.S_CGM, b ? r.bounce.S : nan, b ? r.bounce.dS : nan,
                     b ? r.bounce.derrick : nan, b ? r.bounce.dilatation : nan, b ? r.bounce.lnphi0 / std::log(10.0) : nan,
                     b ? -r.bounce.lnR / std::log(10.0) : nan, b ? r.log10_p_reference : nan, r.log10_p_FL_literature,
                     r.log10_p_FL_bounce_scale, b && r.bounce.reduced_accuracy ? 1 : 0, msg.c_str());
        std::fflush(f);
    }
    std::fclose(f);
    return 0;
}
