// Minimal test helpers (no external framework). Each test executable returns non-zero on failure.
#pragma once

#include <cmath>
#include <cstdio>
#include <string>

namespace smvac_test {

inline int& failures() { static int n = 0; return n; }

inline void check(bool cond, const std::string& what) {
    std::printf("[%s] %s\n", cond ? " ok " : "FAIL", what.c_str());
    if (!cond) ++failures();
}

// |a - b| <= abs_tol + rel_tol |b|
inline void close(double a, double b, double rel_tol, double abs_tol, const std::string& what) {
    const bool ok = std::isfinite(a) && std::abs(a - b) <= abs_tol + rel_tol * std::abs(b);
    char buf[512];
    std::snprintf(buf, sizeof buf, "%s: got %.15g, expected %.15g (rel %.2e)", what.c_str(), a, b,
                  b != 0 ? std::abs(a - b) / std::abs(b) : std::abs(a - b));
    check(ok, buf);
}

inline int finish() {
    std::printf("%s (%d failure%s)\n", failures() ? "FAILED" : "PASSED", failures(), failures() == 1 ? "" : "s");
    return failures() ? 1 : 0;
}

}  // namespace smvac_test
