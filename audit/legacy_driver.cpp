// Audit driver: evaluates the legacy classify_stability() (apps/solver_numerical.cpp, compiled with -Dmain=legacy_main) at (Mh, Mt) pairs read from stdin.
#include <tuple>
#include <cstdio>
std::tuple<int,double,double,double> classify_stability(double Mh, double Mt);
int main() {
    double mh, mt;
    std::printf("Mh,Mt,Stability,S_exact,S_approx\n");
    while (std::scanf("%lf %lf", &mh, &mt) == 2) {
        auto r = classify_stability(mh, mt);
        std::printf("%.17g,%.17g,%d,%.17g,%.17g\n", mh, mt, std::get<0>(r), std::get<1>(r), std::get<2>(r));
    }
    return 0;
}
