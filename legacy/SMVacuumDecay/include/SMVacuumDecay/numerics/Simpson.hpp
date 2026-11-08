#pragma once

namespace SMVacuumDecay {
namespace numerics {

template <typename Func>
double simpson_integrate(double a, double b, int N, Func f) {
    if (N % 2 != 0) N++; // N must be even
    double h = (b - a) / N;
    double sum = (f(a) + f(b)) / 3.0;
    
    for (int i = 1; i < N; ++i) {
        double x = a + i * h;
        double weight = (i % 2 == 1) ? 4.0 / 3.0 : 2.0 / 3.0;
        sum += weight * f(x);
    }
    
    return sum * h;
}

}
}
