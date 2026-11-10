#pragma once
#include <cmath>

namespace SMVacuumDecay {
namespace numerics {

template <typename State, typename DerivFunc>
State rk4_single_step(const State& y, double t, double dt, DerivFunc dydt) {
    State k1 = dydt(y, t);
    State k2 = dydt(y + k1 * (0.5 * dt), t + 0.5 * dt);
    State k3 = dydt(y + k2 * (0.5 * dt), t + 0.5 * dt);
    State k4 = dydt(y + k3 * dt, t + dt);
    
    return y + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
}

template <typename State, typename DerivFunc, typename ErrorFunc, typename ValidFunc>
bool rk4_adaptive_step(State& y, double& t, double& dt, DerivFunc dydt, ErrorFunc get_error, ValidFunc is_valid, double TOL = 1e-13) {
    while(true) {
        State y1 = rk4_single_step(y, t, dt, dydt);
        State y_half = rk4_single_step(y, t, dt/2.0, dydt);
        State y2 = rk4_single_step(y_half, t + dt/2.0, dt/2.0, dydt);
        
        double error = get_error(y1, y2);
        if (!std::isfinite(error) || !is_valid(y2)) {
            return false;
        }
        if (error < TOL) {
            y = y2;
            t += dt;
            if (error < TOL/10.0) dt *= 1.5;
            return true;
        } else if (dt < 1e-5) {
            return false;
        } else {
            dt *= 0.5;
        }
    }
}

} // namespace numerics
} // namespace SMVacuumDecay
