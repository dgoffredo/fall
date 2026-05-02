#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numbers>
#include <ostream>

#include "density.h"

using real = long double;

real mass_below(const EarthRecord *lower_bound, real r) {
    if (lower_bound == std::end(earth_table)) {
        return (lower_bound - 1)->m;
    }

    const real R = (lower_bound - 1)->r;
    const real M = (lower_bound - 1)->m;
    const real ρ = (lower_bound - 1)->ρ;
    const real μ = lower_bound->μ;

    constexpr auto π = std::numbers::pi_v<real>;
    return M + π*μ*(r*r*r*r - R*R*R*R) + 4*π/3*(ρ - μ*R)*(r*r*r - R*R*R);
}

real acceleration(const EarthRecord *lower_bound, real x) {
    const real M = mass_below(lower_bound, std::abs(x));
    constexpr real G = 6.67430e-11;
    return std::copysign(M*G / (x*x), -x);
}

const EarthRecord *radius_lower_bound(real r) {
    // linear fit from gnuplot
    constexpr real a = 7.85178e-05;
    constexpr real b = -0.449567;
    const auto guess = std::min(std::size_t(std::ceil(a*r + b)), std::size(earth_table));

    const EarthRecord *record = std::begin(earth_table) + guess;

    // In case we overshot.
    while (record != std::begin(earth_table) && !((record - 1)->r < r)) {
        --record;
    }

    // In case we undershot.
    while (record != std::end(earth_table) && record->r < r) {
        ++record;
    }

    return record;
}

real acceleration(real x) {
    const EarthRecord *lower_bound = radius_lower_bound(std::abs(x));
    return acceleration(lower_bound, x);
}

struct Vec {
    real x;
    real xp;
};

// Arithmetic needed to form the Runge-Kutta stage combinations.
Vec operator+(Vec a, Vec b) { return {a.x + b.x, a.xp + b.xp}; }
Vec operator-(Vec a, Vec b) { return {a.x - b.x, a.xp - b.xp}; }
Vec operator*(real s, Vec v) { return {s * v.x, s * v.xp}; }

// ODE right-hand side: dy/dt = (xp, a(x))
Vec deriv(Vec y) { return {y.xp, acceleration(y.x)}; }

struct State {
    real t;
    real dt;
    Vec current;
};

real error_squared(Vec current, Vec coarse, Vec fine) {
    constexpr real ε_x = 1e-6;
    constexpr real ε_xp = 1e-6;
    const auto err_x = (coarse.x - fine.x) / (ε_x  + std::abs(current.x));
    const auto err_xp = (coarse.xp - fine.xp) / (ε_xp  + std::abs(current.xp));
    return err_x*err_x + err_xp*err_xp;
}

struct Step {
    State next;
    real error_squared;
};

std::ostream& operator<<(std::ostream& out, Step step) {
    const auto [state, error_squared] = step;
    const auto [t, dt, current] = state;
    const auto [x, xp] = current;
    return out << t << ' ' << x << ' ' << xp << ' ' << dt << ' ' << error_squared << '\n';
}

// Dormand-Prince (DOPRI5) adaptive step.
//
// Runge-Kutta methods generalize Euler's method by evaluating the derivative
// at several intermediate points ("stages") within [t, t+dt], then taking a
// weighted average as the step. More stages → higher accuracy order.
//
// DOPRI5 uses 6 stages and produces two solutions of different orders from
// the same evaluations (an "embedded pair"): a 5th-order result y5 (what we
// advance to) and a 4th-order result y4 (used only to estimate the error).
// Their difference y5 - y4 is a cheap local truncation error estimate.
Step dopri5_step(State state) {
    const auto [t, dt, y] = state;

    // Stages k1–k6: each is a derivative evaluation at a point c_i * dt
    // ahead of t, using a weighted sum of prior stages to get there.
    // Weights come from the Dormand-Prince Butcher tableau.
    const Vec k1 = deriv(y);
    const Vec k2 = deriv(y + dt * ((1.0L/5)*k1));
    const Vec k3 = deriv(y + dt * ((3.0L/40)*k1        + (9.0L/40)*k2));
    const Vec k4 = deriv(y + dt * ((44.0L/45)*k1       + (-56.0L/15)*k2      + (32.0L/9)*k3));
    const Vec k5 = deriv(y + dt * ((19372.0L/6561)*k1  + (-25360.0L/2187)*k2 + (64448.0L/6561)*k3 + (-212.0L/729)*k4));
    const Vec k6 = deriv(y + dt * ((9017.0L/3168)*k1   + (-355.0L/33)*k2     + (46732.0L/5247)*k3 + (49.0L/176)*k4   + (-5103.0L/18656)*k5));

    // 5th-order solution. (k2's weight is 0 in this combination.)
    const Vec y5 = y + dt * ((35.0L/384)*k1 + (500.0L/1113)*k3 + (125.0L/192)*k4 + (-2187.0L/6784)*k5 + (11.0L/84)*k6);

    // 7th stage at the 5th-order result — needed for the error estimate.
    // DOPRI5 has the FSAL ("First Same As Last") property: k7 here equals
    // k1 for the next step, saving one evaluation per accepted step.
    // We don't exploit that optimization here.
    const Vec k7 = deriv(y5);

    // err_vec ≈ y5 - y4: difference between 5th- and 4th-order solutions.
    // Coefficients e_i = b5_i - b4_i from the Butcher tableau.
    const Vec err_vec = dt * ((71.0L/57600)*k1 + (-71.0L/16695)*k3 + (71.0L/1920)*k4 + (-17253.0L/339200)*k5 + (22.0L/525)*k6 + (-1.0L/40)*k7);

    // Normalize error relative to solution magnitude.
    // y5 - err_vec = y4 (the coarser 4th-order solution).
    const real err2 = error_squared(y, y5 - err_vec, y5);

    // Step-size control: scale dt by (tolerance/err)^(1/5).
    // The exponent 1/5 = 1/(p+1) for the 4th-order component: the local
    // error scales as dt^5, so reducing error by (tolerance/err) requires
    // scaling dt by (tolerance/err)^(1/5).
    //
    // dt_max prevents the step from growing so large that 10,000 output
    // lines cover an unreasonable span of time.  Near the surface the
    // gravitational acceleration is nearly constant, so DOPRI5's 5th-order
    // error is almost zero and any step would be accepted — without this
    // cap the step balloons to ~100 s in a handful of iterations, causing
    // `head -10000` to span ~165 oscillation periods instead of ~2.
    constexpr real tolerance = 1e-6;
    constexpr real safety = 0.9L;
    constexpr real dt_max = 1.0L;
    const real err_norm = std::sqrt(err2);

    if (err_norm > tolerance) {
        // Too much error: shrink dt and retry.
        const real scale = safety * std::pow(tolerance / err_norm, 0.2L);
        return dopri5_step({t, dt * std::max(scale, 0.1L), y});
    }

    // Accepted. Grow dt if the error is comfortably below tolerance.
    const real scale = err_norm > 0 ? safety * std::pow(tolerance / err_norm, 0.2L) : 5.0L;
    return {{t + dt, std::min(dt * std::min(scale, 5.0L), dt_max), y5}, err2};
}

void simulate(std::ostream& out, State state) {
    out << Step{.next=state, .error_squared=0};
    for (;;) {
        const Step step = dopri5_step(state);
        state = step.next;
        out << step;
    }
}

int main() {
    const State start{
        .t = 0,
        .dt = 0.001,
        .current = Vec{.x = 6371000, .xp = 0 }
    };
    simulate(std::cout, start);
}


