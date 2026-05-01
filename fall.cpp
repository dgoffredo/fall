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

Step one_step(State current, Vec advanced) {
    const auto [t, dt, curr] = current;
    const auto [x, xp] = curr;
    const real xpp = acceleration(x);
    // Calculate the new (x, xp) when advanced by dt/2.
    // Then advance that result another dt/2.
    // Compare this result with `advanced` by calculating an "error".
    // If the error is acceptable, return the two-step result.
    // If the error is too large, recurse with dt/2 (using some values we just calculated).
    // If the error is too small, return the two-step but adjust dt -> 2*dt.
    const Vec half{
        .x = x + xp * dt/2,
        .xp = xp + xpp * dt/2
    };
    const Vec half_twice{
        .x = half.x + half.xp * dt/2,
        .xp = half.xp + acceleration(half.x) * dt/2
    };
    const real err_2 = error_squared(current.current, advanced, half_twice);
    constexpr real threshold = 1e-6;
    if (err_2 > threshold*threshold) {
        return one_step({t, dt/2, curr}, half);
    }
    if (err_2 < threshold*threshold/16) {
        return {{t+dt, dt*2, half_twice}, err_2};
    }
    return {{t+dt, dt, half_twice}, err_2};
}

void simulate(std::ostream& out, State state) {
    // initial point: `state` with zero error by definition
    out << Step{.next=state, .error_squared=0};

    auto& [t, dt, curr]  = state;
    auto& [x, xp] = curr;
    for (;;) {
        Vec advanced{
            .x = x + xp * dt,
            .xp = xp + acceleration(x) * dt
        };
        const Step step = one_step(state, advanced);
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

/* notes

dt, t, x, xp, xpp,
or in general:
dt, t, (xi, xip for i in x, y, z, ...)

1. calculate new x and xp from dt and force function
2. calculate two-step dt/2 values for x and xp
3. calculate the "error" between the two (x, xp) vectors
4. if the error is below threshold/4, double dt and goto 1
5. if the error is above threshold, halve dt and goto 1
6. otherwise commit new t, x, and xp

one_step(t, dt, x, xp) -> (new: t, dt, x, xp, error)

suppose we decide we need to take dt -> dt/2
    - we already calculated two dt/2 steps, so we can reuse the first as the "big step" in the recursion
    - so should the "result of half step" be a parameter?

suppose we decide we can afford to take dt -> 2*dt
    - we already calculated the dt step, so we can reuse that as the "little step" in the recursion

is it possible that we oscillate between {dt -> 2*dt} and {dt -> dt/2} infinitely?
    - not sure, but we can prevent it by falling back to the smaller dt instead of recalculating
    - for that matter, we could oscillate the other way as well
    - we need a state: "just increased dt," "just decreased dt," and "neither"
    - alternatively, could binary search for a dt that yields an error in the range [threshold/4, threshold]
        - that's probably best

there are optimizations here, but let's start with the simplest code I can manage, at the price of less efficiency

*/

