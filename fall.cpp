#include <cmath>
#include <iostream>
#include <ostream>

using real = long double;

real accel(real r) {
    // copsign(mag, sgn) -> number with magnitude of mag and sign of sgn
    return std::copysign(1 / (r*r), -r);
}

struct Vec {
    real x;
    real xp;
    // TODO: Could store accel(x) here, but might be faster to recalculate all the time.
};

struct State {
    real t;
    real dt;
    Vec current;
};

real error(Vec current, Vec coarse, Vec fine) {
    constexpr real ε_x = 1e-6;
    constexpr real ε_xp = 1e-6;
    return std::hypot(
        (coarse.x - fine.x) / (ε_x  + current.x),
        (coarse.xp - fine.xp) / (ε_xp  + current.xp));
}

struct Step {
    State next;
    real error;
};

std::ostream& operator<<(std::ostream& out, Step step) {
    const auto [state, error] = step;
    const auto [t, dt, current] = state;
    const auto [x, xp] = current;
    return out << t << ' ' << x << ' ' << xp << ' ' << dt << ' ' << error << '\n';
}

Step one_step(State current, Vec advanced) {
    const auto [t, dt, curr] = current;
    const auto [x, xp] = curr;
    const real xpp = accel(x);
    // Calculate the new (x, xp) when advanced by dt/2.
    // Then advance that result another dt/2.
    // Compare this result with `advanced` by calculating an "error".
    // If the error is acceptable, return advanced.
    // If the error is too large, recurse with dt/2 (using some values we just calculated).
    // If the error is too small, return advanced but adjust dt -> 2*dt.
    const Vec half{
        .x = x + xp * dt/2,
        .xp = xp + xpp * dt/2
    };
    const Vec half_twice{
        .x = half.x + half.xp * dt/2,
        .xp = half.xp + accel(half.x) * dt/2
    };
    const real err = error(current.current, advanced, half_twice);
    constexpr real threshold = 1e-6;
    if (err > threshold) {
        return one_step({t, dt/2, curr}, half);
    }
    if (err < threshold/4) {
        return {{t+dt/2, dt*2, half}, err};
    }
    return {{t+dt/2, dt, half}, err};
}

void simulate(std::ostream& out, real min_x, State state) {
    // initial point: `state` with zero error by definition
    out << Step{.next=state, .error=0};

    auto& [t, dt, curr]  = state;
    auto& [x, xp] = curr;
    while (x >= min_x) {
        // TODO: store the accel value
        Vec advanced{
            .x = x + xp * dt,
            .xp = xp + accel(x) * dt
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
        .current = Vec{.x = 1.0, .xp = 0.0 }
    };
    simulate(std::cout, 1e-6, start);
}

/* more notes

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

/* notes

dr = (r_t + r_tt * dt) * dt = r_t * dt  +  r_tt * dt * dt

suppose there's a max_abs_dr > 0
Then we want to choose a dt such that abs(dr) is at most max_abs_dr

call dt {such that abs(dr) = max_abs_dr} max_dt > 0

max_abs_dr = abs[r_t * dt_max  +  r_tt * dt_max * dt_max]

quadratic. let's ignore abs for now.
a = r_tt
b = r_t
c = -max_abs_dr

dt_max = [-r_t +- sqrt(r_t * r_t  +  4 * r_tt * max_abs_dr)] / (2 * r_tt)

*/

