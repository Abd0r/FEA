// =============================================================================
// FEA_fzc_selector_v3.cpp -- stated four-DB selector geometry and rate model
//
// Geometry (declared lattice units, not nanometres):
//   selector square sites at (+0.5,+0.5), (-0.5,+0.5), (-0.5,-0.5), (+0.5,-0.5)
//   P+ occupies the ( +0.5,+0.5) / (-0.5,-0.5) diagonal
//   P- occupies the ( -0.5,+0.5) / (+0.5,-0.5) diagonal
//   one command electron sits on a dual rail at (-d, +0.5) or (-d, -0.5)
//
// The bias is computed from 1/r Coulomb sums. Haider 2009 reports coupling
// changes below about 2 nm; this file does not map lattice units to nm.
//
// Dynamics: a six-outcome master equation. Clock energy is an explicit
// waveform. A sensor reads stored occupation with a declared error rate.
// Parameters are uncalibrated. No timing, energy, or retention claim follows.
// =============================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace selector {

struct Vec2 {
    double x = 0.0;
    double y = 0.0;
};

struct Outcome {
    double lead = 1.0;
    double reflected = 0.0;
    double transmitted = 0.0;
    double actuator = 0.0;
    double stored = 0.0;
    double reservoir = 0.0;
    double reservoir_energy = 0.0;

    double probability() const { return lead + reflected + transmitted + actuator + stored + reservoir; }
};

struct Pulse {
    double amplitude = 0.0;
    double duration = 1.0;
};

static void require(bool ok, const std::string& what) {
    if (!ok) throw std::runtime_error("ASSERTION FAILED: " + what);
}

static double distance(Vec2 a, Vec2 b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

static double coulomb(Vec2 a, Vec2 b) {
    return 1.0 / distance(a, b);
}

static const Vec2 kPlusA{0.5, 0.5};
static const Vec2 kPlusB{-0.5, -0.5};
static const Vec2 kMinusA{-0.5, 0.5};
static const Vec2 kMinusB{0.5, -0.5};

static double polarization_energy(Vec2 command, bool plus) {
    if (plus) return coulomb(command, kPlusA) + coulomb(command, kPlusB);
    return coulomb(command, kMinusA) + coulomb(command, kMinusB);
}

// Positive means the command electron lowers P+ relative to P-.
static double polarization_bias(Vec2 command) {
    return polarization_energy(command, false) - polarization_energy(command, true);
}

static Vec2 command_site(double separation, bool upper_rail) {
    return Vec2{-separation, upper_rail ? 0.5 : -0.5};
}

static double clock_value(const Pulse& pulse, double time) {
    if (pulse.amplitude <= 0.0 || time < 0.0 || time > pulse.duration) return 0.0;
    const double phase = 3.14159265358979323846 * time / pulse.duration;
    return pulse.amplitude * std::sin(phase);
}

static bool clock_falling(const Pulse& pulse, double time, double dt) {
    return clock_value(pulse, time) > clock_value(pulse, time + dt);
}

struct RateParams {
    double theta = 1.0;
    double barrier0 = 1.5;
    double k0 = 2.0;
    double k_leak = 0.02;
    double k_relax = 8.0;
    double k_escape = 0.05;
    double store_margin = 0.5;
    double dt = 0.001;
    int steps = 1500;
};

static Outcome integrate(double bias, const Pulse& pulse, const RateParams& params) {
    Outcome state;
    for (int step = 0; step < params.steps; ++step) {
        const double time = step * params.dt;
        const double clock = clock_value(pulse, time);
        const double barrier = std::max(0.0, params.barrier0 - bias - clock);
        const double k_enter = params.k0 * std::exp(-barrier / params.theta);
        const double k_reflect = params.k_leak;
        const double k_transmit = params.k_leak;
        const bool favored = bias > params.store_margin && clock_falling(pulse, time, params.dt);
        const double k_relax = favored ? params.k_relax : 0.0;
        const double k_escape = params.k_escape * std::exp(-std::abs(bias) / params.theta);

        const double leave_lead = params.dt * (k_enter + k_reflect + k_transmit) * state.lead;
        const double enter = params.dt * k_enter * state.lead;
        const double reflect = params.dt * k_reflect * state.lead;
        const double transmit = params.dt * k_transmit * state.lead;
        const double relax = params.dt * k_relax * std::max(state.actuator, 0.0);
        const double escape = params.dt * k_escape * std::max(state.actuator, 0.0);

        state.lead -= leave_lead;
        state.actuator += enter - relax - escape;
        state.reflected += reflect;
        state.transmitted += transmit;
        state.stored += relax;
        state.reservoir += escape;
        state.reservoir_energy += relax + escape;
    }
    return state;
}

static bool sensor_confirms(const Outcome& state, double threshold, double error_rate, bool flip) {
    const bool raw = state.stored >= threshold;
    if (error_rate <= 0.0) return raw;
    return flip ? !raw : raw;
}

static void geometry_separation() {
    std::cout << "\n[SCENARIO 1] stated geometry: separation changes selector bias\n";
    const double near_bias = polarization_bias(command_site(1.0, true));
    const double far_bias = polarization_bias(command_site(6.0, true));
    const double opposite = polarization_bias(command_site(1.0, false));
    require(near_bias > 0.5, "near upper-rail command must favor P+");
    require(opposite < -0.5, "near lower-rail command must favor P-");
    require(std::abs(far_bias) < 0.5 * std::abs(near_bias), "far command must weaken the bias");
    std::cout << "  d=1 upper bias=" << near_bias << ", lower bias=" << opposite << ", d=6 bias=" << far_bias << "\n";
}

static void clocked_valid_command() {
    std::cout << "\n[SCENARIO 2] valid rail plus falling clock can store; invalid rail cannot\n";
    RateParams params;
    Pulse pulse{1.2, 1.0};
    const double valid_bias = polarization_bias(command_site(1.0, true));
    const double invalid_bias = polarization_bias(command_site(1.0, false));
    const Outcome valid = integrate(valid_bias, pulse, params);
    const Outcome invalid = integrate(invalid_bias, pulse, params);
    require(std::abs(valid.probability() - 1.0) < 1e-6, "valid trajectory must conserve probability");
    require(std::abs(invalid.probability() - 1.0) < 1e-6, "invalid trajectory must conserve probability");
    require(valid.stored > 0.5 && valid.reservoir_energy > 0.0, "valid command must reach stored occupation with reservoir exchange");
    require(invalid.stored < 0.1, "opposite rail must not store under the same clock");
    require(sensor_confirms(valid, 0.5, 0.0, false), "ideal sensor confirms the valid stored state");
    require(!sensor_confirms(invalid, 0.5, 0.0, false), "ideal sensor rejects the invalid rail");
    std::cout << "  valid stored=" << valid.stored << " reservoir_energy=" << valid.reservoir_energy
              << " invalid stored=" << invalid.stored << "\n";
}

static void clock_off_blocks_storage() {
    std::cout << "\n[SCENARIO 3] the same valid geometry stores nothing without clock energy\n";
    RateParams params;
    Pulse off{0.0, 1.0};
    const Outcome dark = integrate(polarization_bias(command_site(1.0, true)), off, params);
    require(std::abs(dark.probability() - 1.0) < 1e-6, "clock-off trajectory must conserve probability");
    require(dark.stored < 0.1, "stored pattern bias without clock energy must not restore");
    require(!sensor_confirms(dark, 0.5, 0.0, false), "sensor must not confirm a clock-off run");
    std::cout << "  clock amplitude=0 stored=" << dark.stored << " transmitted=" << dark.transmitted << "\n";
}

static void sensor_error_and_far_geometry() {
    std::cout << "\n[SCENARIO 4] sensor flip can change CONFIRM; far geometry fails storage\n";
    RateParams hot = {};
    hot.theta = 8.0;
    Pulse pulse{1.2, 1.0};
    const Outcome noisy = integrate(polarization_bias(command_site(1.0, true)), pulse, hot);
    const Outcome far = integrate(polarization_bias(command_site(6.0, true)), pulse, {});
    const bool raw = sensor_confirms(noisy, 0.5, 0.0, false);
    const bool flipped = sensor_confirms(noisy, 0.5, 0.1, true);
    require(raw, "ideal sensor must confirm the computed stored occupation");
    require(!flipped, "declared sensor flip must reject that same occupation");
    require(far.stored < 0.5, "far geometry must fail the stored-occupation gate under this parameter set");
    std::cout << "  hot theta=8 stored=" << noisy.stored << ", far d=6 stored=" << far.stored << "\n";
}

} // namespace selector

int main() {
    using namespace selector;
    try {
        std::cout << "FEA V3 stated four-DB selector model\n";
        std::cout << "Lattice units, rates, clock amplitude, and sensor error are declared parameters.\n";
        geometry_separation();
        clocked_valid_command();
        clock_off_blocks_storage();
        sensor_error_and_far_geometry();
        std::cout << "\nPASS: geometry, clock, and sensor gates held for this parameter set.\n";
        std::cout << "NEXT EVIDENCE GATE: map lattice units to a screened Si geometry and replace 1/r rates with a calibrated Hamiltonian.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
