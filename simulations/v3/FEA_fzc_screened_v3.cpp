// =============================================================================
// FEA_fzc_screened_v3.cpp -- screened selector bias and thermal escape
//
// Same declared four-DB geometry as FEA_fzc_selector_v3.cpp. Bare 1/r is
// replaced by a Yukawa form exp(-r/lambda)/r. lambda is a declared screening
// length in lattice units, not a measured silicon value.
//
// Stored occupation then decays as exp(-k t), with k = k0 exp(-E_bind/theta).
// E_bind and theta are declared parameters. This is not a retention time.
// =============================================================================

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace screened {

struct Vec2 { double x = 0.0; double y = 0.0; };

static const Vec2 kPlusA{0.5, 0.5};
static const Vec2 kPlusB{-0.5, -0.5};
static const Vec2 kMinusA{-0.5, 0.5};
static const Vec2 kMinusB{0.5, -0.5};
static const double kStoreMargin = 0.5;

static void require(bool ok, const std::string& what) {
    if (!ok) throw std::runtime_error("ASSERTION FAILED: " + what);
}

static double distance(Vec2 a, Vec2 b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

static double yukawa(Vec2 a, Vec2 b, double lambda) {
    const double r = distance(a, b);
    return std::exp(-r / lambda) / r;
}

static double bias(Vec2 command, double lambda) {
    const double plus = yukawa(command, kPlusA, lambda) + yukawa(command, kPlusB, lambda);
    const double minus = yukawa(command, kMinusA, lambda) + yukawa(command, kMinusB, lambda);
    return minus - plus;
}

static Vec2 command_site(double separation, bool upper) {
    return Vec2{-separation, upper ? 0.5 : -0.5};
}

static double clock_at(double amplitude, double time) {
    if (amplitude <= 0.0 || time < 0.0 || time > 1.0) return 0.0;
    return amplitude * std::sin(3.14159265358979323846 * time);
}

static double stored_after_clock(double bias_value, double amplitude) {
    double lead = 1.0;
    double actuator = 0.0;
    double stored = 0.0;
    const double dt = 0.001;
    for (int step = 0; step < 1500; ++step) {
        const double time = step * dt;
        const double clock = clock_at(amplitude, time);
        const double barrier = std::max(0.0, 1.5 - bias_value - clock);
        const double enter = dt * 2.0 * std::exp(-barrier) * lead;
        const double leak = dt * 0.04 * lead;
        const bool falling = clock_at(amplitude, time) > clock_at(amplitude, time + dt);
        const double relax = dt * ((bias_value > kStoreMargin && falling) ? 8.0 : 0.0) * std::max(actuator, 0.0);
        const double escape = dt * 0.05 * std::exp(-std::abs(bias_value)) * std::max(actuator, 0.0);
        lead -= enter + leak;
        actuator += enter - relax - escape;
        stored += relax;
    }
    return stored;
}

static double retained(double stored, double theta, double e_bind, double hold_time, double k0 = 1.0) {
    const double rate = k0 * std::exp(-e_bind / theta);
    return stored * std::exp(-rate * hold_time);
}

static void screening_changes_bias() {
    std::cout << "\n[SCENARIO 1] screening can erase a bias that bare coupling keeps\n";
    const Vec2 near = command_site(1.0, true);
    const double weak = bias(near, 100.0);
    const double strong = bias(near, 0.35);
    require(weak > kStoreMargin, "long screening length must keep the near-rail bias above the store margin");
    require(strong < kStoreMargin, "short screening length must drop the same geometry below the store margin");
    std::cout << "  lambda=100 bias=" << weak << ", lambda=0.35 bias=" << strong << "\n";
}

static void thermal_escape() {
    std::cout << "\n[SCENARIO 2] stored occupation escapes when theta rises\n";
    const double stored = stored_after_clock(bias(command_site(1.0, true), 100.0), 1.2);
    require(stored > kStoreMargin, "clocked near-rail geometry must compute a stored occupation before the hold");
    const double e_bind = 2.0;
    const double hold = 5.0;
    const double cold = retained(stored, 0.25, e_bind, hold);
    const double hot = retained(stored, 4.0, e_bind, hold);
    require(cold > kStoreMargin, "low theta must retain the prior valid stored occupation in this parameter set");
    require(hot < 0.1, "high theta must empty that occupation during the same hold");
    std::cout << "  theta=0.25 retained=" << cold << ", theta=4 retained=" << hot << "\n";
}

static void cascade_needs_margin() {
    std::cout << "\n[SCENARIO 3] a second stage fails when screened hop bias falls below margin\n";
    const double source = bias(command_site(1.0, true), 100.0);
    const double hop = 3.0;
    const double delivered_long = source * std::exp(-hop / 100.0);
    const double delivered_short = source * std::exp(-hop / 0.35);
    require(delivered_long > kStoreMargin, "weak screening must still deliver a usable stage-B bias");
    require(delivered_short < kStoreMargin, "strong screening must block the cascade");
    std::cout << "  hop=3 lambda=100 delivered=" << delivered_long << ", lambda=0.35 delivered=" << delivered_short << "\n";
}

} // namespace screened

int main() {
    using namespace screened;
    try {
        std::cout << "FEA V3 screened selector and thermal-escape model\n";
        std::cout << "Screening length and binding energy are declared parameters, not silicon measurements.\n";
        screening_changes_bias();
        thermal_escape();
        cascade_needs_margin();
        std::cout << "\nPASS: screening and thermal-escape gates held for this parameter set.\n";
        std::cout << "NEXT EVIDENCE GATE: replace lambda and E_bind with a screened Si Hamiltonian and a measured escape barrier.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
