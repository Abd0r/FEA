// =============================================================================
// FEA_multifire_v3.cpp -- M7 single-pass capture and correlated multi-FIRE
//
// Reviewers 5, 6 and 7 all flag the same gap: V2 reports ideal on-resonance
// absorption near unity, then uses a single-pass capture P_abs about 0.46 from
// the wavepacket simulation, then assumes successive FIRE pulses are
// independent Bernoulli trials to reach 99.99% fidelity. That independence is
// an assumption, not a measured result. This module computes the independent
// result, then adds a common-mode block fraction and shows that correlation
// places a hard ceiling no number of retries can cross.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace multifire {

using fea::check_unit;
using fea::params;
using fea::require;
using fea::p_abs_single;

// Single-pass capture, taken from V2's wavepacket result rather than the ideal
// A(E_F) = 1. V2 never derived the gap quantitatively; Reviewer CF14 #5 and
// Reviewer 7 #5 both demand it be explained. p itself comes from fea_params,
// shared with M9, M10, M11 and M13 through one definition.

// Independent Bernoulli: probability at least one of N attempts captures.
static double independent_success(int n, double p) {
    return 1.0 - std::pow(1.0 - p, n);
}

// Correlated model: a common-mode fraction f of attempts is blocked for reasons
// shared across pulses (a stuck actuator, a permanently off-resonance Block, a
// neighbour that pins the state). Retries cannot clear a common-mode block.
// Conditional on not being common-mode blocked, attempts are independent.
static double correlated_success(int n, double p, double common_mode_fraction) {
    const double available = 1.0 - common_mode_fraction;
    return available * independent_success(n, p);
}

static int fires_for_target(double target, double p, double common_mode_fraction = 0.0) {
    for (int n = 1; n <= 10000; ++n) {
        if (correlated_success(n, p, common_mode_fraction) >= target) return n;
    }
    return -1;
}

static void scenario_single_pass_gap() {
    std::cout << "\n[SCENARIO 1] ideal absorption and simulated single-pass capture disagree\n";
    const double ideal = 1.0;      // A(E_F) at exact resonance, V2 Section 3
    const double simulated = p_abs_single(); // wavepacket SIM 4 result, V2 Section 3
    const double gamma = params().device.gamma_v1_hardcoded_meV;
    (void)gamma;
    // ideal and simulated are V2's OWN two values (Section 3: A(E_F)=1 ideal,
    // SIM 4 wavepacket P_abs=0.46). The three requires that used to sit here
    // gated those literals (ideal > 0.99, 0.40 < p < 0.50, gap > 0.5) and could
    // not fail unless someone edited a literal, so they are removed. What is
    // actually COMPUTED in this module is downstream: fires needed to close the
    // gap, the correlation ceiling, and the implied program pass rate.
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  ideal A(E_F)              : " << ideal << "\n";
    std::cout << "  wavepacket single pass    : " << simulated << "\n";
    std::cout << "  unexplained gap           : " << (ideal - simulated)
              << " (" << std::setprecision(1) << (ideal / simulated * 100.0) << "% of ideal)\n";
    std::cout << "  V2 attributed this to unspecified finite-pulse effects with no mechanism.\n";
    std::cout << "  label: both values are V2's own. The gap stays open until a finite-pulse\n";
    std::cout << "  calculation derives 0.46 from the same Hamiltonian that gives A(E_F) = 1.\n";
}

static void scenario_independent_redundancy() {
    std::cout << "\n[SCENARIO 2] independent-Bernoulli redundancy reproduces V2's 18-fire figure\n";
    const double p = p_abs_single();
    const double target_9999 = 0.9999;
    const int n = fires_for_target(target_9999, p);
    const double at_18 = independent_success(18, p);
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  P_abs = " << p << ", target " << target_9999 << "\n";
    std::cout << "  minimum fires for target     : " << n << "\n";
    std::cout << "  success at 18 fires          : " << at_18 << "\n";
    std::cout << "  success at 1 fire            : " << independent_success(1, p) << "\n";
    require(at_18 >= target_9999,
            "V2's 18-fold redundancy must land above its own 99.99% write-fidelity claim");
    require(n > 0 && n < 18,
            "the true minimum must be below 18 so V2's stated redundancy can be checked against it");
    std::cout << "  V2's 18 fires reach " << std::setprecision(4) << at_18
              << ", which covers its 99.99% claim.\n";
    std::cout << "  but " << std::setprecision(0) << n
              << " fires already suffice, so 18 is conservative or came from a different target.\n";
    std::cout << "  either way the whole result rests on the independence assumption.\n";
}

static void scenario_correlation_ceiling() {
    std::cout << "\n[SCENARIO 3] a common-mode block fraction puts a hard ceiling on fidelity\n";
    const double p = p_abs_single();
    std::cout << std::fixed << std::setprecision(4);
    const std::vector<double> fractions{0.0, 0.01, 0.02, 0.05, 0.10};
    for (const double f : fractions) {
        const double many = correlated_success(10000, p, f);
        const int need = fires_for_target(0.9999, p, f);
        std::cout << "  common-mode fraction " << std::setprecision(2) << f
                  << std::setprecision(4) << "  ->  ceiling " << many;
        if (need > 0) std::cout << "   99.99% reachable in " << need << " fires";
        else std::cout << "   99.99% UNREACHABLE at any fire count";
        std::cout << "\n";
    }
    const double ceiling_5pct = correlated_success(10000, p, 0.05);
    require(ceiling_5pct < 0.96, "a 5% common-mode fraction must cap fidelity below 96%");
    require(fires_for_target(0.9999, p, 0.05) < 0, "99.99% must be unreachable once common-mode blocking exists");
    std::cout << "  a 5% shared failure fraction caps fidelity at " << std::setprecision(1)
              << (ceiling_5pct * 100.0) << "%, and no number of retries crosses it.\n";
    std::cout << "  V2's 99.99% claim is therefore conditional on exactly zero correlation,\n";
    std::cout << "  which has not been measured for the 5-atom cluster.\n";
}

static void scenario_program_pass_rate() {
    std::cout << "\n[SCENARIO 4] program pass rate must follow from the capture model\n";
    const double p = p_abs_single();
    // V2 reports 96% for a 16-element vector-add and 98% for a 16-element dot
    // product. Under a per-operation capture success rate s, an m-operation
    // program passes with s^m. Solve for the s each figure implies.
    const double vector_ops = 16.0;
    const double dot_ops = 16.0;
    const double s_vector = std::pow(0.96, 1.0 / vector_ops);
    const double s_dot = std::pow(0.98, 1.0 / dot_ops);
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "  V2 program pass: 96% vector-add, 98% dot-product, both 16 elements\n";
    std::cout << "  implied per-op success for 96% over 16 ops : " << s_vector << "\n";
    std::cout << "  implied per-op success for 98% over 16 ops : " << s_dot << "\n";

    // Now ask what multi-FIRE gives, independently and with correlation.
    const int fires_used = 18;
    const double s_independent = independent_success(fires_used, p);
    const double s_corr = correlated_success(fires_used, p, 0.02);
    std::cout << "  model per-op success at 18 fires, independent : " << s_independent << "\n";
    std::cout << "  model per-op success at 18 fires, 2% correlated: " << s_corr << "\n";
    std::cout << "  program pass at 16 ops, independent            : "
              << std::setprecision(4) << std::pow(s_independent, 16) << "\n";
    std::cout << "  program pass at 16 ops, 2% correlated           : "
              << std::pow(s_corr, 16) << "\n";

    // The reported rates are far below what 18 independent fires predicts.
    require(s_vector < s_independent,
            "V2's implied per-op success must fall short of its own 18-fire independent model");
    require(s_dot < s_independent,
            "the dot-product implied success must also fall short of the independent model");
    std::cout << "  finding: V2's reported 96-98% program rates are far WORSE than its own\n";
    std::cout << "  18-fire independent model predicts, yet the paper never links the two.\n";
    std::cout << "  either the program runs use fewer fires than 18, or another error source\n";
    std::cout << "  exists that the capture model does not contain. Reviewer 6 #7 is correct:\n";
    std::cout << "  program success and SECDED are asserted, not derived from one model.\n";
}

static void scenario_latent_costs() {
    std::cout << "\n[SCENARIO 5] multi-FIRE latency is not a single FIRE latency\n";
    const double t_fire_ps = params().timing.t_fire_ps;
    const double p = p_abs_single();
    // Expected fires to first success under independence.
    const double expected_fires = 1.0 / p;
    const double single_pass_latency = t_fire_ps;
    const double multi_latency = expected_fires * t_fire_ps;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  single-pass latency (1 FIRE)      : " << single_pass_latency << " ps\n";
    std::cout << "  expected fires to first success   : " << expected_fires << "\n";
    std::cout << "  expected write latency            : " << multi_latency << " ps\n";
    std::cout << "  V2 headline single-FIRE 0.87 ns figure is the idealised case only.\n";
    require(multi_latency > 2.0 * single_pass_latency,
            "expected multi-FIRE latency must exceed the single-FIRE figure");
    std::cout << "  the cycle time in M9 must use the expected latency, not the ideal one.\n";
}

} // namespace multifire

int main() {
    using namespace multifire;
    try {
        std::cout << "FEA V3 M7 multi-FIRE reliability\n";
        std::cout << "Independent Bernoulli reproduced first, then a declared common-mode fraction added.\n";
        scenario_single_pass_gap();
        scenario_independent_redundancy();
        scenario_correlation_ceiling();
        scenario_program_pass_rate();
        scenario_latent_costs();
        std::cout << "\nPASS: independence reproduced, correlation ceiling exposed, program rates unreconciled.\n";
        std::cout << "LABEL: derived from V2's own numbers, estimated correlation, open physical validation.\n";
        std::cout << "NEXT EVIDENCE GATE: measure attempt-to-attempt correlation for repeated FIRE on one Block.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
