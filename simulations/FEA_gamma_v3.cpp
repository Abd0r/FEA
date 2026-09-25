// =============================================================================
// FEA_gamma_v3.cpp -- M3 one linewidth derivation, then everything uses it
//
// Reviewers 5, 6 and 7 all found that V2 derived Gamma ~ 45 meV from the lead
// self-energy, then wrote that absorption used Gamma = 8 meV. This module
// derives Gamma once from the stated hoppings, uses that single value for
// absorption, thermal capture, and contrast, and fails if the two V2 literals
// are allowed to coexist as if they were the same quantity.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace gamma_model {

using fea::check_unit;
using fea::params;
using fea::require;

// Gamma from the lead self-energy Sigma = tc^2 g_L(E_F), following V2's stated
// form Gamma = -2 Im[Sigma] = 4 tc^2 / (2 t) at the band centre of the
// semi-infinite 1D DBW lead (hopping t). The cross has two leads, so the
// two-lead linewidth is twice the one-lead value.
static double gamma_one_lead_meV() {
    const auto& d = params().device;
    const double tc_eV = d.t_cluster_eV;
    const double t_eV = d.t_hop_eV;
    const double g_eV = (4.0 * tc_eV * tc_eV) / (2.0 * t_eV);
    return g_eV / fea::kMebitEV;
}

static double gamma_two_lead_meV() { return 2.0 * gamma_one_lead_meV(); }

// Breit-Wigner absorption on resonance and detuned, at a stated TOTAL Gamma
// (Gamma = Gamma_L + Gamma_R, as derived from the lead self-energy above).
//
// Convention: for a single resonant level the Lorentzian half-width at half
// maximum is Gamma/2, not Gamma. The denominator therefore carries
// (Gamma/2)^2, which makes A(E0) = 1 on resonance and gives HWHM = Gamma/2
// = 22.5 meV for the two-lead value of 45 meV. Using Gamma^2 here instead
// would double the line width and was the source of a factor-of-two error
// reported in external review (Jauho, Wingreen and Meir, cond-mat/9404027).
static double absorption(double detuning_meV, double gamma_meV) {
    const double half = 0.5 * gamma_meV;          // HWHM
    return (half * half) / (detuning_meV * detuning_meV + half * half);
}

// Thermal average of absorption over a Fermi window at temperature T.
static double thermal_absorption(double gamma_meV, double temperature_K) {
    const auto& d = params().device;
    const double kT_meV = (fea::kB * temperature_K) / fea::kEV / fea::kMebitEV;
    const double half_window_meV = 100.0;
    const int steps = 20001;
    double sum = 0.0;
    double weight_sum = 0.0;
    for (int i = 0; i < steps; ++i) {
        const double e = -half_window_meV + 2.0 * half_window_meV * i / (steps - 1);
        // Fermi weighting around E_F; normalized rather than absolute occupation.
        const double x = e / kT_meV;
        const double weight = 1.0 / (1.0 + std::exp(x));
        sum += weight * absorption(e, gamma_meV);
        weight_sum += weight;
    }
    (void)d;
    return weight_sum > 0.0 ? sum / weight_sum : 0.0;
}

// Off-resonance contrast: state 0 is detuned by the gate swing.
static double off_state_absorption(double gamma_meV) {
    return absorption(300.0, gamma_meV);
}

static void scenario_single_derivation() {
    std::cout << "\n[SCENARIO 1] Gamma is derived once from the stated hoppings\n";
    const double one = gamma_one_lead_meV();
    const double two = gamma_two_lead_meV();
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  tc = " << params().device.t_cluster_eV * 1e3 << " meV, t = "
              << params().device.t_hop_eV * 1e3 << " meV\n";
    std::cout << "  one-lead Gamma  = 4tc^2/(2t) = " << one << " meV\n";
    std::cout << "  two-lead Gamma  = 2 x one     = " << two << " meV\n";
    check_unit("derived one-lead Gamma", one, "meV", params().device.gamma_lead_meV, 0.05);
    check_unit("derived two-lead Gamma", two, "meV", params().device.gamma_two_lead_meV, 0.05);
    require(std::abs(two - 45.0) < 1.0, "two-lead Gamma must reproduce V2's derived 45 meV");
    std::cout << "  both of V2's stated numbers are reproducible: 22.5 meV per lead, 45 meV two-lead.\n";
}

static void scenario_v2_literals_are_different_quantities() {
    std::cout << "\n[SCENARIO 2] 8 meV and 45 meV cannot both be 'the' linewidth\n";
    const double v2_hardcoded = params().device.gamma_v1_hardcoded_meV;
    const double derived = gamma_two_lead_meV();
    require(std::abs(v2_hardcoded - derived) > 10.0,
            "the v1 hardcoded value and the derived value must be distinguishable");
    std::cout << "  v1 hardcoded Gamma = " << v2_hardcoded << " meV\n";
    std::cout << "  derived Gamma      = " << derived << " meV  (factor "
              << std::setprecision(1) << (derived / v2_hardcoded) << " apart)\n";
    std::cout << "  V2 Methods used 8 meV for absorption while Results claimed 45 meV for numerics.\n";
    std::cout << "  V3 rule: one derivation, one value, and the sensitivity of every downstream\n";
    std::cout << "  result to Gamma must be printed, not hidden.\n";
}

static void scenario_capture_depends_on_gamma() {
    std::cout << "\n[SCENARIO 3] capture probability is a strong function of the chosen Gamma\n";
    const double derived = gamma_two_lead_meV();
    const double v1 = params().device.gamma_v1_hardcoded_meV;
    const std::vector<double> sweep{v1, 16.0, 22.5, derived, 60.0};
    const double on_derived = thermal_absorption(derived, 300.0);
    const double on_v1 = thermal_absorption(v1, 300.0);
    std::cout << std::fixed << std::setprecision(4);
    for (const double g : sweep) {
        std::cout << "  Gamma=" << std::setprecision(1) << g << std::setprecision(4)
                  << " meV  thermal <A>=" << thermal_absorption(g, 300.0)
                  << "  off-state A(300 meV)=" << off_state_absorption(g) << "\n";
    }
    require(std::abs(on_derived - on_v1) > 0.01,
            "thermal absorption must actually change between the two V2 Gamma values");
    const double contrast_derived = off_state_absorption(derived) > 0
                                        ? thermal_absorption(derived, 300.0) / off_state_absorption(derived)
                                        : 0.0;
    std::cout << "  on/off contrast at derived Gamma = " << std::setprecision(1)
              << contrast_derived << "x\n";
    std::cout << "  conclusion: V2's headline capture figure is not Gamma-invariant,\n";
    std::cout << "  so the manuscript must pick one Gamma and report the sweep with it.\n";
}

} // namespace gamma_model

int main() {
    using namespace gamma_model;
    try {
        std::cout << "FEA V3 M3 linewidth consistency\n";
        std::cout << "Gamma derived once from tc and t. Not a fit to a preferred answer.\n";
        scenario_single_derivation();
        scenario_v2_literals_are_different_quantities();
        scenario_capture_depends_on_gamma();
        std::cout << "\nPASS: one linewidth derivation holds, and the V2 literals are exposed as distinct.\n";
        std::cout << "LABEL: derived Gamma, estimated hopping parameters, open capture calibration.\n";
        std::cout << "NEXT EVIDENCE GATE: replace the effective-mass lead Green's function with a real DBW band structure.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
