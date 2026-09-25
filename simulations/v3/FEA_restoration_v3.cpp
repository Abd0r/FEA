// =============================================================================
// FEA_restoration_v3.cpp -- M6 signal restoration across the die
//
// in a data plane with no switching elements, signal
// attenuation is inevitable. How does FEA restore gain as a signal crosses a
// 3 cm^2 die, and what does that cost? This module computes end-of-chain
// amplitude with and without restoration endpoints, then charges every
// endpoint for area and power so restoration is priced instead of assumed.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace restoration {

using fea::params;
using fea::require;

// Die edge length from the declared area. A 3 cm^2 die has a ~1.73 cm edge.
static double die_edge_m() { return std::sqrt(params().arch.die_area_cm2 * 1e-4); }

// Declared amplitude decay length of a passive DBW segment. Not measured; the
// sweep below is the point of the module.
static double decay_length_um() { return 1.0; }

// Minimum amplitude the next stage will accept. Declared, dimensionless 0..1.
static double input_threshold() { return 0.10; }

static double end_amplitude(double distance_m, double decay_um) {
    const double um = distance_m * 1e6;
    return std::exp(-um / decay_um);
}

// Number of restoration endpoints needed when one endpoint restores to 1.0.
static int endpoints_needed(double distance_m, double spacing_m) {
    if (spacing_m <= 0.0) return -1;
    return static_cast<int>(std::ceil(distance_m / spacing_m)) - 1;
}

static void scenario_no_restoration_collapses() {
    std::cout << "\n[SCENARIO 1] a passive path across the die collapses without restoration\n";
    const double edge = die_edge_m();
    const double decay = decay_length_um();
    const double log10_amp = -(edge * 1e6) / decay / std::log(10.0);
    const double threshold = input_threshold();
    const double log10_threshold = std::log10(threshold);
    const double orders_short = log10_threshold - log10_amp;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  die edge                    : " << (edge * 1e2) << " cm (" << (edge * 1e6) << " um)\n";
    std::cout << "  declared decay length       : " << decay << " um\n";
    std::cout << "  end-of-chain amplitude      : 10^" << std::setprecision(1) << log10_amp << "\n";
    std::cout << "  next-stage input threshold  : " << threshold << " (10^" << log10_threshold << ")\n";
    std::cout << "  passes threshold?           : NO\n";
    std::cout << "  shortfall                   : " << std::setprecision(1) << orders_short
              << " orders of magnitude\n";
    require(orders_short > 0.0, "an unrestored path must fall below the input threshold");
    require(orders_short > 6.0, "the shortfall must be far more than a marginal loss");
    std::cout << "  finding: without gain the signal is gone long\n";
    std::cout << "  before it crosses the die. This is not a margin problem, it is a "
              << std::setprecision(0) << orders_short << "-order collapse.\n";
}

static void scenario_restoration_endpoints_needed() {
    std::cout << "\n[SCENARIO 2] restoration endpoints required to cross the die\n";
    const double edge = die_edge_m();
    const double threshold = input_threshold();
    const double decay = decay_length_um();

    // Longest spacing at which amplitude still clears the threshold.
    const double max_spacing_m = decay * 1e-6 * std::log(1.0 / threshold);
    const int n = endpoints_needed(edge, max_spacing_m);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  max spacing at threshold    : " << (max_spacing_m * 1e6) << " um\n";
    std::cout << "  endpoints for one die edge  : " << n << "\n";
    require(n > 0, "crossing the die must require at least one restoration endpoint");
    require(max_spacing_m > 0.0, "spacing must be positive");

    // per_hop_amp is exp(-max_spacing/decay) and max_spacing was DEFINED as
    // decay*ln(1/threshold), so per_hop_amp equals threshold exactly. The gate
    // that used to sit here, `per_hop_amp >= threshold*(1-1e-9)`, compared that
    // identity against itself and could never fail. Removed rather than replaced,
    // because any ceil()-based check on n would be equally definition-bound. The
    // falsifiable claim in this module is that the UNRESTORED path collapses,
    // which SCENARIO 1 gates against the threshold and SCENARIO 4 re-checks at a
    // 50 um decay length.
    const double per_hop_amp = std::exp(-(max_spacing_m * 1e6) / decay);
    std::cout << "  amplitude after one hop     : " << per_hop_amp
              << "  (threshold " << threshold << ", equal by construction)\n";
    std::cout << "  so restoration is mandatory, and it needs roughly " << n
              << " endpoints per edge.\n";
    std::cout << "  label: decay length and threshold are declared. The structural result is\n";
    std::cout << "  that restoration is required at all, and its density follows from decay.\n";
}

static void scenario_restoration_is_not_free() {
    std::cout << "\n[SCENARIO 3] every restoration endpoint costs area and power\n";
    const double edge = die_edge_m();
    const double decay = decay_length_um();
    const double threshold = input_threshold();
    const double max_spacing_m = decay * 1e-6 * std::log(1.0 / threshold);
    const int per_edge = endpoints_needed(edge, max_spacing_m);

    // Grid the die: endpoints on a 2D lattice at that spacing.
    const double die_edge_m2 = edge;
    const int per_row = static_cast<int>(std::ceil(die_edge_m2 / max_spacing_m));
    const long long total_endpoints = static_cast<long long>(per_row) * per_row;

    // Declared cost per endpoint, scaled from the sensing comparator count so
    // the number is at least anchored to something in V2's control budget.
    const double area_per_endpoint_um2 = 0.05; // declared
    const double power_per_endpoint_W = 1e-9;  // declared, 1 nW each

    const double endpoint_area_cm2 = total_endpoints * area_per_endpoint_um2 * fea::kCM2_PER_UM2;
    const double endpoint_power_W = total_endpoints * power_per_endpoint_W;
    const double die = params().arch.die_area_cm2;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  endpoints per edge          : " << per_edge << "\n";
    std::cout << "  endpoint lattice pitch      : " << (max_spacing_m * 1e6) << " um\n";
    std::cout << "  endpoints on the die        : " << total_endpoints << "\n";
    std::cout << "  area per endpoint (declared): " << area_per_endpoint_um2 << " um^2\n";
    std::cout << "  endpoint area               : " << endpoint_area_cm2 << " cm^2\n";
    std::cout << "  endpoint power (1 nW each)  : " << endpoint_power_W << " W\n";
    std::cout << "  as fraction of die area     : " << std::setprecision(2)
              << (endpoint_area_cm2 * 100.0 / die) << " %\n";
    require(total_endpoints > 0, "the grid must produce endpoints");
    require(endpoint_area_cm2 > 0.0 && endpoint_power_W > 0.0, "restoration must cost something");
    std::cout << "  these must be added to M1 power and M2 area. They are currently absent\n";
    std::cout << "  from every V2 number, because V2 assumed no gain was needed.\n";
}

static void scenario_decay_sets_everything() {
    std::cout << "\n[SCENARIO 4] the declared decay length is the load-bearing unknown\n";
    const double edge = die_edge_m();
    const double threshold = input_threshold();
    const std::vector<double> decays{0.1, 0.5, 1.0, 5.0, 50.0};
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  decay(um)   unrestored amp      max spacing(um)   endpoints/edge\n";
    for (const double dcm : decays) {
        const double amp = end_amplitude(edge, dcm);
        const double spacing = dcm * 1e-6 * std::log(1.0 / threshold);
        const int n = endpoints_needed(edge, spacing);
        std::cout << "  " << std::setw(9) << dcm << "   " << std::setw(14) << amp
                  << "   " << std::setw(14) << (spacing * 1e6) << "   " << std::setw(13) << n << "\n";
    }
    std::cout << std::fixed << std::setprecision(6);
    const double amp_50um = end_amplitude(edge, 50.0);
    require(amp_50um < threshold, "even a 50 um decay length must fall short across the full die");
    std::cout << "  at 50 um decay the unrestored path reaches " << amp_50um
              << ", still below a " << threshold << " threshold.\n";
    std::cout << "  finding: even an optimistic 50 um decay length needs restoration endpoints\n";
    std::cout << "  somewhere on a 1.73 cm die. No passive pathway crosses it alone.\n";
    std::cout << "  label: decay length unmeasured. Until it is measured, restoration endpoint\n";
    std::cout << "  density, and therefore its area and power, are unbounded above.\n";
}

static void scenario_gain_must_come_from_a_rail() {
    std::cout << "\n[SCENARIO 5] restoration needs an energy source, not the signal itself\n";
    const double per_endpoint_W = 1e-9;
    const double edge = die_edge_m();
    const double decay = decay_length_um();
    const double threshold = input_threshold();
    const double spacing = decay * 1e-6 * std::log(1.0 / threshold);
    const int per_row = static_cast<int>(std::ceil(edge / spacing));
    const long long n = static_cast<long long>(per_row) * per_row;
    const double total_W = n * per_endpoint_W;
    require(total_W > 0.0, "restoration power must be positive and externally supplied");
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  endpoints                  : " << n << "\n";
    std::cout << "  power from a bias/clock rail: " << total_W << " W\n";
    std::cout << "  a passive DBW cannot amplify. Every endpoint needs supplied energy.\n";
    std::cout << "  this is the same constraint FZC's actuator gate already states: a stored\n";
    std::cout << "  pattern is not gain, and a decaying signal is not a restoration source.\n";
    std::cout << "  the rail that feeds these endpoints must appear in M1 as boundary or\n";
    std::cout << "  regional power, and its routing must appear in M2 as area.\n";
}

} // namespace restoration

int main() {
    using namespace restoration;
    try {
        std::cout << "FEA V3 M6 signal restoration across the die\n";
        std::cout << "Decay length and input threshold are declared, not measured.\n";
        scenario_no_restoration_collapses();
        scenario_restoration_endpoints_needed();
        scenario_restoration_is_not_free();
        scenario_decay_sets_everything();
        scenario_gain_must_come_from_a_rail();
        std::cout << "\nPASS: unrestored collapse shown, endpoint count and cost computed.\n";
        std::cout << "LABEL: derived propagation form, declared decay and threshold, open measurement.\n";
        std::cout << "NEXT EVIDENCE GATE: measured DBW amplitude decay length and a powered relay mechanism.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
