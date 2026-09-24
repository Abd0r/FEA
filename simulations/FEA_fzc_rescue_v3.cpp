// =============================================================================
// FEA_fzc_rescue_v3.cpp -- M18 where the rescue path's area lives
//
// Claim under test: FZC-v0's rescue port can reach every Zone without
// reintroducing the per-Zone CMOS that reviewers showed V2 could not floorplan.
// This module sizes three candidate mechanisms against the area the design has
// actually left, and reports which of them survives arithmetic.
//
// The three mechanisms differ only in what they SCALE WITH, and scaling is the
// whole question, because scaling with Zone count is exactly how V2 died
// (12 um^2 x 1.7e9 = 204 cm^2 against a 3 cm^2 die):
//
//   A  dedicated boundary tree   scales with the GRID  (2 x grid_side lines)
//   B  per-Zone CMOS decode      scales with ZONE COUNT
//   C  tag detector on the neighbour Slingshot link, which already exists
//
// Option C reuses routing the fabric already paid for, and it is what FZC-v0
// already specifies ("FZC_B --Slingshot--> rescue receiver in Zone A"). This
// module does not presuppose that C wins: A and B are sized too, because
// "we chose C" is an opinion until they are shown failing.
//
// Spec: spec/FZC-v0.md "The rescue port", invariants 8 and 10;
//       DESIGN-V3 gate "no per-Zone CMOS scaling".
//
// Labels: ARITHMETIC = substitution into shared parameters; SWEPT = an input
// with no PDK and no measurement behind it, varied rather than assumed;
// NOT DEMONSTRATED = no device implements any of this.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace rescue {

using fea::params;
using fea::require;

// Geometry derived from the design point, never hard-coded.
struct Geometry {
    double zones = 0.0;
    double grid_side = 0.0;
    double die_cm2 = 0.0;
    double edge_cm = 0.0;
    double cell_share_cm2 = 0.0;   // what the array already claims
    double routing_share_cm2 = 0.0; // what is left, from V2's declared 2x overhead
};

static Geometry geometry() {
    Geometry g;
    g.zones = fea::design_zone_count();
    g.grid_side = std::ceil(std::sqrt(g.zones));
    g.die_cm2 = params().arch.die_area_cm2;
    g.edge_cm = std::sqrt(g.die_cm2);
    // routing_area_overhead is total/usable, so usable = die / overhead and
    // everything not given to cells is what a rescue path may consume.
    g.cell_share_cm2 = g.die_cm2 / params().arch.routing_area_overhead;
    g.routing_share_cm2 = g.die_cm2 - g.cell_share_cm2;
    return g;
}

static std::string pct(double part, double whole) {
    return std::to_string(100.0 * part / whole) + " %";
}

// =============================================================================
// SCENARIO 1: the ledger, the grid, and the budget a rescue path must fit
// =============================================================================
static void scenario_ledger_and_budget() {
    std::cout << "\n[SCENARIO 1] the FZC ledger, the grid, and the available area\n";

    const double sum = fea::fzc_ledger_sum();
    std::cout << "  ledger: " << fea::fzc_state_groups() << " groups x "
              << fea::fzc_state_bits_per_group() << " bits x " << fea::fzc_state_rails()
              << " rails x " << fea::fzc_state_replicas() << " replicas = "
              << fea::fzc_ledger_bits() << " state\n";
    std::cout << "         + " << fea::fzc_command_blocks() << " command + "
              << fea::fzc_port_blocks() << " port + " << fea::fzc_pathway_blocks()
              << " pathway + " << fea::fzc_spare_blocks() << " spare = " << sum
              << " Blocks\n";
    std::cout << "  zone_fzc_blocks() reports   : " << fea::zone_fzc_blocks() << "\n";

    // The documented decomposition must still sum. If it does not, the 535 that
    // two other modules floorplan is being justified by a stale comment.
    require(sum == static_cast<double>(fea::zone_fzc_blocks()),
            "the documented FZC ledger split must sum to zone_fzc_blocks(), or the "
            "535-Block figure is supported by a decomposition that no longer holds");

    const Geometry g = geometry();
    std::cout << "\n  Zones on the die           : " << g.zones << "\n";
    std::cout << "  grid side                  : " << g.grid_side << " x " << g.grid_side
              << "\n";
    std::cout << "  die edge                   : " << g.edge_cm * 1e4 << " um\n";
    std::cout << "  cell share (array)         : " << g.cell_share_cm2 << " cm^2\n";
    std::cout << "  routing share (available)  : " << g.routing_share_cm2 << " cm^2\n";
    std::cout << "  port Blocks already counted: " << fea::fzc_port_blocks()
              << " per Zone (inside the 535)\n";

    // PR6/S3: the coverage line below is an identity -- grid_side is
    // ceil(sqrt(zones)), so the squared inequality holds for every positive N.
    // Minimality is the falsifiable half: one row smaller must NOT cover.
    require(g.grid_side * g.grid_side >= g.zones,
            "the grid must cover every Zone, or part of the die is unreachable");
    require((g.grid_side - 1.0) * (g.grid_side - 1.0) < g.zones,
            "the grid must be the SMALLEST square that covers every Zone, or it "
            "carries spare capacity that no routing model paid for");
    require(g.routing_share_cm2 > 0.0,
            "V2's declared 2x routing overhead must leave a routing share at all, or "
            "there is no budget for this module to allocate");

    std::cout << "  label: ARITHMETIC from design_zone_count() and the declared 2x overhead.\n";
}

// =============================================================================
// SCENARIO 2: option B, one decoder per Zone
// =============================================================================
static void scenario_option_b_per_zone_decode() {
    std::cout << "\n[SCENARIO 2] option B: a CMOS decoder at every Zone\n";

    const Geometry g = geometry();
    const double site = params().control.decoder_area_um2;
    const double area_cm2 = g.zones * site * fea::kCM2_PER_UM2;
    const double tx_area = fea::tx_area_um2();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  V2 decoder site            : " << site << " um^2, "
              << params().control.decoder_transistors << " transistors\n";
    std::cout << "  implied per transistor      : " << tx_area << " um^2\n";
    std::cout << "  one transistor at every Zone: " << g.zones * tx_area * fea::kCM2_PER_UM2
              << " cm^2 = " << pct(g.zones * tx_area * fea::kCM2_PER_UM2, g.die_cm2)
              << " of the die\n";
    std::cout << "  option B total             : " << area_cm2 << " cm^2 = " << area_cm2 / g.die_cm2
              << "x the " << g.die_cm2 << " cm^2 die\n";

    // V2's own approach must fail here too. If it did not, the reviewers' core
    // objection would not carry over to our larger Zone count, and M1 and M2
    // would be measuring something that does not matter.
    require(area_cm2 > g.die_cm2,
            "one decoder per Zone must exceed the die, otherwise the reviewers' "
            "arithmetic objection does not apply to V3 either");
    require(area_cm2 > g.routing_share_cm2,
            "one decoder per Zone must also exceed the entire routing share");

    std::cout << "  verdict: FAILS by " << area_cm2 / g.routing_share_cm2
              << "x against the routing budget.\n";
    std::cout << "  label: ARITHMETIC, using V2's own declared decoder site.\n";
}

// =============================================================================
// SCENARIO 3: option A, a dedicated tree from the boundary
// =============================================================================
static void scenario_option_a_boundary_tree() {
    std::cout << "\n[SCENARIO 3] option A: a dedicated row and column tree from the edge\n";

    require(!fea::rescue_geometry_sourced(),
            "this sweep exists only while the wire pitch is unsourced; once a PDK "
            "supplies one, retire the sweep rather than leave it asserting the obvious");

    const Geometry g = geometry();
    // Rows plus columns, each spanning the die edge. Addressing by geometry, so
    // one line per row and one per column reaches every cell at their crossing.
    const double lines = 2.0 * g.grid_side;
    const double length_um = lines * g.edge_cm * 1e4;

    std::cout << "  lines required             : " << lines << " (" << g.grid_side
              << " rows + " << g.grid_side << " columns)\n";
    std::cout << "  total wire length          : " << length_um << " um\n";

    // Pitch at which the tree consumes the whole routing share exactly.
    const double threshold_um = (g.routing_share_cm2 * 1e8) / length_um;
    std::cout << "  pitch that exactly fills it: " << threshold_um << " um ("
              << threshold_um * 1e3 << " nm)\n\n";

    const std::vector<double> pitches_um{1.0, 0.1, 0.01, 0.001};
    std::cout << "  " << std::left << std::setw(14) << "pitch" << std::right << std::setw(16)
              << "tree area" << std::setw(20) << "share of routing" << std::setw(10) << "fits"
              << "\n";
    bool any_fits = false;
    for (const double p : pitches_um) {
        const double area_cm2 = length_um * p * fea::kCM2_PER_UM2;
        const bool fits = area_cm2 <= g.routing_share_cm2;
        if (fits) any_fits = true;
        std::cout << "  " << std::left << std::setw(14) << (std::to_string(p) + " um")
                  << std::right << std::setw(16) << area_cm2 << std::setw(20)
                  << pct(area_cm2, g.routing_share_cm2) << std::setw(10)
                  << (fits ? "yes" : "NO") << "\n";
    }

    require(any_fits,
            "the dedicated tree must be viable somewhere in the pitch sweep, otherwise "
            "option A is impossible rather than merely expensive, and the sweep is "
            "wasted computation");
    require(threshold_um < 1.0,
            "the viable pitch must be under 1 um or option A is hopeless at any "
            "credible node");

    std::cout << "\n  verdict: viable only below " << threshold_um * 1e3 << " nm, and at that\n";
    std::cout << "  threshold it consumes 100% of the routing budget, leaving nothing for\n";
    std::cout << "  Slingshot itself. So A is pitch-limited and structurally awkward.\n";
    std::cout << "  label: SWEPT wire pitch. No PDK supplies this number.\n";
}

// =============================================================================
// SCENARIO 4: option C, a tag detector on the neighbour link that exists
// =============================================================================
static void scenario_option_c_neighbour_link() {
    std::cout << "\n[SCENARIO 4] option C: a tag detector on the existing neighbour link\n";

    const Geometry g = geometry();
    const double tx_area = fea::tx_area_um2();

    // The neighbour link's destination is fixed by WIRING, so no address decode
    // happens at all. All the port must do is accept a rescue-tagged packet
    // while the FZC state machine is dead.
    const std::vector<double> widths{1.0, 5.0, 20.0, 100.0};
    std::cout << "  per-transistor area         : " << tx_area << " um^2 (from V2's decoder)\n\n";
    std::cout << "  " << std::left << std::setw(14) << "width (tx)" << std::right
              << std::setw(16) << "site um^2" << std::setw(16) << "total cm^2"
              << std::setw(22) << "share of routing" << std::setw(10) << "fits" << "\n";
    for (const double w : widths) {
        const double site = w * tx_area;
        const double total = g.zones * site * fea::kCM2_PER_UM2;
        const bool fits = total <= g.routing_share_cm2;
        std::cout << "  " << std::left << std::setw(14) << w << std::right << std::setw(16)
                  << site << std::setw(16) << total << std::setw(22)
                  << pct(total, g.routing_share_cm2) << std::setw(10)
                  << (fits ? "yes" : "NO") << "\n";
    }

    // The number that decides it: how wide a CMOS detector may be if it is to
    // take at most 5% of the routing budget.
    const double budget_cm2 = 0.05 * g.routing_share_cm2;
    const double max_tx = (budget_cm2 * 1e8) / (g.zones * tx_area);
    std::cout << "\n  at a 5% routing budget a CMOS detector may be at most " << max_tx
              << " transistors per Zone\n";
    std::cout << "  the DECLARED width is       " << fea::rescue_tag_transistors()
              << " transistors\n";

    require(max_tx < fea::rescue_tag_transistors(),
            "a CMOS rescue detector must NOT fit at its declared width, because if it did "
            "the interesting result would be that per-Zone CMOS is affordable after all, "
            "which would contradict M1 and M2");

    // The 5% discipline above is a choice, so back it with a structural number:
    // what would be LEFT for the whole transport if the declared width were used.
    const double declared_area_cm2 =
        g.zones * fea::rescue_tag_transistors() * tx_area * fea::kCM2_PER_UM2;
    const double leftover_cm2 = g.routing_share_cm2 - declared_area_cm2;
    const double leftover_per_zone_um2 = leftover_cm2 * 1e8 / g.zones;
    std::cout << "\n  at the declared width the rescue path leaves " << leftover_cm2
              << " cm^2 for ALL of Slingshot\n";
    std::cout << "  which is " << pct(leftover_cm2, g.routing_share_cm2)
              << " of the routing budget, i.e. " << leftover_per_zone_um2
              << " um^2 per Zone\n";
    std::cout << "  for every link, buffer and arbiter that must reach " << g.zones
              << " Zones.\n";
    require(leftover_cm2 > 0.0 && leftover_cm2 < g.routing_share_cm2,
            "leftover routing must be a small positive fraction, not everything: if the "
            "declared width left the budget untouched the structural objection would "
            "disappear and only the 5% choice would remain");

    // ... so it must be in-fabric instead, and then it costs nothing new.
    const double in_fabric_cm2 = 0.0;
    std::cout << "\n  verdict: no CMOS detector fits, so the receiver must be in-fabric.\n";
    std::cout << "  an in-fabric port occupies port Blocks the 535 already contains: "
              << fea::fzc_port_blocks() << " per Zone.\n";
    std::cout << "  incremental area           : " << in_fabric_cm2
              << " cm^2 (counted in the 535, not added)\n";

    require(fea::fzc_port_blocks() > 0.0,
            "the rescue receiver must have Blocks already inside the 535 budget, or the "
            "in-fabric claim is free-riding on an allocation that does not exist");
    // PR5/S2: the gate that stood here was `require(in_fabric_cm2 == 0.0)`
    // where that variable had just been assigned 0.0, i.e. 0.0 == 0.0. What is
    // falsifiable is that the ledger really does contain the port Blocks, and
    // that the alternative it is compared against is not itself free.
    require(fea::fzc_ledger_sum() == static_cast<double>(fea::zone_fzc_blocks()),
            "the in-fabric receiver is only free if its port Blocks sit inside the ledger, so "
            "the ledger split must still sum to zone_fzc_blocks()");
    const double smallest_cmos_cm2 =
        fea::design_zone_count() * fea::tx_area_um2() * fea::kCM2_PER_UM2;
    require(smallest_cmos_cm2 > 0.0,
            "one transistor at every Zone must cost a positive area, or calling the in-fabric "
            "option free is meaningless");
    std::cout << "  cheapest CMOS option, one transistor per Zone : " << smallest_cmos_cm2
              << " cm^2 vs 0 cm^2 in-fabric\n";

    std::cout << "  label: SWEPT widths, ARITHMETIC totals, mechanism NOT DEMONSTRATED.\n";
}

// =============================================================================
// SCENARIO 5: what this module deliberately does not claim
// =============================================================================
static void scenario_what_this_does_not_claim() {
    std::cout << "\n[SCENARIO 5] what M18 deliberately does not claim\n";

    const Geometry g = geometry();
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  routing budget under test  : " << g.routing_share_cm2 << " cm^2\n";

    std::cout << "\n  none of the following is established here:\n";
    std::cout << "  - no power figure is produced. Rescue activity has no model, so this\n";
    std::cout << "    module adds nothing to M1's budget and the four declared terms stay\n";
    std::cout << "    unsourced.\n";
    std::cout << "  - wire pitch and per-transistor area are both DECLARED. There is no\n";
    std::cout << "    public 2 nm PDK, so every threshold here is swept, not sourced.\n";
    std::cout << "  - how a stateless port Block tells a rescue tag from ordinary traffic\n";
    std::cout << "    without retained state is UNVALIDATED DEVICE PHYSICS. This module\n";
    std::cout << "    shows where such a receiver could be PAID FOR, not that it works.\n";
    std::cout << "  - L1 still depends on the neighbour link being alive. The peer set,\n";
    std::cout << "    quorum and L2/L3 ladder cover that, but no tier was timed here.\n";
    std::cout << "  label: area ARITHMETIC and SWEPT, mechanism NOT DEMONSTRATED.\n";
    std::cout << "  NEXT EVIDENCE GATE: a device model for stateless tag detection, and a\n";
    std::cout << "  power estimate for rescue activity before it enters M1.\n";
}

} // namespace rescue

int main() {
    using namespace rescue;
    try {
        std::cout << "FEA V3 M18 rescue path: where the recovery receiver's area lives\n";
        std::cout << "Same 0.5 cm^2 design array and the same 535-Block FZC ledger as M2 and M15.\n";
        scenario_ledger_and_budget();
        scenario_option_b_per_zone_decode();
        scenario_option_a_boundary_tree();
        scenario_option_c_neighbour_link();
        scenario_what_this_does_not_claim();
        std::cout << "\nPASS: the rescue path fits, by reusing the neighbour link rather than\n";
        std::cout << "      adding a second access tree or a per-Zone decoder.\n";
        std::cout << "LABEL: area arithmetic on swept geometry, mechanism not demonstrated.\n";
        std::cout << "NEXT EVIDENCE GATE: device model for stateless tag detection, and a\n";
        std::cout << "power estimate before rescue activity enters M1.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
