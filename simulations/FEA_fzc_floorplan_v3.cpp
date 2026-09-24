// =============================================================================
// FEA_fzc_floorplan_v3.cpp -- FZC block-budget ledger
//
// Counts Blocks required by the v0 contract. It does not draw a layout and it
// does not convert Blocks into nanometres, cm2, watts, or terabytes.
// Every width below is a declared allocation parameter.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace floorplan {

// Zone dimensions come from fea_params so this ledger cannot drift from the
// Zone restatement (65536 data + 535 FZC). They were local copies of 65536/256.
static const int kZoneBlocks = static_cast<int>(fea::zone_data_blocks());
static const int kZoneSide = static_cast<int>(fea::params().arch.zone_side_blocks);
static constexpr int kProvisionalFzc = 512;

struct Budget {
    int state_groups = 7;
    int bits_per_group = 8;
    int rails = 2;
    int replicas = 3;
    int commands = 10;
    int opcode_bits = 4;
    int selector_sites = 4;
    int sensor_sites = 1;
    int ports = 4;
    int blocks_per_port = 4;
    int pathway_blocks = 4;
    double spare_fraction = 0.10;
    int cmos_per_zone = 0;
    int edge_pathway_width = 1;
};

struct Ledger {
    int state = 0;
    int commands = 0;
    int ports = 0;
    int pathways = 0;
    int spares = 0;
    int fzc = 0;
    int edge_pathways = 0;
    int usable = 0;
};

static void require(bool ok, const std::string& what) {
    if (!ok) throw std::runtime_error("ASSERTION FAILED: " + what);
}

static bool complete(const Budget& b) {
    return b.state_groups > 0 && b.bits_per_group > 0 && b.rails == 2 && b.replicas >= 3 &&
           b.commands > 0 && b.opcode_bits > 0 && b.selector_sites >= 4 && b.sensor_sites > 0 &&
           b.ports >= 4 && b.blocks_per_port > 0 && b.pathway_blocks > 0 && b.spare_fraction > 0.0 &&
           b.cmos_per_zone == 0 && b.edge_pathway_width > 0;
}

static Ledger account(const Budget& b) {
    Ledger out;
    out.state = b.state_groups * b.bits_per_group * b.rails * b.replicas;
    out.commands = b.commands * (b.opcode_bits * b.rails + b.selector_sites + b.sensor_sites);
    out.ports = b.ports * b.blocks_per_port;
    out.pathways = b.pathway_blocks;
    const int before_spares = out.state + out.commands + out.ports + out.pathways;
    out.spares = static_cast<int>(std::ceil(b.spare_fraction * before_spares));
    out.fzc = before_spares + out.spares;
    // Each zone edge is shared by two zones, so count half of the perimeter ring.
    out.edge_pathways = 2 * kZoneSide * b.edge_pathway_width;
    out.usable = kZoneBlocks - out.fzc - out.edge_pathways;
    return out;
}

static void missing_category_fails() {
    std::cout << "\n[SCENARIO 1] a budget with no pathways or per-Zone CMOS is invalid\n";
    Budget missing = {};
    missing.pathway_blocks = 0;
    Budget cmos = {};
    cmos.cmos_per_zone = 1;
    require(!complete(missing), "pathways are a required budget term");
    require(!complete(cmos), "V3 must not replicate CMOS once per Zone");
    std::cout << "  omitted pathways rejected; per-Zone CMOS rejected\n";
}

static void provisional_512_fails_default() {
    std::cout << "\n[SCENARIO 2] declared default allocation exceeds the provisional 512\n";
    const Budget b;
    require(complete(b), "default budget must include every required category");
    const Ledger n = account(b);
    require(n.fzc > kProvisionalFzc, "default declared allocation must show that 512 is not sufficient");
    require(n.usable > 0 && n.usable < kZoneBlocks, "usable Blocks must be positive and below the raw Zone size");
    // PR6/S4: this module keeps its OWN copy of the FZC ledger (the Budget
    // struct) and gated it only against its own local kProvisionalFzc, so
    // editing the canonical allocation left it printing a stale 535 and PASSing.
    //
    // The first attempt at this gate compared n.fzc to zone_fzc_blocks() and
    // STILL survived a perturbation, because zone_fzc_blocks() is itself the
    // literal `return 535` -- a second independent definition of the same
    // number as fzc_ledger_sum() (336+130+16+4+49), agreeing by authorship
    // rather than by construction. All three paths are therefore tied together:
    // this module's ledger, the computed sum, and the documented literal.
    require(static_cast<double>(n.fzc) == fea::fzc_ledger_sum(),
            "the floorplan's local ledger must equal the computed canonical sum, "
            "or this module is accounting for a different Zone than the suite");
    require(static_cast<double>(fea::zone_fzc_blocks()) == fea::fzc_ledger_sum(),
            "zone_fzc_blocks()'s documented literal must equal the computed sum of "
            "its own parts, or the suite carries two canonical answers for one number");
    require(static_cast<double>(n.state) == fea::fzc_ledger_bits(),
            "the local state line must equal the canonical state ledger");
    std::cout << "  state=" << n.state << " commands=" << n.commands << " ports=" << n.ports
              << " pathways=" << n.pathways << " spares=" << n.spares << " fzc=" << n.fzc << "\n";
    std::cout << "  edge pathways=" << n.edge_pathways << " usable=" << n.usable
              << " of " << kZoneBlocks << "\n";
}

static void smaller_state_fits_but_is_not_layout() {
    std::cout << "\n[SCENARIO 3] fewer declared state bits can fit, which shows the target is conditional\n";
    Budget small;
    small.bits_per_group = 4;
    const Ledger n = account(small);
    require(complete(small), "reduced state budget is still complete");
    require(n.fzc < kProvisionalFzc, "4 bits per state group fits inside 512 under these other declarations");
    require(n.fzc != account(Budget{}).fzc, "the 512 result must change when the declared state width changes");
    std::cout << "  4 bits/group fzc=" << n.fzc << " fits in 512; this is not a drawn layout\n";
}

static void edge_width_dominates_fzc() {
    std::cout << "\n[SCENARIO 4] shared edge pathways can exceed the FZC itself\n";
    Budget wide;
    wide.edge_pathway_width = 4;
    const Ledger base = account(Budget{});
    const Ledger n = account(wide);
    require(n.edge_pathways > n.fzc, "a 4-Block shared edge ring must cost more Blocks than the default FZC");
    require(n.usable < base.usable, "wider pathways must reduce usable Blocks");
    std::cout << "  width=4 edge pathways=" << n.edge_pathways << " fzc=" << n.fzc << " usable=" << n.usable << "\n";
}

} // namespace floorplan

int main() {
    using namespace floorplan;
    try {
        std::cout << "FEA V3 FZC block-budget ledger\n";
        std::cout << "Counts are declared Block allocations, not a physical floorplan or density claim.\n";
        missing_category_fails();
        provisional_512_fails_default();
        smaller_state_fits_but_is_not_layout();
        edge_width_dominates_fzc();
        std::cout << "\nPASS: block-budget gates held.\n";
        std::cout << "NEXT EVIDENCE GATE: replace declared widths with a drawn layout and fabrication pitch.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
