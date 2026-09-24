// =============================================================================
// FEA_layout_v3.cpp -- M15 fabric layout parameterised by die area
//
// Everything so far assumed a fixed 3 cm^2 die. This module takes die area as
// an input and derives the layout that follows from it: Zone grid dimensions,
// FZC footprint, inter-zone routing, boundary ring, clock/bias grid, and the
// worst-case Slingshot hop count. It also checks whether V2's stated "2x
// routing overhead" is self-consistent in linear terms.
//
// Zone geometry comes from fea_params: 256 x 256 Blocks at 1.15 nm pitch.
// No term here is a drawn layout. Widths are declared.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace layout {

using fea::params;
using fea::require;

// Declared design fractions. Not measured, not from a layout.
static double routing_area_overhead() { return 2.0; }   // V2 states 2x routing overhead
static int fzc_blocks_per_zone() { return fea::zone_fzc_blocks(); }
static int clock_grain_zones() { return 256; }          // declared: one bias feed per 256 Zones

// Boundary ring WIDTH is a declared PHY/pad ring thickness, not a back-solve.
// A previous revision used 0.02166 cm, chosen so a 3 cm^2 die landed on exactly
// 5.000%, which made M15 agree with M2 because M15 had been fitted to M2. The
// The ring width is one declared value in fea_params, shared with M2 and the
// power budget, so all three agree because they use one number rather than by
// fitting one declaration to another.
static double boundary_ring_width_cm() { return fea::boundary_ring_width_cm(); }

// V2's stated Block count, line 131, lives in fea_params so this module and
// M15's tiling cross-check use the same independent input.
static double cells_from_v2_cm2() { return fea::v2_stated_blocks() * fea::raw_block_area_cm2(); }

struct Layout {
    double die_cm2 = 0.0;
    double die_side_cm = 0.0;
    double zone_side_cm = 0.0;
    double zone_area_cm2 = 0.0;
    double zones = 0.0;
    double grid_side = 0.0;
    double active_side_cm = 0.0;
    double routing_slack = 0.0;
    double fzc_blocks = 0.0;
    double fzc_area_cm2 = 0.0;
    double boundary_area_cm2 = 0.0;
    double boundary_fraction = 0.0;
    double clock_feeds = 0.0;
    double worst_case_hops = 0.0;
    double usable_cm2 = 0.0;
};

static Layout build(double die_cm2) {
    Layout L;
    L.die_cm2 = die_cm2;
    L.die_side_cm = std::sqrt(die_cm2);

    const auto& a = params().arch;
    L.zone_side_cm = a.zone_side_blocks * a.block_pitch_nm * 1e-7; // nm -> cm, DATA array side
    // Phase 1: the Zone carries the 256 x 256 data array (65536 = 1024 Words)
    // PLUS the 535-Block FZC strip, for 66071 Blocks total. This keeps V2's
    // "a Zone addresses 1,024 Words" true instead of carving FZC out of the data.
    L.zone_area_cm2 = fea::zone_total_blocks() * fea::raw_block_area_cm2();
    const double zone_pitch_cm = std::sqrt(L.zone_area_cm2);

    // Zones that fit once the declared routing overhead is paid.
    L.zones = die_cm2 / (L.zone_area_cm2 * routing_area_overhead());
    L.grid_side = std::sqrt(L.zones);
    L.active_side_cm = L.grid_side * zone_pitch_cm;
    L.routing_slack = L.die_side_cm / L.active_side_cm;

    L.fzc_blocks = fzc_blocks_per_zone() * L.zones;
    L.fzc_area_cm2 = L.fzc_blocks * fea::raw_block_area_cm2();

    // Boundary ring of declared fixed width on a square die. Fixed width means
    // the FRACTION amortises as 1/sqrt(area), which is the point of scenario 2.
    const double w = boundary_ring_width_cm();
    L.boundary_area_cm2 = 4.0 * w * L.die_side_cm - 4.0 * w * w;
    L.boundary_fraction = L.boundary_area_cm2 / die_cm2;

    // One bias/clock feed per declared grain of Zones, laid out on the die.
    L.clock_feeds = std::ceil(L.zones / clock_grain_zones());

    // Worst-case XY route across the Zone grid: 2 hops per grid step.
    L.worst_case_hops = 2.0 * (L.grid_side - 1.0);

    // Usable area after routing overhead and the boundary ring.
    L.usable_cm2 = die_cm2 - (L.zones * L.zone_area_cm2) - L.boundary_area_cm2;
    return L;
}

static void print_row(const char* tag, const Layout& L) {
    std::cout << std::fixed;
    std::cout << "  " << std::left << std::setw(9) << tag << std::right
              << std::setprecision(4) << std::setw(9) << L.die_side_cm
              << std::setprecision(0) << std::setw(13) << L.zones
              << std::setprecision(0) << std::setw(11) << L.grid_side
              << std::setprecision(4) << std::setw(9) << L.active_side_cm
              << std::setprecision(4) << std::setw(9) << L.routing_slack
              << std::setprecision(4) << std::setw(10) << L.fzc_area_cm2
              << std::setprecision(3) << std::setw(9) << (L.boundary_fraction * 100.0)
              << "%  " << std::setprecision(0) << std::setw(9) << L.worst_case_hops << "\n";
}

static void scenario_zone_geometry() {
    std::cout << "\n[SCENARIO 1] Zone geometry follows from Block pitch, not from chip size\n";
    const auto& a = params().arch;
    const double zone_side_nm = a.zone_side_blocks * a.block_pitch_nm;
    const double zone_area_cm2 = (zone_side_nm * 1e-7) * (zone_side_nm * 1e-7);
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Block pitch                 : " << a.block_pitch_nm << " nm\n";
    std::cout << "  Blocks per Zone side        : " << a.zone_side_blocks << "\n";
    std::cout << "  Zone side                   : " << std::setprecision(1) << zone_side_nm
              << " nm (" << std::setprecision(5) << (zone_side_nm * 1e-4) << " um)\n";
    std::cout << "  Zone area                   : " << std::scientific << std::setprecision(3)
              << zone_area_cm2 << " cm^2" << std::fixed << "\n";
    std::cout << "  Blocks per Zone             : " << (a.zone_side_blocks * a.zone_side_blocks)
              << "\n";
    require(zone_side_nm > 200.0 && zone_side_nm < 400.0,
            "a 256-Block Zone at 1.15 nm pitch must be a few hundred nm across");
    require(zone_area_cm2 > 0.0, "Zone area must be positive");
    std::cout << "  a Zone is a ~300 nm object. Chip size does not change Zone geometry, only\n";
    std::cout << "  how many of them are tiled and how far apart the far ones sit.\n";
}

static void scenario_sweep_die_area() {
    std::cout << "\n[SCENARIO 2] layout as a function of die area\n";
    const std::vector<double> areas{1.0, 3.0, 10.0};
    std::cout << "  " << std::left << std::setw(9) << "die cm^2" << std::right
              << std::setw(9) << "side cm" << std::setw(13) << "zones"
              << std::setw(11) << "grid" << std::setw(9) << "array cm"
              << std::setw(9) << "slack" << std::setw(10) << "FZC cm^2"
              << std::setw(9) << "bnd" << std::setw(9) << "max hops" << "\n";
    std::vector<Layout> layouts;
    for (const double a : areas) {
        const Layout L = build(a);
        layouts.push_back(L);
        print_row(std::to_string(static_cast<int>(a)).c_str(), L);
    }

    const Layout& one = layouts[0];
    const Layout& three = layouts[1];
    const Layout& ten = layouts[2];

    require(three.zones / one.zones > 2.9 && three.zones / one.zones < 3.1,
            "Zone count must scale linearly with die area");
    require(ten.worst_case_hops / one.worst_case_hops > 3.0,
            "worst-case hops must grow with the grid, not stay flat");

    // A previous revision asserted routing_slack == sqrt(overhead) here. That
    // was an identity: zones is DEFINED as die/(zone_area*overhead), which forces
    // slack to sqrt(overhead) for any value of overhead, so the check re-derived
    // its own input and could never fail. Replaced by a genuine cross-check:
    // does this module's Zone tiling reproduce the Block count V2 independently
    // states? Two different routes to the same cell area.
    const double data_tiled = three.zones * fea::zone_data_blocks() * fea::raw_block_area_cm2();
    const double total_tiled = three.zones * three.zone_area_cm2;
    const double cells_stated = cells_from_v2_cm2();
    const double tiling_error = std::abs(data_tiled - cells_stated) / cells_stated;
    std::cout << "\n" << std::setprecision(4);
    std::cout << "  data cells from this module's Zone tiling : " << data_tiled << " cm^2\n";
    std::cout << "  cells from V2's stated 1.13e14 Blocks    : " << cells_stated << " cm^2\n";
    std::cout << "  disagreement between the two routes       : " << std::setprecision(3)
              << (tiling_error * 100.0) << " %\n";
    std::cout << "  Zone total including the FZC strip       : " << std::setprecision(4)
              << total_tiled << " cm^2 of " << three.die_cm2 / routing_area_overhead()
              << " cm^2 budgeted\n";
    require(tiling_error < 0.01,
            "die-geometry occupancy must reproduce V2's separately stated Block count within 1%");
    // routing_slack is sqrt(routing_area_overhead) BY ALGEBRA: zones is DEFINED as
    // die/(zone_area*routing), which forces slack to sqrt(routing) for any value.
    // It is an identity, so it must not be gated. An earlier revision gated
    // |slack - sqrt(2)| < 1e-6 (dead) and a later one gated 1.0 < slack < 2.0,
    // which was the same identity inside a wider bracket (also dead).
    std::cout << "  routing slack                : " << std::setprecision(4)
              << three.routing_slack << " = sqrt(" << routing_area_overhead()
              << ") BY CONSTRUCTION, an identity, not a measurement.\n";
    std::cout << "  note: data_tiled cancels block_area (zones already divide by it), so the\n";
    std::cout << "  check above is die-geometry occupancy versus V2's stated count x block\n";
    std::cout << "  area, not two pitch-derived routes. Claiming independent pitch routes was\n";
    std::cout << "  wrong; the comparison is still a real die-versus-claim test.\n";

    std::cout << "  zones scale linearly: " << std::setprecision(1)
              << (three.zones / one.zones) << "x area -> " << (three.zones / one.zones)
              << "x zones, and " << (ten.zones / one.zones) << "x at 10 cm^2.\n";
    std::cout << "  worst-case Slingshot hops grow sub-linearly: " << std::setprecision(3)
              << (ten.worst_case_hops / one.worst_case_hops) << "x for " << (ten.die_cm2 / one.die_cm2)
              << "x the area, because hops follow the grid side, not the area.\n";
    std::cout << "  boundary ring share FALLS as the die grows: " << std::setprecision(1)
              << (one.boundary_fraction * 100.0) << "% at 1 cm^2 -> "
              << (three.boundary_fraction * 100.0) << "% at 3 cm^2 -> "
              << (ten.boundary_fraction * 100.0) << "% at 10 cm^2.\n";
    require(ten.boundary_fraction < one.boundary_fraction,
            "the perimeter-scaled boundary ring must amortise over a larger die");
    // The width is now ONE declared input in fea_params, shared with M2 and the
    // power budget, so all three agree by construction of a single number
    // instead of by fitting one to another. What remains to check is that the
    // ring never swallows the die it sits on.
    require(fea::boundary_ring_area_cm2(three.die_cm2) < 0.25 * three.die_cm2,
            "the declared boundary ring must stay under a quarter of the die");
    std::cout << "  declared ring " << std::setprecision(3)
              << (fea::boundary_ring_width_cm() * 1000.0) << " mm costs "
              << (fea::boundary_ring_fraction_at(3.0) * 100.0) << "% of V2's 3 cm^2 die but "
              << (fea::boundary_ring_fraction_at(0.5) * 100.0)
              << "% of our 0.5 cm^2 design die.\n";
    std::cout << "  the width is fixed because boundary pads do not shrink, so a SMALLER die\n";
    std::cout << "  pays a LARGER share for I-O. Shrinking the die to ease fabrication costs\n";
    std::cout << "  this, and it is budgeted rather than hidden.\n";
    std::cout << "  so a bigger chip amortises the boundary ring, which is why the ratio is\n";
    std::cout << "  perimeter over area: it falls as 1/sqrt(area).\n";
}

static void scenario_clock_and_hops() {
    std::cout << "\n[SCENARIO 3] clock feeds and worst-case route length\n";
    const Layout L = build(3.0);
    std::cout << std::fixed << std::setprecision(0);
    std::cout << "  Zones                       : " << L.zones << "\n";
    std::cout << "  Zone grid side              : " << L.grid_side << " Zones\n";
    std::cout << "  clock feeds at 1 per " << clock_grain_zones() << " Zones  : " << L.clock_feeds << "\n";
    std::cout << "  worst-case XY route         : " << L.worst_case_hops << " Slingshot hops\n";
    require(L.clock_feeds > 1.0, "a large die needs more than one clock feed");
    require(L.worst_case_hops > L.grid_side, "worst-case route must exceed one grid side");
    std::cout << "\n  a Zone at one corner talking to the opposite corner crosses about "
              << L.worst_case_hops << " hops.\n";
    std::cout << "  that number is what turns M9's local cycle into a system latency, and it\n";
    std::cout << "  scales as sqrt(die area), not linearly with Zone count.\n";
    std::cout << "  label: geometry derived from pitch, widths and grain DECLARED. Not a drawing.\n";
}

static void scenario_capacity_options() {
    std::cout << "\n[SCENARIO 4] capacity options: every structure must live inside the die\n";
    const Layout L = build(3.0);
    const double die = 3.0;
    const double cells_stated = cells_from_v2_cm2(); // V2's own Block count, independent
    const double fzc = L.fzc_area_cm2;
    const double pdn10 = 0.10 * die;
    const double pdn5 = 0.05 * die;
    const double boundary = L.boundary_area_cm2;
    const double v2_capacity_TB = 14.1;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  V2: Block = 3a x 3a = 1.15 x 1.15 nm^2, practical density 3.77e13 cm^-2\n";
    std::cout << "  after 2x routing. Cells taken from V2's stated 1.13e14 Block count:\n";
    std::cout << "    Block cells             : " << cells_stated << " cm^2\n";
    std::cout << "    routing (V2's own 2x)   : " << (cells_stated * (routing_area_overhead() - 1.0))
              << " cm^2\n";
    std::cout << "    area V2 claims outright : " << (cells_stated * routing_area_overhead())
              << " cm^2  = " << std::setprecision(1)
              << (cells_stated * routing_area_overhead() / die * 100.0)
              << "% of a " << die << " cm^2 die\n\n";

    // Binding constraint, checked FEASIBILITY rather than identity: does V2's
    // stated cell requirement still fit once support is reserved from the die?
    // total is NOT forced to equal die, because cells_stated comes from V2.
    auto required_area = [&](double support, double routing_factor) {
        return cells_stated * routing_factor + support;
    };
    auto capacity_for = [&](double support, double routing_factor) {
        const double cells = (die - support) / routing_factor;
        // Geometry-derived capacity: how many Blocks those cells hold at the
        // physical Block area, converted to bytes. NOT proportional to V2's
        // claimed 14.1 TB, which is what makes the option-A gate below a real
        // cross-check instead of v2_capacity_TB * k ~= v2_capacity_TB.
        const double blocks = cells / fea::raw_block_area_cm2();
        return std::make_pair(cells, blocks / 8.0 / 1e12);
    };

    struct Option {
        const char* name;
        double support;
        double routing_factor;
        const char* basis;
    };
    const Option opts[] = {
        {"A  V2 as published", 0.0, routing_area_overhead(), "assumes PDN, boundary and control are free"},
        {"B  FZC only", fzc, routing_area_overhead(), "reserve control, still free PDN and I/O"},
        {"C  FZC + PDN 5%", fzc + pdn5, routing_area_overhead(), "cheaper power grid, no boundary ring"},
        {"D  FZC + PDN 10% + boundary", fzc + pdn10 + boundary, routing_area_overhead(), "every declared structure paid"},
        {"E  D but routing 1.5x", fzc + pdn10 + boundary, 1.5, "needs justification for lower routing"},
    };

    std::cout << "  option                          support cm^2   routing   cells cm^2   capacity   area needed   fits?\n";
    double previous_same_routing = 1e9;
    double previous_routing = -1.0;
    double previous_required = -1.0;
    for (const Option& o : opts) {
        const auto pr = capacity_for(o.support, o.routing_factor);
        const double needed = required_area(o.support, o.routing_factor);
        const bool fits = needed <= die + 1e-9;
        std::cout << "  " << std::left << std::setw(31) << o.name << std::right
                  << std::setprecision(4) << std::setw(12) << o.support
                  << std::setprecision(1) << std::setw(9) << o.routing_factor
                  << std::setprecision(4) << std::setw(12) << pr.first
                  << std::setprecision(2) << std::setw(9) << pr.second << " TB"
                  << std::setprecision(4) << std::setw(13) << needed << "   "
                  << (fits ? "yes" : "NO") << "\n";
        // Capacity falls as more support is reserved at a fixed routing factor.
        if (std::abs(o.routing_factor - previous_routing) < 1e-12) {
            require(pr.second < previous_same_routing,
                    "at fixed routing, reserving more support must not raise capacity");
            require(needed > previous_required,
                    "reserving more support must demand more area, not less");
        }
        previous_same_routing = pr.second;
        previous_routing = o.routing_factor;
        previous_required = needed;
        std::cout << "      basis: " << o.basis << "\n";
    }

    const auto optA = capacity_for(0.0, routing_area_overhead());
    const auto optD = capacity_for(fzc + pdn10 + boundary, routing_area_overhead());
    const auto optE = capacity_for(fzc + pdn10 + boundary, 1.5);
    const double areaA = required_area(0.0, routing_area_overhead());
    const double areaD = required_area(fzc + pdn10 + boundary, routing_area_overhead());
    std::cout << "\n" << std::setprecision(4);
    require(areaA <= die,
            "V2's stated cell count with zero support reserved must still fit the die");
    require(areaD > die,
            "V2's stated cell count plus every declared support structure must NOT fit the die");
    std::cout << "  A needs " << areaA << " cm^2 of a " << die
              << " cm^2 die -> fits, with " << (die - areaA) << " cm^2 left over.\n";
    std::cout << "  D needs " << areaD << " cm^2 of a " << die
              << " cm^2 die -> DOES NOT FIT, short by " << (areaD - die)
              << " cm^2 (" << std::setprecision(1) << ((areaD / die - 1.0) * 100.0) << "%).\n\n";
    // PR5/m8: the gate that stood here compared optA against 14.1 TB, but both
    // come from the same stated cell count (1.13e14 bits / 8 = 14.125 TB), so it
    // re-derived its own definition and only a header edit could fail it. The
    // independent evidence in this scenario is the AREA pair above, which uses
    // cells against a literal die size. What replaces it tests that reserving
    // support actually costs capacity, which fails if any support term drops out.
    require(optA.second - optD.second > 1.0,
            "paying for every declared structure must cost at least 1 TB of capacity, "
            "otherwise a support term is not entering the calculation");
    require(optD.second > 11.0 && optD.second < 12.5,
            "paying for every declared structure must land near 12 TB");
    require(optE.second > v2_capacity_TB,
            "relaxing routing below V2's own 2x would exceed the published claim");

    std::cout << "  D is the physically correct reading: the die is a fixed " << die
              << " cm^2, so cells,\n";
    std::cout << "  routing, FZC, PDN and the boundary ring must all come out of it. That\n";
    std::cout << "  takes capacity from " << optA.second << " TB to " << optD.second << " TB.\n";
    std::cout << "  A is not an option in the physical sense. It is what you get only if PDN,\n";
    std::cout << "  boundary I/O and control genuinely cost zero area, which they do not.\n";
    std::cout << "  E is the only route back above 14.1 TB, and it requires justifying a routing\n";
    std::cout << "  factor below the 2x V2 itself states.\n";
    std::cout << "  FZC alone needs " << std::setprecision(4) << fzc
              << " cm^2 against V2's residual of " << std::setprecision(4) << (die - areaA)
              << " cm^2, so it barely tips option B over (" << std::setprecision(1)
              << ((fzc + cells_stated * routing_area_overhead() - die) * 100.0 / die)
              << "% of die). But FZC is only " << std::setprecision(1)
              << (fzc / (fzc + pdn10 + boundary) * 100.0)
              << "% of the support term: PDN and the boundary ring are "
              << ((pdn10 + boundary) / (fzc + pdn10 + boundary) * 100.0)
              << "% of it. FZC is the smallest term, and it is not what drives the\n";
    std::cout << "  shortfall. V2's published numbers simply leave almost no headroom.\n";
    std::cout << "  label: cell area from V2's stated Block count, feasibility checked against the\n";
    std::cout << "  die. The three identities that used to sit here were removed: slack and the\n";
    std::cout << "  budget total were true by construction, not by measurement.\n";
}


// =============================================================================
// SCENARIO 8: capacity options at OUR 0.5 cm^2 design die
//
// PR8/C04: fig7(c) plotted five design-die capacities (2.363 / 2.353 / 2.235 /
// 1.832 / 2.443 TB) that NO captured output produced -- scenario 4 evaluates
// options only at the 3 cm^2 audit basis (14.18 / 14.12 / 13.41 / 11.99 /
// 15.99 TB). The same option structure is therefore solved at the design die
// here, from the canonical Block area, so the figure has a producing line and
// a single source of truth.
static void scenario_design_die_options() {
    std::cout << "\n[SCENARIO 8] capacity options at the 0.5 cm^2 design die\n";

    const Layout L = build(params().arch.die_area_cm2);
    const double die = params().arch.die_area_cm2;
    const double fzc = L.fzc_area_cm2;
    const double pdn5 = 0.05 * die;
    const double pdn10 = 0.10 * die;
    const double boundary = L.boundary_area_cm2;

    auto capacity_TB = [&](double support, double routing_factor) {
        const double cells = (die - support) / routing_factor;
        return (cells / fea::raw_block_area_cm2()) / 8.0 / 1e12;
    };

    struct Option {
        const char* code;
        const char* name;
        double support;
        double routing;
    };
    const Option opts[] = {
        {"A", "no support", 0.0, routing_area_overhead()},
        {"B", "+FZC", fzc, routing_area_overhead()},
        {"C", "+FZC+PDN5", fzc + pdn5, routing_area_overhead()},
        {"D", "+FZC+PDN10+boundary", fzc + pdn10 + boundary, routing_area_overhead()},
        {"E", "D at routing 1.5x", fzc + pdn10 + boundary, 1.5},
    };

    // fixed(4) would print the Block area (~1.3e-14 cm^2) as "0.0000", which
    // reads like a bug, so that one line gets scientific notation.
    std::cout << "  die " << die << " cm^2, Block area " << std::scientific
              << std::setprecision(4) << fea::raw_block_area_cm2() << " cm^2"
              << std::fixed << std::setprecision(4) << ", routing "
              << routing_area_overhead() << "x\n";
    std::cout << "  support: FZC " << fzc << ", PDN5 " << pdn5 << ", PDN10 " << pdn10
              << ", boundary ring " << boundary << " cm^2\n\n";
    std::cout << "  code  option                     support cm^2  routing   capacity TB\n";

    double first = -1.0, last_same = 1e9;
    for (const Option& o : opts) {
        const double tb = capacity_TB(o.support, o.routing);
        std::cout << "  " << std::left << std::setw(6) << o.code << std::setw(27) << o.name
                  << std::right << std::setw(11) << o.support << std::setw(9) << o.routing
                  << std::setw(14) << tb << "\n";
        if (first < 0.0) first = tb;
        if (o.routing == routing_area_overhead()) {
            require(tb < last_same + 1e-9,
                    "reserving more support must never increase capacity, or a term "
                    "is not entering the subtraction");
            last_same = tb;
        }
    }

    const double capA = capacity_TB(0.0, routing_area_overhead());
    const double capD = capacity_TB(fzc + pdn10 + boundary, routing_area_overhead());
    const double capE = capacity_TB(fzc + pdn10 + boundary, 1.5);
    std::cout << std::setprecision(6);
    require(capA > 0.0 && capD > 0.0 && capE > 0.0,
            "every option must yield a positive capacity");
    require(capD < capA,
            "paying for every declared structure must cost capacity at the design die");
    require(capE > capA,
            "relaxing routing below the stated 2x must exceed the no-support figure, "
            "which is what marks it as needing justification");
    require(capD > 1.0,
            "the fully-paid design die must still hold more than a terabyte, or the "
            "0.5 cm^2 target is not viable as stated");
    std::cout << "  A " << capA << " TB, D " << capD << " TB, E " << capE << " TB\n";
    std::cout << "  label: ARITHMETIC from the canonical Block area at the design die.\n";
}

} // namespace layout

int main() {
    using namespace layout;
    try {
        std::cout << "FEA V3 M15 fabric layout as a function of die area\n";
        std::cout << "Die area is an input. Geometry comes from Block pitch. Widths are declared.\n";
        scenario_zone_geometry();
        scenario_sweep_die_area();
        scenario_clock_and_hops();
        scenario_capacity_options();
        scenario_design_die_options();
        std::cout << "\nPASS: Zone grid, FZC, boundary ring, clock feeds and hop count all scale correctly.\n";
        std::cout << "LABEL: derived geometry, declared widths and clock grain, no drawn layout.\n";
        std::cout << "NEXT EVIDENCE GATE: an actual placement of Zones, FZC, ports and boundary PHYs.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
