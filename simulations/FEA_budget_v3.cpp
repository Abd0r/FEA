// =============================================================================
// FEA_budget_v3.cpp -- M1 reconciled control-plane power and area budget
//
// Computes every V2 headline power term from fea_params.h instead of restating
// it. Reviewers found that V2 wrote 3.3 W where the arithmetic gives 3300 W,
// and 0.14 uW where 15 fJ x 9.19 GHz gives 137.85 uW per Zone. This module
// reproduces the reviewer arithmetic, then reports the corrected totals.
//
// It does not claim a tape-out. It claims only that the arithmetic is checked.
// =============================================================================

#include "fea_params.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace budget {

using fea::check_unit;
using fea::params;
using fea::require;
using fea::zone_count_stated;

struct Term {
    std::string name;
    double watts = 0.0;
    double area_cm2 = 0.0;
    std::string note;
};

// Zones on the die, derived from V2's own stated count. Kept as a parameter
// rather than recomputed from pitch so the comparison to the reviewer's own
// V2 states ~1.7e9 Zones (line 130). One definition, from fea_params.
static double zone_count() { return zone_count_stated(); }

static double f_sys() {
    const auto& t = params().timing;
    const double cycle_s = (t.t_arm_ps + t.t_fire_ps + t.t_confirm_ps) * fea::kPS;
    return 1.0 / cycle_s;
}

static std::vector<Term> build_terms() {
    const auto& c = params().control;
    const auto& a = params().arch;
    std::vector<Term> terms;

    const double zones = zone_count();
    const double f = f_sys();

    // Decoder: dynamic energy per event times clock times number of Zones.
    const double decoder_per_zone_W = c.decoder_event_J * f;
    check_unit("decoder power per Zone", decoder_per_zone_W, "W", 1.3785e-4);
    terms.push_back({"decoder", decoder_per_zone_W * zones, 0.0,
                     "15 fJ x 9.19 GHz x 1.7e9 Zones"});

    // PLL: count shrinks by sharing group K, then power per instance.
    const double pll_count = zones / c.pll_group_K;
    check_unit("PLL count", pll_count, "count", 6.64e6, 0.10);
    const double pll_W = pll_count * c.pll_power_W;
    check_unit("aggregate PLL power", pll_W, "W", 3320.0, 0.10);
    terms.push_back({"pll", pll_W, 0.0, "1.7e9 / 256 PLLs x 0.5 mW"});

    // Sequencer and charge sensing: declared per-Zone figure.
    const double sense_W = c.sense_zone_W * zones;
    check_unit("sequencer and sensing power", sense_W, "W", 170.0, 0.10);
    terms.push_back({"sequencer-sensing", sense_W, 0.0, "0.1 uW/Zone x 1.7e9 Zones"});

    // Data plane: power density times data-plane area.
    const double data_W = c.data_plane_mW_per_cm2 * 1e-3 * c.data_plane_area_cm2;
    // Expected order: 26.47 mW/cm^2 x 0.5 cm^2 = 0.013235 W. The literal was
    // 0.0794, the 3 cm^2 figure, and went stale the moment the design die
    // changed. This gate guards arithmetic slip, not the power density itself.
    check_unit("data-plane power", data_W, "W", 0.013235, 0.10);
    terms.push_back({"data-plane", data_W, 0.0,
                     "26.47 mW/cm2 x " + std::to_string(c.data_plane_area_cm2) + " cm2"});

    // Decoder area: per-Zone area times Zone count.
    const double decoder_area_cm2 = c.decoder_area_um2 * fea::kCM2_PER_UM2 * zones;
    check_unit("decoder area", decoder_area_cm2, "cm^2", 204.0, 0.10);
    terms[0].area_cm2 = decoder_area_cm2;
    terms[0].note += "; area 1.7e9 x 12 um^2";

    (void)a;
    return terms;
}

static void scenario_unit_checks() {
    std::cout << "\n[SCENARIO 1] every V2 control-plane term is recomputed and unit-checked\n";
    const auto terms = build_terms();
    double total_W = 0.0;
    double total_area = 0.0;
    std::cout << std::fixed << std::setprecision(3);
    for (const Term& t : terms) {
        std::cout << "  " << std::left << std::setw(20) << t.name << std::right
                  << std::setw(14) << t.watts << " W";
        if (t.area_cm2 > 0.0) std::cout << "   area " << std::setprecision(2) << t.area_cm2 << " cm^2";
        std::cout << std::setprecision(3) << "   (" << t.note << ")\n";
        total_W += t.watts;
        total_area += t.area_cm2;
    }
    std::cout << "  " << std::left << std::setw(20) << "TOTAL" << std::right
              << std::setw(14) << total_W << " W";
    if (total_area > 0.0) std::cout << "   area " << std::setprecision(2) << total_area << " cm^2";
    std::cout << "\n";
    require(total_W > 1000.0, "reconciled per-Zone control-plane power must expose the kilowatt-scale V2 error");
    require(total_W > 100.0 * 3.8, "corrected total must be far above the V2 3.8 W claim");
    std::cout << "  V2 claimed 3.8 W. The per-Zone CMOS control plane this replaces costs "
              << std::setprecision(1) << total_W << " W.\n";
    std::cout << "  reading: this is the cost of the approach FZC eliminates, not a V3 budget line.\n";
}

static void scenario_decoder_area_exceeds_die() {
    std::cout << "\n[SCENARIO 2] decoder area alone exceeds the die\n";
    const auto terms = build_terms();
    const double decoder_area = terms[0].area_cm2;
    // V2's claim: 12 um^2 per Zone x V2's own 1.7e9 Zones = 204 cm^2. That was
    // written for V2's 3 cm^2 die, so it is compared against V2's reference die.
    // Using OUR 0.5 cm^2 here would divide V2's zone count by our area and
    // produce a ratio true of neither die.
    const double die = params().arch.v2_reference_die_cm2;
    require(decoder_area > die, "decoder area must exceed V2's 3 cm^2 die to expose the V2 error");
    require(decoder_area > 50.0 * die, "decoder overcommit must be more than 50x, matching the reviewer figure");
    std::cout << "  decoder area " << std::fixed << std::setprecision(1) << decoder_area
              << " cm^2 vs die " << die << " cm^2 -> overcommit "
              << std::setprecision(0) << (decoder_area / die) << "x\n";
    std::cout << "  a per-Zone CMOS decoder cannot fit. This is why V3 moves control into FZC.\n";
}

static void scenario_v2_literal_values_fail_checks() {
    std::cout << "\n[SCENARIO 3] the literal V2 values fail their own unit checks\n";
    const auto& c = params().control;
    const double f = f_sys();
    const double correct_decoder = c.decoder_event_J * f;
    const double v2_decoder = 0.14e-6;
    const double correct_pll = (zone_count() / c.pll_group_K) * c.pll_power_W;
    const double v2_pll = 3.3;

    bool decoder_ratio_bad = std::abs(correct_decoder / v2_decoder - 1.0) > 0.05;
    bool pll_ratio_bad = std::abs(correct_pll / v2_pll - 1.0) > 0.05;
    require(decoder_ratio_bad, "V2 0.14 uW decoder figure must not survive a unit check");
    require(pll_ratio_bad, "V2 3.3 W PLL figure must not survive a unit check");
    std::cout << "  V2 decoder 0.14 uW vs correct " << std::fixed << std::setprecision(2)
              << (correct_decoder * 1e6) << " uW -> factor " << std::setprecision(0)
              << (correct_decoder / v2_decoder) << "\n";
    std::cout << "  V2 PLL 3.3 W vs correct " << std::setprecision(0) << correct_pll
              << " W -> factor " << (correct_pll / v2_pll) << "\n";
}

static void scenario_reviewer_hand_calculations() {
    std::cout << "\n[SCENARIO 4] reproduce the reviewers' own hand calculations from V2 text\n";
    const auto& c = params().control;
    const double zones = zone_count();

    // Reviewer 7 and Reviewer 6: take V2's stated per-Zone figure as written,
    // then multiply by the Zone count. This is how they got 238 W and 170 W.
    const double v2_decoder_per_zone_W = 0.14e-6;
    const double reviewer_decoder_W = v2_decoder_per_zone_W * zones;
    check_unit("reviewer decoder aggregate", reviewer_decoder_W, "W", 238.0, 0.10);

    const double reviewer_sense_W = c.sense_zone_W * zones;
    check_unit("reviewer sensing aggregate", reviewer_sense_W, "W", 170.0, 0.10);

    const double reviewer_pll_W = (zones / c.pll_group_K) * c.pll_power_W;
    check_unit("reviewer PLL aggregate", reviewer_pll_W, "W", 3300.0, 0.05);

    std::cout << "  decoder (V2 stated 0.14 uW/Zone x 1.7e9): " << std::fixed
              << std::setprecision(1) << reviewer_decoder_W << " W   (reviewer said ~238 W)\n";
    std::cout << "  sensing/sequencer (0.1 uW/Zone x 1.7e9):  "
              << std::setprecision(1) << reviewer_sense_W << " W   (reviewer said ~170 W)\n";
    std::cout << "  PLL (1.7e9/256 x 0.5 mW):                  "
              << std::setprecision(1) << reviewer_pll_W << " W   (reviewer said ~3300 W)\n";

    // There are two self-consistent readings of V2's decoder text and both fail.
    const double correct_per_zone = c.decoder_event_J * f_sys();
    const double as_arithmetic_W = correct_per_zone * zones;
    require(as_arithmetic_W > reviewer_decoder_W,
            "correcting 0.14 uW to 137.85 uW makes the decoder aggregate worse, not better");
    std::cout << "  decoder if 15 fJ x 9.19 GHz is taken literally: "
              << std::setprecision(1) << (correct_per_zone * 1e6) << " uW/Zone, "
              << as_arithmetic_W << " W aggregate\n";
    std::cout << "  either reading is orders of magnitude above V2's stated 0.4 W decoder+sequencer.\n";
    std::cout << "  V2 total was 3.8 W. Lowest defensible control-plane floor is "
              << std::setprecision(0)
              << (reviewer_decoder_W + reviewer_sense_W + reviewer_pll_W) << " W.\n";
    require(reviewer_decoder_W + reviewer_sense_W + reviewer_pll_W > 3000.0,
            "even the reviewers' most generous reading must exceed the V2 3.8 W claim");
}

static void scenario_fzc_boundary_alternative() {
    std::cout << "\n[SCENARIO 5] V3 boundary-CMOS alternative replaces per-Zone control\n";
    // V3 keeps CMOS only at the chip boundary, not once per Zone.
    const double zones = zone_count();
    const auto& a = params().arch;
    const double edge_length_m = 2.0 * std::sqrt(a.die_area_cm2 * 1e-4);
    require(edge_length_m > 0.0 && zones > 0.0, "geometry inputs must be positive");

    // Per-Zone decoder area eliminated entirely.
    const double per_zone_decoder_area =
        params().control.decoder_area_um2 * fea::kCM2_PER_UM2 * zones;
    require(per_zone_decoder_area > 0.0, "per-Zone decoder area must be a positive quantity to eliminate");

    std::cout << "  eliminated per-Zone decoder area: " << std::fixed << std::setprecision(1)
              << per_zone_decoder_area << " cm^2\n";
    std::cout << "  eliminated per-Zone decoder power: " << std::setprecision(1)
              << (params().control.decoder_event_J * f_sys() * zones) << " W\n";
    std::cout << "  eliminated PLL instances: " << std::setprecision(2)
              << (zones / params().control.pll_group_K) << " (replaced by boundary bias/clock)\n";
    std::cout << "  note: boundary CMOS still costs area and power. M2 and M1 must size it.\n";
}

// The lowest defensible per-Zone CMOS control-plane total, using V2's own
// stated per-zone figures rather than its literal arithmetic. Reviewers 6 and 7
// both land in this vicinity. Used only for the contrast line in SCENARIO 7.
static double reviewer_floor_W() {
    const auto& c = params().control;
    const double zones = zone_count();
    const double decoder = 0.14e-6 * zones;
    const double pll = (zones / c.pll_group_K) * c.pll_power_W;
    const double sense = c.sense_zone_W * zones;
    return decoder + pll + sense;
}

static void scenario_total_power() {
    std::cout << "\n[SCENARIO 7] whole-chip total power: floor, open terms, and no false precision\n";
    const auto& c = params().control;
    const auto& g = params().gaps;

    // per_zone_cmos is 0 STRUCTURALLY: FZC replaces the per-Zone CMOS controller,
    // so this term is excluded from the V3 floor by design. It stays in the sum as
    // an explicit zero so the exclusion is visible. The gate that used to sit here,
    // `require(per_zone_cmos == 0.0)`, asserted that constant and could never fail.
    const double per_zone_cmos = 0.0; // eliminated by construction
    const double data_plane = c.data_plane_mW_per_cm2 * 1e-3 * c.data_plane_area_cm2;
    const double restoration = fea::restoration_power_W();
    const double refresh = (fea::payload_bits() / fea::refresh_interval_s()) *
                           g.refresh_energy_per_bit_J;

    // V2 omitted the last four. Each now carries a DECLARED value with a stated
    // basis (fea_params, and swept in scenario 8), but none is SOURCED, so none
    // counts toward the floor. Floor and declared total are reported apart.
    const double declared_boundary = fea::boundary_ring_W();
    const double declared_clock = fea::clock_distribution_W();
    const double declared_io = fea::external_io_W();
    const double sized_sum = per_zone_cmos + data_plane + restoration + refresh;
    const double declared_pdn =
        (sized_sum + declared_boundary + declared_clock + declared_io) * fea::pdn_loss_fraction();
    const double declared_four = declared_boundary + declared_io + declared_clock + declared_pdn;

    struct Open { const char* name; double value; bool sized; };
    const Open terms[] = {
        {"per-Zone CMOS (FZC eliminates)", per_zone_cmos, true},
        {"data plane (V2 figure)", data_plane, true},
        {"restoration endpoints (M6)", restoration, true},
        {"refresh (M13, declared energy)", refresh, true},
        {"boundary CMOS ring (declared)", declared_boundary, false},
        {"external I/O (declared)", declared_io, false},
        {"clock and bias distribution (declared)", declared_clock, false},
        {"power-delivery losses (declared)", declared_pdn, false},
    };

    double known = 0.0;
    int open_count = 0;
    int total_count = 0;
    std::cout << std::scientific << std::setprecision(4);
    for (const Open& t : terms) {
        ++total_count;
        if (t.sized) known += t.value;
        else ++open_count;
        std::cout << "  " << std::left << std::setw(34) << t.name << std::right
                  << std::setw(12) << t.value << " W   "
                  << (t.sized ? "known" : "DECLARED, unsourced") << "\n";
    }
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  ----------------------------------------------\n";
    std::cout << "  FLOOR (" << (total_count - open_count) << " of " << total_count
              << " terms)          " << std::setw(12) << known << " W\n";
    std::cout << "  unsourced (declared, not measured)    " << std::setw(4) << open_count << "\n";
    std::cout << "  declared sum of those four            " << std::setw(12) << declared_four << " W\n";
    std::cout << "  DECLARED TOTAL (all 8 terms)         " << std::setw(12) << (known + declared_four)
              << " W\n";
    std::cout << "  upper bound                           see SCENARIO 8 for the swept interval\n";

    require(known > 0.0, "the floor must be positive");
    require(open_count > 0, "the budget must keep reporting open terms rather than hiding them");
    // Two counting paths must agree: this term list and fea_params' own flags.
    require(open_count == fea::open_power_term_count(),
            "M1 and M11 must agree on how many power terms are still unsourced");
    require(declared_four > 0.0,
            "declaring the four terms must produce a positive sum, not the old four zeros");

    // Scale comparison, not an invented absolute: what if each open term were
    // as large as the data plane?
    const double hypothetical = known + open_count * data_plane;
    std::cout << "  if each open term equalled the data plane: " << hypothetical << " W\n";
    std::cout << "  STATUS: NOT VALIDATED. " << open_count << " of " << total_count
              << " terms declared but unsourced.\n";
    std::cout << "  what is settled: the total does not contain the "
              << std::setprecision(0) << (reviewer_floor_W()) << " W per-Zone control plane.\n";
    std::cout << "  label: floor DERIVED from stated terms; four terms DECLARED with a basis but\n";
    std::cout << "  UNSOURCED, so the total stays unvalidated and the M11 ratio stays withheld.\n";
}

static void scenario_v3_control_budget() {
    std::cout << "\n[SCENARIO 6] the V3 control budget is FZC plus boundary CMOS only\n";
    const auto& c = params().control;
    const double zones = zone_count();

    // In V3 the per-Zone CMOS terms are zero by construction: FZC replaces them.
    const double eliminated_decoder = c.decoder_event_J * f_sys() * zones;
    const double eliminated_pll = (zones / c.pll_group_K) * c.pll_power_W;
    const double eliminated_sense = c.sense_zone_W * zones;
    const double eliminated_total = eliminated_decoder + eliminated_pll + eliminated_sense;

    // What V3 actually pays: native FZC control sits in the data plane, plus a
    // declared boundary ring for boot, bias/clock, diagnostics, and external I-O.
    const double data_plane = c.data_plane_mW_per_cm2 * 1e-3 * c.data_plane_area_cm2;
    const double boundary_fraction = fea::boundary_ring_fraction_at(params().arch.die_area_cm2);
    const double boundary_power = 0.0; // sizing is an open gate, declared not estimated

    std::cout << "  eliminated per-Zone terms: " << std::fixed << std::setprecision(1)
              << eliminated_decoder << " W decoder + " << eliminated_pll << " W PLL + "
              << eliminated_sense << " W sensing\n";
    std::cout << "  eliminated PLL instances: " << std::setprecision(0)
              << (zones / c.pll_group_K) << "\n";
    std::cout << "  eliminated decoder area:  " << std::setprecision(1)
              << (c.decoder_area_um2 * fea::kCM2_PER_UM2 * zones) << " cm^2\n";
    std::cout << "  V3 control lines that remain:\n";
    std::cout << "    FZC native control     : in data plane, " << std::setprecision(3)
              << data_plane << " W (FZC Block count from fzc-floorplan ledger)\n";
    std::cout << "    boundary CMOS ring     : " << std::setprecision(0)
              << (boundary_fraction * 100) << "% of die declared; power "
              << std::setprecision(1) << boundary_power << " W is an OPEN gate\n";
    std::cout << "    per-Zone CMOS         : 0 by construction\n";

    require(eliminated_total > 3000.0,
            "the eliminated per-Zone cost must be the dominant term that FZC removes");
    std::cout << "  V3 does not inherit the kilowatt-scale per-Zone cost. It inherits only\n";
    std::cout << "  a boundary ring and a native FZC footprint, both of which still need sizing.\n";
    std::cout << "  label: eliminated terms derived, boundary power open, FZC count from M2 ledger.\n";
}

// Reviewer 4 points 4, 5 and 6 asked for control-plane overhead, I/O energy and
// a normalized comparison. V2 omitted I/O, PDN and clock distribution outright.
// This gives each previously OPEN term a stated basis plus a declared intensity
// sweep, so the budget becomes a bounded interval instead of an unsized floor.
// Every intensity is DECLARED. None is measured.
static void scenario_power_economics() {
    std::cout << "\n[SCENARIO 8] V3 power economics: the four OPEN terms given a stated basis\n";

    const auto& c = params().control;
    const auto& g = params().gaps;
    const auto& a = params().arch;

    const double data_plane = c.data_plane_mW_per_cm2 * 1e-3 * c.data_plane_area_cm2;
    const double restoration = fea::restoration_power_W();
    const double refresh = (fea::payload_bits() / fea::refresh_interval_s()) *
                           g.refresh_energy_per_bit_J;
    const double known = data_plane + restoration + refresh;

    // External bandwidth is a choice of WHAT YOU CONNECT, not a fabric property.
    // Each row has a cited interface bandwidth and a cited or flagged pJ/bit.
    struct Attachment { const char* name; double GBps; double pJ; const char* note; };
    // NVMe PCIe 5.0 SSD: ~12 GB/s sequential (Crucial, MSI). pJ not cited.
    // LPDDR5X 8-channel: 17.1 GB/s per channel (Micron LPDDR5X datasheet) x 8.
    //   pJ taken as DIMM-class 12, cited from 3D-PATH.
    // HBM3 one stack: JEDEC max 819 GB/s. pJ taken as HBM 4, cited from 3D-PATH.
    const Attachment atts[3] = {
        {"NVMe SSD (PCIe 5.0)", 12.0, 12.0, "12 GB/s cited, pJ DECLARED"},
        {"LPDDR5X 8-channel", 136.0, 12.0, "17.1 GB/s/ch cited x8, pJ cited 12 (DIMM)"},
        {"HBM3 one stack", 819.0, 4.0, "819 GB/s JEDEC cited, pJ cited 4 (HBM)"},
    };
    // M10 finding: V2's claimed chip figure equals ONE Zone's rate, so the
    // fabric itself cannot deliver more than this until shared routing is fixed.
    const double fabric_ceiling_GBps = 1064.0;
    auto io_of = [&](double gbps, double pj) {
        const double capped = std::min(gbps, fabric_ceiling_GBps);
        return capped * 1e9 * 8.0 * pj * 1e-12;
    };
    const double io_lo = io_of(atts[0].GBps, atts[0].pJ);
    const double io_md = io_of(atts[1].GBps, atts[1].pJ);
    const double io_hi = io_of(atts[2].GBps, atts[2].pJ);

    // CITED clock share of dynamic power for clocked CMOS: 30-45%
    // (IJSR CSEIT23112577, 2025). Grade C venue, so the low end is the floor.
    // FZC removes per-Zone clocking, so this share applies ONLY to the boundary
    // ring, not to the whole fabric.
    const double clock_lo = 0.30, clock_md = 0.375, clock_hi = 0.45;

    // DECLARED fabric-side clock/bias. FZC replaced per-Zone clock generation,
    // so no citation exists for this path.
    const double fabric_lo = 0.05, fabric_md = 0.10, fabric_hi = 0.20;
    // DECLARED PDN loss fraction. Needs a floorplan to close.
    const double pdn_lo = 0.03, pdn_md = 0.07, pdn_hi = 0.12;
    // DECLARED boundary CMOS density. Needs the ring RTL to close.
    const double boundary_area = fea::boundary_ring_area_cm2(a.die_area_cm2);  // derived edge ring
    const double density_lo = 1.0, density_md = 10.0, density_hi = 50.0; // W/cm^2

    auto boundary_of = [&](double d) { return boundary_area * d; };
    auto total = [&](double density, double gbps, double pj, double cf, double ff, double pf) {
        const double b = boundary_of(density);
        double s = known + b + b * cf + data_plane * ff + io_of(gbps, pj);
        s += s * pf;  // power-delivery loss
        return s;
    };

    const double b_lo = boundary_of(density_lo);
    const double b_md = boundary_of(density_md);
    const double b_hi = boundary_of(density_hi);
    const double bclk_lo = b_lo * clock_lo;
    const double bclk_md = b_md * clock_md;
    const double bclk_hi = b_hi * clock_hi;
    const double fclk_lo = data_plane * fabric_lo;
    const double fclk_md = data_plane * fabric_md;
    const double fclk_hi = data_plane * fabric_hi;
    const double sub_lo = known + b_lo + bclk_lo + fclk_lo + io_lo;
    const double sub_md = known + b_md + bclk_md + fclk_md + io_md;
    const double sub_hi = known + b_hi + bclk_hi + fclk_hi + io_hi;
    const double pdn_lo_W = sub_lo * pdn_lo;
    const double pdn_md_W = sub_md * pdn_md;
    const double pdn_hi_W = sub_hi * pdn_hi;
    const double low = sub_lo + pdn_lo_W;
    const double mid = sub_md + pdn_md_W;
    const double high = sub_hi + pdn_hi_W;
    const double v2_claimed = 3.8;
    const double v2_corrected = reviewer_floor_W() + data_plane;
    (void)total; // arithmetic shown term by term above

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  external I/O, by attachment              \n";
    for (int i = 0; i < 3; ++i) {
        const double effective = std::min(atts[i].GBps, fabric_ceiling_GBps);
        std::cout << "    " << std::left << std::setw(22) << atts[i].name << std::right
                  << std::setprecision(0) << std::setw(7) << atts[i].GBps << " GB/s -> "
                  << std::setw(7) << effective << " effective x " << std::setprecision(1)
                  << atts[i].pJ << " pJ/bit = " << std::setprecision(3)
                  << io_of(atts[i].GBps, atts[i].pJ) << " W  [" << atts[i].note << "]\n";
    }
    std::cout << "  fabric delivery ceiling (M10): " << std::setprecision(0)
              << fabric_ceiling_GBps << " GB/s, so no attachment can draw more.\n\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  previously OPEN term   stated basis                         low       mid       high\n";
    std::cout << "  boundary CMOS          " << boundary_area << " cm^2 x density  "
              << std::setw(9) << b_lo << std::setw(10) << b_md << std::setw(10) << b_hi
              << "   W  DECLARED (1/10/50 W/cm^2)\n";
    std::cout << "  external I/O           by attached device          "
              << std::setw(9) << io_lo << std::setw(10) << io_md << std::setw(10) << io_hi
              << "   W  bandwidth CITED per attachment, pJ cited or flagged\n";
    std::cout << "  boundary-ring clock    frac of boundary ring   "
              << std::setw(9) << bclk_lo << std::setw(10) << bclk_md << std::setw(10) << bclk_hi
              << "   W  CITED 30/37.5/45% of dynamic power\n";
    std::cout << "  fabric clock and bias  frac of data plane       "
              << std::setw(9) << fclk_lo << std::setw(10) << fclk_md << std::setw(10) << fclk_hi
              << "   W  DECLARED 5/10/20%\n";
    std::cout << "  PDN losses             frac of delivered       "
              << std::setw(9) << pdn_lo_W << std::setw(10) << pdn_md_W << std::setw(10) << pdn_hi_W
              << "   W  DECLARED 3/7/12%\n";
    std::cout << "  known (data+restore+refresh)                              "
              << std::setw(9) << known << "   W  computed\n\n";
    std::cout << "  V3 TOTAL  low  : " << low << " W\n";
    std::cout << "  V3 TOTAL  mid  : " << mid << " W\n";
    std::cout << "  V3 TOTAL  high : " << high << " W\n\n";
    std::cout << "  contrast:\n";
    std::cout << "    V2 as written                 " << std::setw(9) << v2_claimed
              << " W   (wrong arithmetic, omits I/O, PDN, clock)\n";
    std::cout << "    V2 control plane, corrected   " << std::setw(9) << v2_corrected
              << " W   (reviewers' own figures + data plane)\n";
    std::cout << "    V3 mid, with I/O actually paid " << std::setw(9) << mid << " W\n";

    require(high > mid && mid > low, "the economics must open into a bounded interval");
    require(low > known, "adding the four terms must raise the floor");
    require(mid > v2_claimed,
            "V3's honest mid must exceed V2's 3.8 W claim, because V2 never paid for I/O");
    require(v2_corrected / mid > 100.0,
            "V3 must still be far below the corrected per-Zone CMOS control plane");
    require(high / low < 50.0,
            "pinning bandwidth to real attachments must narrow the interval below 50x");
    require(io_hi < 100.0, "an HBM attachment at cited pJ/bit must stay well under 100 W");
    std::cout << "\n" << std::setprecision(1);
    std::cout << "  sourced or attachment-pinned:\n";
    std::cout << "    I/O bandwidth        three cited interface rates\n";
    std::cout << "    I/O energy per bit   cited for LPDDR5X (12) and HBM3 (4)\n";
    std::cout << "    boundary-ring clock  cited 30-45% of dynamic power\n";
    std::cout << "  still DECLARED, no source found:\n";
    std::cout << "    SSD pJ/bit, fabric clock share, boundary CMOS density, PDN fraction\n";
    std::cout << "  interval span: " << (high / low) << "x, down from 326x before any sourcing.\n";
    std::cout << "  V3 mid is " << (mid / v2_claimed) << "x V2's CLAIMED 3.8 W, because V2 left\n";
    std::cout << "  I/O, PDN and clock distribution out of its own budget.\n";
    std::cout << "  V3 mid is " << (v2_corrected / mid) << "x BELOW corrected V2, which is the real\n";
    std::cout << "  advancement: FZC removes a kilowatt-scale control plane.\n";
    std::cout << "  finding: the defensible V3 power story is NOT lower than V2's claim. It is that\n";
    std::cout << "  V2's claim was unsound, and V3 removes the term that made it unsound.\n";
    std::cout << "  label: arithmetic derived; bandwidth and ring clock CITED; SSD pJ, fabric\n";
    std::cout << "  clock share, boundary density and PDN fraction DECLARED, not measured.\n";
}

} // namespace budget

int main() {
    using namespace budget;
    try {
        std::cout << "FEA V3 M1 reconciled control-plane budget\n";
        std::cout << "All terms computed from fea_params.h. Not a tape-out estimate.\n";
        scenario_unit_checks();
        scenario_decoder_area_exceeds_die();
        scenario_v2_literal_values_fail_checks();
        scenario_reviewer_hand_calculations();
        scenario_fzc_boundary_alternative();
        scenario_total_power();
        scenario_v3_control_budget();
        scenario_power_economics();
        std::cout << "\nPASS: V2 control-plane arithmetic reconciled and unit-checked.\n";
        std::cout << "LABEL: derived arithmetic, estimated device parameters, proposed control migration.\n";
        std::cout << "NEXT EVIDENCE GATE: size boundary CMOS, I-O, PDN, and clock-distribution power explicitly.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
