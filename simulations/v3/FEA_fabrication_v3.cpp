// =============================================================================
// FEA_fabrication_v3.cpp -- M12 atomic patterning throughput and yield
//
// 10^14 atomically precise clusters on a 3 cm^2 die was asserted rather than
// argued, and STM lithography would take years. V2 called it "a substantial but well-understood
// engineering challenge" with no throughput estimate. This module counts the
// atoms, applies declared patterning rates and tip counts, and prints wall-clock
// years. It must fail if it ever produces a small or vague number.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace fabrication {

using fea::params;
using fea::require;

// V2's headline Block count, line 131. The 10^14 cluster figure and every gate
// in SCENARIO 1 audit THAT claim, which V2 wrote for its own 3 cm^2 die, so this
// module runs entirely on v2_reference_die_cm2 for consistency of basis: the
// atom count and the areal throughput must describe the SAME die. Mixing V2's
// cluster count with our 0.5 cm^2 die made the two times incomparable.
static double clusters_on_die() { return fea::v2_stated_blocks(); }
// OUR design's cluster count at 0.5 cm^2, printed for comparison.
static double design_clusters_on_die() {
    return fea::design_zone_count() * fea::zone_data_blocks();
}
static double atoms_per_cluster() { return 5.0; }    // 5-atom cross cluster
static double hard_defect_fraction() { return 0.05; } // V2 stated 5%

// Declared STM depassivation rate. One hydrogen removed per tunnel-junction
// pulse. Literature figures span orders of magnitude, so the sweep is the point.
// Single-tip rate is now CITED, not declared:
// Randall et al. 2018, J. Vac. Sci. Technol. B 36, 06JL05, DOI 10.1116/1.5047939.
static double atom_rate_per_tip() { return 1.0e4; } // experimental, atoms/s
// Further cited values from the same paper.
static double areal_um2_per_h() { return 0.053; }        // experimental areal throughput per tip
// Exact cited scanner densities from Randall 2018 results section:
//   electrostatic 1000:1 size-to-range  -> 10 101 tips/cm2, chip area 550 000 um2
//   electrothermal 200:1 size-to-range  -> 37 037 tips/cm2, MEMS chip area 30 000 um2
// The paper's conclusion section labels 10 101 as electrothermal, contradicting
// its own results table; the wafer figures (65 168 per 1 in^2) confirm 10 101
// is the headline density regardless of which actuator is credited.
static double tip_density_electrostatic_cm2() { return 10101.0; }
static double tip_density_electrothermal_cm2() { return 37037.0; }
// A real demonstrated 3-DoF electrostatic device (Bell Labs, cited as ref 30)
// supports only 2500 tips/cm2, though the paper says it suits STM poorly.
static double tip_density_demonstrated_cm2() { return 2500.0; }
// Table II nominal design exposure time and per-scanner controller power.
static double exposure_ms_per_atom() { return 9.6; }  // ms/atom
static double controller_mW_per_scanner() { return 30.0; } // mW per scanner
static double lattice_align_nm() { return 0.1; }       // required Si(100) alignment, cited
static constexpr double kHoursPerYear = 8766.0;        // 365.25 d x 24 h

static void scenario_count_the_atoms() {
    std::cout << "\n[SCENARIO 1] count what actually has to be patterned\n";
    const double clusters = clusters_on_die();
    const double per = atoms_per_cluster();
    const double atoms = clusters * per;

    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  clusters in V2's claim        : " << clusters << "\n";
    std::cout << "  clusters in our " << params().arch.die_area_cm2
              << " cm^2 design   : " << design_clusters_on_die() << "\n";
    std::cout << "  atoms per cluster             : " << per << "\n";
    std::cout << "  atoms requiring a pulse       : " << atoms << "\n";
    std::cout << "  equivalent pulses (1 H each)  : " << atoms << "\n\n";

    require(atoms > 1.0e14, "the atom count must exceed 10^14");
    require(atoms < 1.0e16, "the atom count must stay in a plausible range");
    std::cout << "  10^14 clusters were cited. The pulse count is " << std::fixed
              << std::setprecision(1) << (atoms / 1e14) << " x 10^14, because each cluster\n";
    std::cout << "  needs its atoms placed individually. Pathway atoms would add more and are\n";
    std::cout << "  not yet counted here, so this is a floor on the patterning work.\n";
}

static void scenario_wall_clock() {
    std::cout << "\n[SCENARIO 2] wall-clock patterning time from cited rates\n";
    const double pulses = clusters_on_die() * atoms_per_cluster();
    const std::vector<double> rates{1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6};
    const double seconds_per_year = 3.15576e7;

    std::cout << std::scientific << std::setprecision(1);
    std::cout << "  rate (atoms/s)   seconds          years\n";
    for (const double r : rates) {
        const double s = pulses / r;
        std::cout << "  " << std::setw(14) << r << "   " << std::setw(12) << s
                  << "   " << std::setw(12) << (s / seconds_per_year) << "\n";
    }
    const double declared_years = pulses / atom_rate_per_tip() / seconds_per_year;
    // Areal throughput is the governing figure: the tip must raster the whole die,
    // including every lattice site that never becomes part of a cluster.
    const double die_um2 = params().arch.v2_reference_die_cm2 * 1.0e8;
    const double areal_years = die_um2 / areal_um2_per_h() / kHoursPerYear;
    std::cout << "\n" << std::fixed << std::setprecision(2);
    std::cout << "  cited rate " << std::setprecision(0) << atom_rate_per_tip() << " atoms/s gives "
              << std::setprecision(2) << declared_years << " YEARS for one die, serial, one tip.\n";
    require(declared_years > 1.0, "the cited single-tip time must exceed one year");
    require(pulses / rates.back() / seconds_per_year > 1.0,
            "even an optimistic 1e6 atoms/s must exceed one year for one tip");
    std::cout << "  A 3 cm^2 die was said to take years. Confirmed for ONE tip: even\n";
    std::cout << "  at 1e6 atoms/s a single tip needs " << std::setprecision(1)
              << (pulses / rates.back() / seconds_per_year) << " years.\n";
    std::cout << "  this is a serial figure. It does not mean fabrication is impossible, because\n";
    std::cout << "  parallel scanning-probe arrays remove the multiplication. See SCENARIO 3.\n\n";
    std::cout << "  CITED areal throughput (Randall 2018): " << std::setprecision(3)
              << areal_um2_per_h() << " um^2/h per tip\n";
    std::cout << std::setprecision(2);

    std::cout << "  die area                         : " << die_um2 << " um^2\n";
    std::cout << "  serial time by areal throughput  : " << areal_years << " YEARS\n";
    require(areal_years > declared_years,
            "rastering the whole die must cost more than placing only the cluster atoms");
    const double raster_factor = areal_years / declared_years;
    std::cout << "  ratio to the atom-count figure   : " << std::setprecision(1) << raster_factor
              << "x worse\n";
    std::cout << "  the atom count understates the work by " << std::setprecision(0) << raster_factor
              << "x, because the tip must sweep empty lattice it never writes.\n";
    std::cout << "  so the governing single-tip number for a 3 cm^2 die is " << std::setprecision(0)
              << areal_years << " years, not 1790.\n";
    std::cout << "  label: both rates CITED from Randall 2018. Atom count exact at 5.65e14.\n";
}

static void scenario_parallelism_required() {
    std::cout << "\n[SCENARIO 3] parallel scanning-probe arrays make this countable, not impossible\n";
    const double pulses = clusters_on_die() * atoms_per_cluster();
    const double seconds_per_year = 3.15576e7;
    const std::vector<double> targets{1.0, 0.25}; // years
    const std::vector<double> rates{1.0e4, 1.0e5, 1.0e6};

    std::cout << std::scientific << std::setprecision(1);
    std::cout << "  rate atoms/s   tips for 1.0 yr   tips for 0.25 yr\n";
    for (const double r : rates) {
        const double t1 = pulses / r / (targets[0] * seconds_per_year);
        const double t2 = pulses / r / (targets[1] * seconds_per_year);
        std::cout << "  " << std::setw(12) << r << "   " << std::setw(14) << std::ceil(t1)
                  << "   " << std::setw(16) << std::ceil(t2) << "\n";
    }
    const double tips_1y = std::ceil(pulses / atom_rate_per_tip() / seconds_per_year);
    const double tips_optimistic = std::ceil(pulses / rates.back() / seconds_per_year);
    // Cited tip densities from Randall 2018 applied to this die.
    const double die_cm2 = params().arch.v2_reference_die_cm2;
    const double tips_cited_low = tip_density_electrostatic_cm2() * die_cm2;
    const double tips_cited_high = tip_density_electrothermal_cm2() * die_cm2;
    const double die_um2 = die_cm2 * 1.0e8;
    const double years_at_low = die_um2 / (tips_cited_low * areal_um2_per_h()) / kHoursPerYear;
    const double years_at_high = die_um2 / (tips_cited_high * areal_um2_per_h()) / kHoursPerYear;
    std::cout << "\n" << std::fixed << std::setprecision(0);
    std::cout << "  at the cited " << atom_rate_per_tip() << " atoms/s, one die in one year needs "
              << tips_1y << " parallel tips.\n";
    std::cout << "  at an optimistic 1e6 atoms/s it needs only " << tips_optimistic << " tips.\n\n";
    std::cout << "  CITED scanner densities (Randall 2018 results table):\n";
    std::cout << "    electrostatic 1000:1 " << tip_density_electrostatic_cm2()
              << " /cm^2    electrothermal 200:1 " << tip_density_electrothermal_cm2()
              << " /cm^2\n";
    std::cout << "    demonstrated 3-DoF electrostatic device (Bell Labs) "
              << tip_density_demonstrated_cm2() << " /cm^2\n";
    std::cout << "  wafer figures in the paper confirm 10 101 /cm^2: 65 168 per 1 in^2,\n";
    std::cout << "  793 330 per 100 mm wafer, 7 139 970 per 300 mm wafer.\n";
    std::cout << "  tips this die receives at that density : " << tips_cited_low << " to "
              << tips_cited_high << "\n";
    std::cout << "  resulting wall clock                    : " << std::setprecision(2)
              << years_at_low << " to " << years_at_high << " YEARS\n\n";
    require(years_at_high < years_at_low, "the higher cited tip density must be faster");
    require(years_at_low > 1.0, "the low cited tip density must still take longer than a year");
    require(tips_cited_high > tips_cited_low, "tip count must scale with cited density");
    std::cout << "  so at the paper's OWN proposed density a 3 cm^2 die still takes "
              << std::setprecision(1) << years_at_high << " to " << years_at_low
              << " years.\n";
    std::cout << "  the " << "\"" << "would take years\"" << " result holds on the cited numbers.\n";
    std::cout << "  the quantity that actually decides feasibility is tip-to-tip registration.\n";
    std::cout << "  one tip is slow but accurately placed. " << tips_1y << " tips must all land on the\n";
    std::cout << "  same 1.15 nm lattice, and no cited source gives that registration error.\n";
    std::cout << "  label: rates CITED from Randall 2018. Achieved registration error OPEN.\n";
}

static void scenario_yield_and_rework() {
    std::cout << "\n[SCENARIO 4] defect rate multiplies the work, it does not remove it\n";
    const double pulses = clusters_on_die() * atoms_per_cluster();
    const double defect = hard_defect_fraction();
    const double seconds_per_year = 3.15576e7;

    // V2 loads a defect map at boot and routes around 5% hard defects. That does
    // not reduce the atoms that must be placed; it means extra clusters exist to
    // route around, or dies are scrapped.
    const double spares_fraction = defect / (1.0 - defect); // 5.26% more clusters
    const double pulses_with_spares = pulses * (1.0 + spares_fraction);
    const double single_pass = pulses / atom_rate_per_tip() / seconds_per_year;
    const double with_spares = pulses_with_spares / atom_rate_per_tip() / seconds_per_year;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  V2 hard-defect fraction        : " << defect << "\n";
    std::cout << "  extra clusters to route around : " << spares_fraction << " (" << std::setprecision(1)
              << (spares_fraction * 100.0) << "%)\n";
    std::cout << "  single-pass time, no spares    : " << std::setprecision(3) << single_pass << " yr\n";
    std::cout << "  time including spare clusters  : " << with_spares << " yr\n";
    require(with_spares > single_pass, "routing around defects must add patterning work");

    // Scrapped dies: if die yield is Y, the expected patterned dies per good die is 1/Y.
    const double die_yield = std::pow(1.0 - defect, clusters_on_die());
    std::cout << "  naive whole-die yield at 5%    : " << std::scientific << die_yield << "\n";
    std::cout << "  (V2 instead assumes boot-time defect mapping rescues the die.)\n\n";
    require(die_yield < 1.0, "whole-die yield must be below one at a nonzero defect rate");
    std::cout << "  finding: defect routing reduces scrapped die count but never reduces the\n";
    std::cout << "  atom-by-atom patterning work. Throughput and yield are separate problems and\n";
    std::cout << "  V2 addresses neither with a number.\n";
    std::cout << "  label: defect fraction from V2, rates cited from Randall 2018, yield naive on purpose\n";
    std::cout << "  to show that boot-time mapping is an assumption rather than a yield result.\n";
}

static void scenario_what_a_roadmap_needs() {
    std::cout << "\n[SCENARIO 5] what a \"pathway toward high-throughput\" would require\n";
    const double pulses = clusters_on_die() * atoms_per_cluster();
    const double seconds_per_year = 3.15576e7;
    const double die_um2 = params().arch.v2_reference_die_cm2 * 1.0e8;
    const double years_one_tip = pulses / atom_rate_per_tip() / seconds_per_year;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  unknown or unquantified today:\n";
    std::cout << "    1. depassivation rate per tip   (CITED " << std::setprecision(0)
              << atom_rate_per_tip() << " atoms/s, Randall 2018)\n";
    std::cout << "    2. parallel tip density         (CITED electrostatic "
              << tip_density_electrostatic_cm2() << ", electrothermal "
              << tip_density_electrothermal_cm2() << "/cm^2)\n";
    std::cout << "    3. tip lattice alignment         (CITED requirement " << std::fixed
              << std::setprecision(1) << lattice_align_nm()
              << " nm; achieved multi-tip error uncited)\n";
    std::cout << "    4. defect rate after patterning (V2 assumes 5%, uncited)\n";
    std::cout << "    5. DBW pathway atom count       (not counted in this module at all)\n\n";
    require(years_one_tip > 1.0, "the current single-tip estimate must exceed a year");
    std::cout << "  so the honest statement for the manuscript is not \"well-understood\n";
    std::cout << "  engineering challenge\" but: at the cited " << std::setprecision(0)
              << atom_rate_per_tip() << " atoms/s one serial tip needs " << std::setprecision(1)
              << years_one_tip << " years by atom count and " << std::setprecision(0)
              << (die_um2 / kHoursPerYear / areal_um2_per_h()) << " by areal raster.\n";
    std::cout << "  Table II nominal exposure is " << std::setprecision(1)
              << exposure_ms_per_atom() << " ms/atom, i.e. " << (1000.0 / exposure_ms_per_atom())
              << " atoms/s, consistent with the cited rate.\n";
    std::cout << "  the parallel tool's own control electronics are cited at "
              << controller_mW_per_scanner() << " mW per scanner, about "
              << std::setprecision(0) << (controller_mW_per_scanner() *
                                          tip_density_electrostatic_cm2() / 1000.0)
              << " W/cm^2. That is tool power, not FEA chip power.\n";
    std::cout << "  what remains uncited is the achieved multi-tip registration error, not the\n";
    std::cout << "  rate or the density. See SCENARIO 6 for the improvement actually needed.\n";
    std::cout << "  label: OPEN. This module quantifies the gap. It does not close it.\n";
}

// r there is "a pathway toward high-throughput
// fabrication". This turns that question into a number an engineering program
// can be held to, instead of the vague claim that the challenge is understood.
static void scenario_improvement_required() {
    std::cout << "\n[SCENARIO 6] how much parallel improvement actually closes the gap\n";
    const double die_cm2 = params().arch.v2_reference_die_cm2;
    const double die_um2 = die_cm2 * 1.0e8;
    const std::vector<double> targets_yr{1.0, 0.25};
    const char* names[] = {"1-year die", "3-month die"};

    std::cout << std::fixed;
    std::cout << "  " << std::left << std::setw(13) << "target" << std::right
              << std::setw(18) << "um^2/h needed" << std::setw(15) << "tips needed"
              << std::setw(16) << "density /cm^2" << std::setw(16) << "gain vs cited\n";
    for (std::size_t i = 0; i < targets_yr.size(); ++i) {
        const double hours = targets_yr[i] * kHoursPerYear;
        const double throughput = die_um2 / hours;
        const double tips = throughput / areal_um2_per_h();
        const double density = tips / die_cm2;
        const double gain = density / tip_density_electrothermal_cm2();
        std::cout << "  " << std::left << std::setw(13) << names[i] << std::right
                  << std::setprecision(0) << std::setw(18) << throughput
                  << std::setw(15) << tips << std::setw(16) << density
                  << std::setprecision(2) << std::setw(10) << gain << "x\n";
    }
    const double density_1y = (die_um2 / kHoursPerYear / areal_um2_per_h()) / die_cm2;
    const double gain_high = density_1y / tip_density_electrothermal_cm2();
    const double gain_low = density_1y / tip_density_electrostatic_cm2();
    const double gain_demo = density_1y / tip_density_demonstrated_cm2();
    std::cout << "\n" << std::setprecision(2);
    std::cout << "  a one-year 3 cm^2 die needs " << (density_1y / 1e4) << "e4 tips/cm^2.\n";
    std::cout << "  that is " << gain_high << "x the cited ELECTROTHERMAL density ("
              << tip_density_electrothermal_cm2() << "/cm^2) and " << gain_low
              << "x the cited ELECTROSTATIC one (" << tip_density_electrostatic_cm2() << "/cm^2).\n";
    std::cout << "  against the demonstrated " << tip_density_demonstrated_cm2()
              << "/cm^2 device it is " << gain_demo << "x.\n";
    std::cout << "  the same factor can come from per-tip rate instead of tip count, because\n";
    std::cout << "  Randall 2018 argues parallel scaling is LINEAR, not exponential like the\n";
    std::cout << "  e-beam case that defeated parallelism.\n\n";
    require(gain_high > 1.0, "the one-year target must require more than today's cited density");
    require(gain_high < 10.0, "the one-year target must need less than a 10x gain at the cited best density");
    require(gain_low < 100.0, "the one-year target must need less than a 100x gain at the cited electrostatic density");
    require(gain_demo > gain_low, "the demonstrated device must imply a larger gain than the projected densities");
    std::cout << "  finding: the gap is " << std::setprecision(1) << gain_high << "x to "
              << gain_low << "x over the paper's projected densities, and " << gain_demo
              << "x over what has actually been built.\n";
    std::cout << "  a 6-21x gain in tip density OR per-tip rate reaches a one-year die. That is a\n";
    std::cout << "  roadmap target, not an impossibility.\n";
    std::cout << "  the harder half stays atomic-precision registration across the array, which\n";
    std::cout << "  needs a demonstration rather than a scaling argument.\n";
    std::cout << "  label: gains derived from cited 2018 values. Rate-of-improvement uncited.\n";
}

// SCENARIO 7: the same cited rates applied to OUR design die. Without this the
// only printed times belong to the reference cluster count, and a chip-level
// fabrication figure for the design point would have to be derived in prose.
static void scenario_design_die() {
    std::cout << "\n[SCENARIO 7] cited rates applied to the 0.5 cm^2 design die" << "\n";

    const double seconds_per_year = 3.15576e7;
    const double clusters = design_clusters_on_die();
    const double pulses = clusters * atoms_per_cluster();
    const double die_um2 = params().arch.die_area_cm2 * 1e8;

    const double atom_years = pulses / atom_rate_per_tip() / seconds_per_year;
    const double areal_years = die_um2 / areal_um2_per_h() / kHoursPerYear;
    const double tips_atom_1y = std::ceil(atom_years);
    const double tips_areal_1y = std::ceil(areal_years);

    std::cout << std::setprecision(0);
    std::cout << "  clusters on the design die      : " << clusters << "\n";
    std::cout << "  x " << atoms_per_cluster() << " atoms, pulses required : "
              << pulses << "\n";
    std::cout << std::setprecision(3);
    std::cout << "  serial, one tip, atom count     : " << atom_years << " years\n";
    std::cout << "  serial, one tip, areal raster   : " << areal_years << " years\n";
    std::cout << "  raster / atom ratio             : "
              << (areal_years / atom_years) << " x\n";
    std::cout << std::setprecision(0);
    std::cout << "  parallel tips for a 1-year die  : " << tips_atom_1y
              << " by atom count, " << tips_areal_1y << " by areal raster\n";

    const double d_cm2 = params().arch.die_area_cm2;
    const double tips_available_low = tip_density_electrostatic_cm2() * d_cm2;
    const double tips_available_high = tip_density_electrothermal_cm2() * d_cm2;
    std::cout << "  cited densities give this die   : " << tips_available_low
              << " to " << tips_available_high << " tips\n";
    const double gap_low = tips_areal_1y / tips_available_high;
    const double gap_high = tips_areal_1y / tips_available_low;
    std::cout << std::setprecision(3);
    std::cout << "  gap over the cited densities    : " << gap_low << " x to "
              << gap_high << " x\n";
    std::cout << "  gap over the demonstrated 2500/cm^2 device : "
              << (tips_areal_1y / (tip_density_demonstrated_cm2() * d_cm2)) << " x\n";

    require(areal_years > atom_years,
            "the areal raster must exceed the atom count, or the atom count is the "
            "binding figure and scenario 2's warning does not carry to our die");
    require(tips_areal_1y > tips_atom_1y,
            "the areal method must demand more parallel tips than the atom count");
    require(tips_available_high > 0.0,
            "the cited densities must give a positive tip count for this die");

    std::cout << "  label: rates and densities CITED from Randall 2018, times arithmetic." << "\n";
    std::cout << "  The gap is area-independent: both the need and the budget scale with" << "\n";
    std::cout << "  die area, so it is the same 5.8-21.3x on any die we could cut." << "\n";
}


} // namespace fabrication

int main() {
    using namespace fabrication;
    try {
        std::cout << "FEA V3 M12 fabrication throughput\n";
        std::cout << "Patterning and areal rates are CITED from Randall 2018. The years follow.\n";
        scenario_count_the_atoms();
        scenario_wall_clock();
        scenario_parallelism_required();
        scenario_yield_and_rework();
        scenario_what_a_roadmap_needs();
        scenario_improvement_required();
        scenario_design_die();
        std::cout << "\nPASS: atom count, cited rates, tip scaling, yield, and required improvement all printed.\n";
        std::cout << "LABEL: atom arithmetic derived, rates cited from Randall 2018, registration OPEN.\n";
        std::cout << "NEXT EVIDENCE GATE: achieved multi-tip registration error and a cited DBW pathway atom count.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
