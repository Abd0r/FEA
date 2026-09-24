// =============================================================================
// FEA_crosstalk_v3.cpp -- M5 neighbouring-Block electrostatic interaction
//
// Reviewer 7 point 3: at 1.15 nm pitch, neighbouring dangling-bond structures
// interact by Coulomb coupling and tunnelling, and the single-Block model may
// not hold. This module places occupied and empty neighbours around a target
// Block, computes the induced charging-energy shift with a declared screened
// potential, and measures whether write selectivity and retention margin survive.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace crosstalk {

using fea::check_unit;
using fea::params;
using fea::require;

// Declared screening length in nm. Silicon surface screening is not a measured
// value here; the sweep below is the point of the module.
static double lambda_nm() { return 1.0; }

// Relative permittivity for the coupling between two surface dangling bonds.
// The previous value of 1.0 was VACUUM, which is physically indefensible for a
// silicon surface and inflated the whole shift by roughly 6x. 6.0 is the
// effective surface value, close to (eps_Si + 1)/2 = 6.35, and is still a
// DECLARED number. Bulk Si is 11.7. The verdict flips across this range, which
// is why scenario 5 sweeps it instead of fixing one value.
static double permittivity_default() { return 6.0; }
static double permittivity_bulk_si() { return 11.7; }

static double neighbour_shift_eV(double distance_nm, bool occupied, double screening_nm,
                                 double relative_permittivity = 6.0) {
    if (!occupied) return 0.0;
    if (distance_nm <= 0.0) throw std::runtime_error("neighbour distance must be positive");
    if (relative_permittivity <= 0.0) throw std::runtime_error("relative permittivity must be positive");
    // One elementary charge at the neighbour, screened by the medium.
    const double e2_over_4pieps_eV_nm = 1.44; // eV nm, vacuum value
    return (e2_over_4pieps_eV_nm / relative_permittivity) *
           std::exp(-distance_nm / screening_nm) / distance_nm;
}

// Nearest-neighbour shell geometry on the square surface lattice.
static std::vector<double> neighbour_distances_nm() {
    const double pitch = params().arch.block_pitch_nm;
    std::vector<double> d;
    d.push_back(pitch);                       // 4 nearest, omitted multiplicity for clarity
    d.push_back(pitch * std::sqrt(2.0));      // 4 diagonal
    d.push_back(2.0 * pitch);                 // 4 next-nearest
    d.push_back(pitch * std::sqrt(5.0));      // 8 further
    return d;
}

// Total shift over all occupied shells for one screening length and permittivity.
static double total_shift_eV(double screening_nm, double relative_permittivity) {
    double s = 0.0;
    const std::vector<double> dists = neighbour_distances_nm();
    const int mult[] = {4, 4, 4, 8};
    for (std::size_t i = 0; i < dists.size(); ++i)
        s += mult[i] * neighbour_shift_eV(dists[i], true, screening_nm, relative_permittivity);
    return s;
}

static void scenario_shift_magnitude() {
    std::cout << "\n[SCENARIO 1] occupied neighbours shift the target charging energy\n";
    const double Ec0 = params().device.Ec_eV;
    const double lam = lambda_nm();
    const double pitch = params().arch.block_pitch_nm;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Ec baseline                : " << Ec0 << " eV\n";
    std::cout << "  declared screening length  : " << lam << " nm\n";
    std::cout << "  Block pitch                : " << pitch << " nm\n\n";

    const std::vector<double> dists = neighbour_distances_nm();
    const char* shell[] = {"nearest (1x pitch)", "diagonal (sqrt2)", "next (2x pitch)", "far (sqrt5)"};
    double total_shift = 0.0;
    for (std::size_t i = 0; i < dists.size(); ++i) {
        const double shift = neighbour_shift_eV(dists[i], true, lam);
        // Multiplicity of each shell on a square lattice.
        const int mult[] = {4, 4, 4, 8};
        total_shift += mult[i] * shift;
        std::cout << "  " << std::left << std::setw(20) << shell[i] << std::right
                  << " shift " << std::setprecision(4) << shift << " eV"
                  << "  (x" << mult[i] << " = " << (mult[i] * shift) << ")\n";
    }
    const double Ec_shifted = Ec0 - total_shift;
    std::cout << "  sum over all occupied shells: " << std::setprecision(4) << total_shift << " eV\n";
    std::cout << "  Ec becomes                  : " << Ec_shifted << " eV  (from " << Ec0 << ")\n";
    require(total_shift > 0.0,
            "occupied neighbours must shift Ec by a definite non-zero amount");
    std::cout << "  label: coupling form declared, screening length and permittivity DECLARED.\n";
}

static void scenario_retention_margin() {
    std::cout << "\n[SCENARIO 2] neighbour occupation moves retention because Ec moves\n";
    const auto& d = params().device;
    const double lam = lambda_nm();
    const double total_shift = [&] {
        double s = 0.0;
        const std::vector<double> dists = neighbour_distances_nm();
        const int mult[] = {4, 4, 4, 8};
        for (std::size_t i = 0; i < dists.size(); ++i) s += mult[i] * neighbour_shift_eV(dists[i], true, lam);
        return s;
    }();

    auto tau = [&](double Ec_eV) {
        const double kT = fea::kB * d.temperature_K / fea::kEV;
        return 1.0 / (d.phonon_attempt_Hz * std::exp(-Ec_eV / kT));
    };
    const double tau_nominal = tau(d.Ec_eV);
    const double tau_neighbours = tau(d.Ec_eV - total_shift);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  tau with isolated Block     : " << (tau_nominal * 1e3) << " ms\n";
    std::cout << "  tau with all neighbours up  : " << (tau_neighbours * 1e3) << " ms\n";
    std::cout << "  ratio                        : " << std::setprecision(2)
              << (tau_neighbours / tau_nominal) << "x\n";
    require(tau_neighbours > 0.0 && tau_nominal > 0.0, "retention must be positive");
    const bool degraded = tau_neighbours < tau_nominal;
    std::cout << (degraded ? "  occupied neighbours shorten retention in this parameterization.\n"
                           : "  occupied neighbours lengthen retention in this parameterization.\n");
    std::cout << "  finding: retention is a function of the occupation of the surrounding array,\n";
    std::cout << "  not of one isolated 5-atom cluster. A single-cluster Ec measurement cannot\n";
    std::cout << "  set the array's retention. Reviewer 7 #3 is correct.\n";
}

static void scenario_write_selectivity() {
    std::cout << "\n[SCENARIO 3] write selectivity must survive a worst-case neighbour pattern\n";
    const double Ec0 = params().device.Ec_eV;
    const double lam = lambda_nm();

    // Best case: all neighbours empty. Worst case: all occupied, and they push
    // Ec down so a write intended for the target also becomes favourable.
    const double full_shift = [&] {
        double s = 0.0;
        const std::vector<double> dists = neighbour_distances_nm();
        const int mult[] = {4, 4, 4, 8};
        for (std::size_t i = 0; i < dists.size(); ++i) s += mult[i] * neighbour_shift_eV(dists[i], true, lam);
        return s;
    }();

    const double Ec_empty_neighbours = Ec0;
    const double Ec_full_neighbours = Ec0 - full_shift;
    // Selectivity margin: how far the target sits from a write threshold.
    // Declared write threshold is half the charging energy barrier.
    const double write_threshold = 0.5 * Ec0;
    const double margin_empty = Ec_empty_neighbours - write_threshold;
    const double margin_full = Ec_full_neighbours - write_threshold;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  declared write threshold    : " << write_threshold << " eV\n";
    std::cout << "  margin, neighbours empty    : " << margin_empty << " eV\n";
    std::cout << "  margin, neighbours occupied : " << margin_full << " eV\n";
    std::cout << "  margin lost to crosstalk    : " << (margin_empty - margin_full) << " eV\n";
    // margin_empty - margin_full is exactly the neighbour shift, which is a sum
    // of exp(-d/lambda)/d terms and therefore positive for any positive
    // permittivity and distance. The gate that used to sit here,
    // `margin_empty > 0.0`, asserted a constant, and its replacement
    // `(margin_empty - margin_full) > 0` asserted a sum of positive terms. Both
    // are ungatable, so both are removed rather than re-wrapped. The falsifiable
    // version of this claim is SCENARIO 5, where the margin SIGN flips at a
    // computable permittivity of 8.8488.
    const bool margin_survives = margin_full > 0.0;
    std::cout << "  worst-case margin survives? : " << (margin_survives ? "yes" : "NO\n");
    if (!margin_survives) {
        std::cout << "  FAIL-CLASS RESULT: crosstalk erases write selectivity under this screening.\n";
    } else {
        std::cout << "  margin remains positive, but it is eroded by crosstalk.\n";
    }
    std::cout << "  label: threshold is declared. The structural point is that margin is a\n";
    std::cout << "  function of the array state, so a coupled-block sweep is mandatory before\n";
    std::cout << "  any write-fidelity number can be quoted.\n";
}

static void scenario_sweep_screening() {
    std::cout << "\n[SCENARIO 4] screening length decides whether crosstalk matters at all\n";
    const double Ec0 = params().device.Ec_eV;
    const std::vector<double> lambdas{0.3, 0.5, 1.0, 2.0, 5.0};
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  lambda(nm)   full-neighbour shift(eV)   Ec(eV)   tau(ms)\n";
    const auto& d = params().device;
    double worst_shift = 0.0;
    for (const double lam : lambdas) {
        double s = 0.0;
        const std::vector<double> dists = neighbour_distances_nm();
        const int mult[] = {4, 4, 4, 8};
        for (std::size_t i = 0; i < dists.size(); ++i) s += mult[i] * neighbour_shift_eV(dists[i], true, lam);
        const double kT = fea::kB * d.temperature_K / fea::kEV;
        const double tau = 1.0 / (d.phonon_attempt_Hz * std::exp(-(Ec0 - s) / kT));
        std::cout << "  " << std::setw(8) << lam << "   " << std::setw(18) << s
                  << "   " << std::setw(6) << (Ec0 - s) << "   " << std::setw(9) << (tau * 1e3) << "\n";
        worst_shift = std::max(worst_shift, s);
    }
    require(worst_shift > 0.0, "the sweep must produce a nonzero shift at some screening length");
    std::cout << "  at weak screening the shift approaches or exceeds Ec itself.\n";
    std::cout << "  finding: without a measured surface screening length, crosstalk is unbounded\n";
    std::cout << "  in either direction. This is the parameter that must come from experiment\n";
    std::cout << "  or an atomistic calculation before single-Block results can be arrayed.\n";
}

// A previous revision fixed relative permittivity at 1.0, which is VACUUM. That
// manufactured the module's headline failure: the write margin was -2.5508 eV
// only because the medium was assumed to provide no screening at all. The sign
// of that margin changes inside the physically plausible range, so it must be
// swept rather than fixed.
static void scenario_permittivity_flips_verdict() {
    std::cout << "\n[SCENARIO 5] permittivity decides the headline verdict\n";
    const double Ec0 = params().device.Ec_eV;
    const double lam = lambda_nm();
    const double write_threshold = 0.5 * Ec0;
    const struct { double eps; const char* label; } cases[] = {
        {1.0, "vacuum, the previous value"},
        {2.0, "low surface bound"},
        {permittivity_default(), "effective surface (eps_Si+1)/2 approx 6.35"},
        {permittivity_bulk_si(), "bulk silicon"},
    };

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  write threshold              : " << write_threshold << " eV\n";
    std::cout << "  eps_r     label                              shift(eV)   margin(eV)  verdict\n";
    double previous_margin = 0.0;
    bool have_previous = false;
    bool flip_seen = false;
    double margin_default = 0.0, margin_bulk = 0.0;
    for (const auto& c : cases) {
        const double shift = total_shift_eV(lam, c.eps);
        const double margin = Ec0 - shift - write_threshold;
        if (std::abs(c.eps - permittivity_default()) < 1e-9) margin_default = margin;
        if (std::abs(c.eps - permittivity_bulk_si()) < 1e-9) margin_bulk = margin;
        std::cout << "  " << std::setw(7) << c.eps << "  " << std::left << std::setw(34) << c.label
                  << std::right << std::setw(9) << shift << std::setw(12) << margin << "    "
                  << (margin > 0.0 ? "PASS" : "FAIL") << "\n";
        if (have_previous && ((previous_margin > 0.0) != (margin > 0.0))) flip_seen = true;
        previous_margin = margin;
        have_previous = true;
    }

    // PR5/n2: the guard that stood here (default != 1.0) tested a pasted
    // constant, which only a header edit could fail. A range check against
    // physically plausible values for a hydrogen-terminated silicon surface is
    // falsifiable: vacuum, or anything above bulk silicon, now fails it.
    require(permittivity_default() > 1.5 && permittivity_default() < 12.0,
            "the default permittivity must sit in the plausible range for a "
            "hydrogen-terminated silicon surface, not at vacuum and not above bulk Si");
    require(flip_seen,
            "the write-margin sign must flip somewhere inside the plausible permittivity range");
    require(margin_default > -1.0,
            "at the effective-surface default the worst-case margin must not be catastrophically negative");
    require(margin_bulk > 0.0,
            "at bulk-silicon permittivity the worst-case write margin must be positive");

    std::cout << "\n" << std::setprecision(4);
    std::cout << "  the verdict is NOT a property of the geometry alone. It flips at about eps_r = "
              << (total_shift_eV(lam, 1.0) / write_threshold) << ".\n";
    std::cout << "  previous revision fixed eps_r = 1.0 (vacuum), which produced the headline\n";
    std::cout << "  FAIL-CLASS result of -2.5508 eV. That was an artifact, not physics.\n";
    std::cout << "  at the effective-surface default the margin is " << margin_default
              << " eV, and at bulk Si it is " << margin_bulk << " eV.\n";
    std::cout << "  so the honest statement is: crosstalk margin is bounded by how strongly the\n";
    std::cout << "  medium screens, and no surface permittivity has been measured for this\n";
    std::cout << "  structure. Reviewer 7 #3 stands, but for the screening constant, not for a\n";
    std::cout << "  predetermined failure.\n";
    std::cout << "  label: all eps_r values DECLARED. Sweep is the result. No verdict without one.\n";
}

} // namespace crosstalk

int main() {
    using namespace crosstalk;
    try {
        std::cout << "FEA V3 M5 neighbouring-Block crosstalk\n";
        std::cout << "Screening length is declared, not measured. The sweep is the result.\n";
        scenario_shift_magnitude();
        scenario_retention_margin();
        scenario_write_selectivity();
        scenario_sweep_screening();
        scenario_permittivity_flips_verdict();
        std::cout << "\nPASS: crosstalk shift, retention coupling, selectivity margin and permittivity sweep computed.\n";
        std::cout << "LABEL: derived coupling form, screening and permittivity both DECLARED, no fixed verdict.\n";
        std::cout << "NEXT EVIDENCE GATE: atomistic or measured surface screening length at 1.15 nm pitch.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
