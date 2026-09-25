// =============================================================================
// FEA_refresh_v3.cpp -- M13 refresh traffic converted to power
//
// r 7 point 4: retention sets a refresh contract
// that V2 called "< 1e-4 overhead". M4 derives the traffic it must pay; at the
// 0.5 cm^2 design point that is 89708 GB/s. This
// module converts that traffic into watts using a declared energy per rewritten
// bit, sweeps that energy because it is unmeasured, and reports the refresh
// power as a range with the unmeasured input named.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace refresh {

using fea::params;
using fea::require;
using fea::p_abs_single;
using fea::zone_count_stated;

static double refresh_power_W(double energy_per_bit_J) {
    const double interval = fea::refresh_interval_s();
    // refresh_bits(), not payload_bits(): the FZC leaks charge too.
    const double bits_per_s = fea::refresh_bits() / interval;
    return bits_per_s * energy_per_bit_J;
}

static void scenario_traffic() {
    std::cout << "\n[SCENARIO 1] retention sets the refresh traffic\n";
    const double interval = fea::refresh_interval_s();
    const double passes_per_s = 1.0 / interval;
    const double bits = fea::refresh_bits();
    const double bytes = bits / 8.0;
    const double GBps = bytes * passes_per_s / 1e9;
    const double tau = fea::kramers_tau_s(params().device.Ec_eV,
                                          params().device.phonon_attempt_Hz,
                                          params().device.temperature_K);
    require(interval > 0.0 && passes_per_s > 0.0, "refresh cadence must be positive");
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  tau (300 K)                 : " << (tau * 1e3) << " ms\n";
    std::cout << "  refresh interval (tau/2, V2) : " << (interval * 1e3) << " ms\n";
    std::cout << "  full passes per second       : " << std::setprecision(1) << passes_per_s << "\n";
    std::cout << "  payload (data only)          : " << std::setprecision(1)
              << (fea::payload_bits() / 1e12) << "e12 bits = "
              << (fea::payload_bits() / 8.0 / 1e12) << " TB\n";
    std::cout << "  plus FZC control state       : " << (fea::fzc_refresh_bits() * fea::design_zone_count() / 1e12)
              << "e12 bits (" << static_cast<long long>(fea::fzc_refresh_bits()) << "/Zone)\n";
    std::cout << "  TOTAL a pass must rewrite    : " << std::setprecision(1) << (bits / 1e12)
              << "e12 bits = " << (bytes / 1e12) << " TB\n";
    std::cout << "  refresh traffic              : " << std::setprecision(0) << GBps << " GB/s\n";
    // Threshold tracks the design die. At 3 cm^2 this was 533024 GB/s and the
    // gate read 1.0e5; at 0.5 cm^2 traffic is 89652 GB/s, still four orders
    // above any real port, so V2's < 1e-4 claim remains worth testing.
    require(GBps > 1.0e4, "refresh traffic must be large enough to make V2's < 1e-4 claim worth testing");
    std::cout << "  V2 reported refresh overhead as 1.1e-4. That is an event ratio, not a\n";
    std::cout << "  bandwidth. This GB/s figure is the STATE MAINTAINED per second; it does not\n";
    std::cout << "  cross the chip port. See SCENARIO 4, where the local duty cycle is derived.\n";
}

static void scenario_power_sweep() {
    std::cout << "\n[SCENARIO 2] refresh power is a range until energy per bit is measured\n";
    const double interval = fea::refresh_interval_s();
    const double bits_per_s = fea::refresh_bits() / interval;

    const std::vector<double> energies{1e-19, 1e-18, 1e-17, 1e-16};
    const char* basis[] = {"optimistic", "declared", "conservative", "pessimistic"};
    std::cout << std::scientific << std::setprecision(1);
    std::cout << "  energy/bit(J)  basis          refresh power(W)\n";
    double lo = 1e300, hi = 0.0;
    for (std::size_t i = 0; i < energies.size(); ++i) {
        const double p = bits_per_s * energies[i];
        lo = std::min(lo, p);
        hi = std::max(hi, p);
        std::cout << "  " << std::setw(12) << energies[i] << "  " << std::left << std::setw(12)
                  << basis[i] << std::right << std::setw(16) << p << "\n";
    }
    require(lo > 0.0 && hi > lo, "the sweep must produce a positive range");
    const double ratio = hi / lo;
    require(ratio > 100.0, "a 1000x sweep in energy per bit must open a wide power range");
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  range: " << lo << " W to " << hi << " W  (span " << std::setprecision(1)
              << ratio << "x)\n";
    std::cout << "  at the declared 1e-18 J the refresh term is " << std::scientific
              << std::setprecision(4) << refresh_power_W(params().gaps.refresh_energy_per_bit_J)
              << " W.\n";
    std::cout << "  label: energy per bit is UNMEASURED. This term cannot be closed by\n";
    std::cout << "  arithmetic alone, and it must not be reported as a point value.\n";
}

static void scenario_refresh_dominates_data_plane() {
    std::cout << "\n[SCENARIO 3] refresh power compared with the data plane\n";
    const auto& c = params().control;
    const auto& a = params().arch;
    const double data_plane = c.data_plane_mW_per_cm2 * 1e-3 * c.data_plane_area_cm2;
    const double refresh = refresh_power_W(params().gaps.refresh_energy_per_bit_J);
    const double restoration = [&] {
        // Mirrors M6: endpoints at the declared spacing and 1 nW each.
        const double edge = std::sqrt(a.die_area_cm2 * 1e-4);
        const double decay_um = 1.0;
        const double threshold = 0.10;
        const double spacing = decay_um * 1e-6 * std::log(1.0 / threshold);
        const int per_row = static_cast<int>(std::ceil(edge / spacing));
        const double n = static_cast<double>(per_row) * per_row;
        return n * 1e-9;
    }();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  data plane (V2)              : " << data_plane << " W\n";
    std::cout << "  restoration endpoints (M6)   : " << restoration << " W\n";
    std::cout << "  refresh (declared energy)    : " << refresh << " W\n";
    std::cout << "  per-Zone CMOS after FZC      : 0.000000 W\n";
    require(data_plane > 0.0 && refresh >= 0.0, "power terms must be positive");
    std::cout << "  refresh is " << std::setprecision(1)
              << (refresh / data_plane * 100.0) << "% of the data plane at the declared energy.\n";
    std::cout << "  label: one of these three is V2's own figure, two are declared.\n";
}

static void scenario_interval_sensitivity() {
    std::cout << "\n[SCENARIO 4] tighter retention margin multiplies refresh power\n";
    const double bits = fea::refresh_bits();
    const double tau = fea::kramers_tau_s(params().device.Ec_eV,
                                          params().device.phonon_attempt_Hz,
                                          params().device.temperature_K);
    const double energy = params().gaps.refresh_energy_per_bit_J;
    // Interval sensitivity. N is the divisor: interval = tau/N, and power is
    // proportional to N because a shorter interval means more passes. V2 states
    // N = 2 at lines 373 and 722. N = 3 is the conservative alternative; N = 1
    // is the least frequent and the riskiest, since it refreshes exactly at tau.
    const std::vector<double> divisors{2.0, 3.0, 1.0, 0.5};
    const char* who[] = {"V2 design (line 373)", "conservative", "refresh at tau", "slack"};
    std::cout << std::scientific << std::setprecision(1);
    std::cout << "  interval = tau/N   interval(s)    refresh power(W)   basis\n";
    for (std::size_t i = 0; i < divisors.size(); ++i) {
        const double interval = tau / divisors[i];
        const double p = (bits / interval) * energy;
        std::cout << "  N = " << std::setw(4) << divisors[i] << "     " << std::setw(12) << interval
                  << "   " << std::setw(14) << p << "   " << who[i] << "\n";
    }
    // The design policy comes from the shared definition so this module has one
    // cadence beside the parameter store; the tau/1 and tau/3 rows are sweep
    // variants that deliberately do not use it.
    const double p_tau_2 = (bits / fea::refresh_interval_s()) * energy;
    const double p_tau_3 = (bits / (tau / 3.0)) * energy;
    const double p_tau_1 = (bits / (tau / 1.0)) * energy;
    // PR5/S1: all four gates that stood here (two orderings, two ratios) were
    // dead: bits, tau and energy cancel in the ratios so p_tau_2/p_tau_1 is 2
    // identically, and the orderings hold whenever energy > 0, so no input could
    // violate them. Replaced by checks on inputs NOT shared with that identity,
    // plus a cross-check against the term recorded from M1.
    require(energy > 0.0 && bits > 0.0 && tau > 0.0,
            "refresh power inputs must be strictly positive, or the divisor sweep is vacuous");
    const double p_design = fea::refresh_power_W();
    require(p_design > 0.0 && p_design < 1.0,
            "refresh power at the declared energy must land inside the stated 1 W band; "
            "an order-of-magnitude jump means the ledger or the interval moved");
    require(std::fabs(p_design - fea::m1_printed_refresh_W()) /
                fea::m1_printed_refresh_W() < 0.05,
            "refresh power must agree with the term recorded from M1 within 5%, otherwise "
            "the budget and the refresh module are computing different things");
    require(fea::fzc_refresh_groups() >= 1.0,
            "the FZC must occupy at least one refresh group, or it is not maintained at all");
    std::cout << "  V2's tau/2 costs " << std::fixed << std::setprecision(1)
              << (p_tau_2 / p_tau_1) << "x tau/1; the conservative tau/3 costs "
              << (p_tau_3 / p_tau_2) << "x V2's policy.\n";
    std::cout << "  so the divisor, the refresh margin and the energy per bit are all load-bearing\n";
    std::cout << "  and all stated. That is why V3 reports a range, not V2's 1.1e-4.\n";
}

static double kT_eV() {
    const auto& d = params().device;
    return fea::kB * d.temperature_K / fea::kEV;
}

static double tau_for_Ec(double Ec_eV) {
    const auto& d = params().device;
    return fea::kramers_tau_s(Ec_eV, d.phonon_attempt_Hz, d.temperature_K);
}

// Full local refresh time for one Zone: one FIRE per Word, 64 bits in parallel.
// Expected FIRE latency, derived exactly as M9 derives it. Extracted so the
// data pass and the FZC self-refresh scenario cannot disagree about it.
static double expected_fire_ps() {
    const auto& d = params().device;
    const auto& a = params().arch;
    const double vg = 2.0 * d.t_hop_eV * fea::kEV * d.a_lattice_m / fea::kHbar;
    return ((a.segment_um * 1e-6) / vg / fea::kPS) / p_abs_single();
}

// One full maintenance pass: every data Word, plus the FZC's own Blocks.
// Data words = Words per Zone = 1024 by definition, plus the FZC groups.
static double local_refresh_time_s() {
    const double data_words = fea::zone_data_blocks() / 64.0;
    return (data_words + fea::fzc_refresh_groups()) * expected_fire_ps() * fea::kPS;
}

// The FZC is built from the same Fusion Blocks as the data, so its control state
// loses charge on the same timescale. payload_bits() excludes it because control
// is not storage, so without this scenario the controller is the one part of a
// Zone nobody maintains. Refreshing all 535 at once would blank the Zone ID,
// routing state and refresh pointer together, so it rotates Word by Word:
// ceil(535/64) = 9 groups, one at a time. No group ever holds the only copy of
// critical state, because the ledger triple-replicates every state group.
static void scenario_fzc_refreshes_itself() {
    std::cout << "\n[SCENARIO 6] who refreshes the FZC: the controller maintains itself\n";
    const double fire_ps = expected_fire_ps();
    const double groups = fea::fzc_refresh_groups();
    const double rotation_s = groups * fire_ps * fea::kPS;
    const double interval = fea::refresh_interval_s();
    const double tau = fea::kramers_tau_s(params().device.Ec_eV, params().device.phonon_attempt_Hz,
                                          params().device.temperature_K);
    const double full_pass_s = local_refresh_time_s();

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  FZC Blocks to maintain        : "
              << static_cast<long long>(fea::zone_fzc_blocks()) << "\n";
    std::cout << "  grouped into Words            : " << std::setprecision(0) << groups
              << " groups of 64, one group per FIRE\n";
    std::cout << std::setprecision(4);
    std::cout << "  expected FIRE latency         : " << fire_ps << " ps\n";
    std::cout << "  time for one full rotation    : " << (rotation_s * 1e9) << " ns\n";
    std::cout << "  refresh interval (tau/2)      : " << (interval * 1e3) << " ms\n";
    std::cout << "  FZC rotation duty             : " << std::scientific << std::setprecision(3)
              << (rotation_s / interval) << "\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  tau (retention)               : " << (tau * 1e3) << " ms\n";
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  longest possible staleness    : one rotation, " << (rotation_s / tau)
              << " of tau\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  full pass incl. FZC           : " << (full_pass_s * 1e9) << " ns\n";

    require(static_cast<int>(groups) * 64 >= static_cast<int>(fea::zone_fzc_blocks()),
            "the Word groups must cover every FZC Block");
    require(static_cast<int>(groups) * 64 - 64 < static_cast<int>(fea::zone_fzc_blocks()),
            "the group count must be minimal, with no padded group hiding spare capacity");
    require(rotation_s < interval, "one FZC rotation must finish inside one refresh interval");
    require(tau > 1000.0 * rotation_s,
            "self-refresh needs a rotation at least 1000x shorter than tau, or the "
            "controller is its own weakest link");
    require(fea::refresh_bits() > fea::payload_bits(),
            "the refresh contract must cover strictly more than the payload, or the FZC is left out");

    std::cout << "  answer: the FZC SCHEDULES refresh, a local FIRE primitive PERFORMS it,\n";
    std::cout << "  and the FZC's own Blocks are part of that same maintenance pass, rotated\n";
    std::cout << "  in " << std::setprecision(0) << groups << " groups so the controller never vanishes.\n";
    std::cout << "  critical state is triple-replicated by the ledger, so no single group\n";
    std::cout << "  ever holds the only copy of Zone ID, routing state or the refresh pointer.\n";
    std::cout << "  the infinite regress stops because staleness is bounded by a ROTATION,\n";
    std::cout << "  not by an interval: " << std::scientific << std::setprecision(3)
              << (rotation_s / tau) << " of tau worst case.\n";
    std::cout << "  at power-on, before any FZC state exists, the boundary boot controller\n";
    std::cout << "  initialises it; after that the fabric maintains itself locally.\n";
    std::cout << "  label: grouping and duty DERIVED from the ledger and M9 latency; the\n";
    std::cout << "  schedule-vs-primitive split is PROPOSED; the heartbeat source for a\n";
    std::cout << "  cold Zone with no valid FZC is OPEN.\n";
}

// CORRECTION. An earlier revision solved for the Ec that fits refresh under
// V2's advertised PORT. That was an artifact of assuming global refresh through
// shared routing, which M10 scenario 4 has withdrawn: refresh is local, so there
// is no port constraint. The real constraint is that one local Zone refresh
// must finish inside one refresh interval. That sets a much lower Ec floor.
// The remaining value of raising Ec is ENERGY, not feasibility.
static void scenario_charging_energy_target() {
    std::cout << "\n[SCENARIO 5] charging energy: local feasibility floor, and energy as the lever\n";
    const auto& d = params().device;
    const double local_time = local_refresh_time_s();
    // Feasibility requires one local Zone refresh to finish inside one interval.
    // The interval is tau/N with N = 2 per V2 lines 373 and 722, so tau must be
    // at least N x local_time.
    const double tau_min = fea::refresh_interval_divisor() * local_time;
    const double rate_max = 1.0 / tau_min;
    const double floor_eV = kT_eV() * std::log(d.phonon_attempt_Hz / rate_max);
    const double energy = params().gaps.refresh_energy_per_bit_J;

    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  one local Zone refresh takes  : " << local_time << " s\n";
    std::cout << "  smallest feasible interval    : " << tau_min << " s (must exceed that)\n";
    std::cout << "  Ec floor for local refresh    : " << std::fixed << std::setprecision(4)
              << floor_eV << " eV\n";
    std::cout << "  V2 design point               : " << d.Ec_eV << " eV, which is "
              << std::setprecision(2) << (d.Ec_eV / floor_eV) << "x above the floor\n\n";
    require(d.Ec_eV > floor_eV,
            "V2's design point must sit above the floor local refresh requires");

    const std::vector<double> ecs{0.40, 0.50, 0.55, 0.60, 0.65, 0.70, 0.82};
    std::cout << "  Ec(eV)   tau        interval     local duty   bit-writes/s   refresh power\n";
    double previous_writes = 1e300;
    for (const double ec : ecs) {
        const double tau = tau_for_Ec(ec);
        const double interval = tau / fea::refresh_interval_divisor();
        const double duty = local_time / interval;
        const double writes = fea::payload_bits() / interval;
        const double power = writes * energy;
        std::cout << "  " << std::fixed << std::setprecision(3) << std::setw(6) << ec
                  << "   " << std::scientific << std::setprecision(2) << std::setw(9) << tau
                  << "   " << std::setw(10) << interval << "   " << std::setw(9) << duty
                  << "   " << std::setw(12) << writes << "   " << std::setw(11) << power << " W";
        if (std::abs(ec - d.Ec_eV) < 1e-9) std::cout << "  <- V2 design";
        if (ec < floor_eV) std::cout << "  <- below floor, refresh cannot finish";
        std::cout << "\n";
        require(writes < previous_writes, "higher Ec must reduce refresh work monotonically");
        previous_writes = writes;
    }

    const double at_design = fea::payload_bits() / (tau_for_Ec(d.Ec_eV) / fea::refresh_interval_divisor());
    const double duty_design = local_time / (tau_for_Ec(d.Ec_eV) / fea::refresh_interval_divisor());
    std::cout << "\n" << std::fixed << std::setprecision(4);
    std::cout << "  at V2's design point: local duty " << std::scientific << duty_design
              << ", refresh power " << (at_design * energy) << " W at the declared energy.\n";
    require(duty_design < 1.1e-4, "at V2's Ec the local duty must sit under V2's own 1.1e-4 claim");

    std::cout << "  WITHDRAWN: an earlier revision solved for Ec = 0.82 eV to fit refresh under\n";
    std::cout << "  V2's advertised port. That target was an artifact of assuming refresh crosses\n";
    std::cout << "  shared routing, which M10 scenario 4 has withdrawn. There is no port constraint.\n\n";

    std::cout << "  so raising Ec is not about feasibility, it is about ENERGY. Each 0.05 eV\n";
    std::cout << "  cuts refresh work by roughly exp(0.05/kT) = " << std::setprecision(1)
              << std::fixed << (std::exp(0.05 / kT_eV())) << "x.\n";
    std::cout << "  from 0.65 to 0.70 eV, refresh power falls about " << std::setprecision(0)
              << (std::exp(0.05 / kT_eV())) << "x for 0.05 eV of extra barrier.\n";
    std::cout << "  two levers V2 does not discuss: refresh only the active region instead of the\n";
    std::cout << "  whole array, or relax tau/2 to tau/1 for 2x less work at 2x the risk.\n";
    std::cout << "  none of this repairs a desorbed atom. Refresh restores charge, not structure.\n";
    std::cout << "  separate finding: FEA-architecture.pdf lists tau = 362 s at Ec = 0.70 eV while\n";
    std::cout << "  its own f_K = 2.76 /s implies 0.362 s. That preprint row has a 1000x slip.\n";
    std::cout << "  label: floor and energy derived from V2's Kramers formula. Ec UNMEASURED.\n";
}

} // namespace refresh

int main() {
    using namespace refresh;
    try {
        std::cout << "FEA V3 M13 refresh power\n";
        std::cout << "Energy per rewritten bit is declared, not measured. The range is the result.\n";
        scenario_traffic();
        scenario_power_sweep();
        scenario_refresh_dominates_data_plane();
        scenario_interval_sensitivity();
        scenario_charging_energy_target();
        scenario_fzc_refreshes_itself();
        std::cout << "\nPASS: refresh power bounded, local feasibility floor and Ec energy lever derived.\n";
        std::cout << "LABEL: derived cadence and traffic, declared energy per bit, open measurement.\n";
        std::cout << "NEXT EVIDENCE GATE: measured write energy per bit for a Fusion Block at 300 K.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
