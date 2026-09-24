// =============================================================================
// fea_params.h -- single source of truth for FEA V3 simulation parameters
//
// No module declares its own copy of a physical or architectural constant.
// Units are carried explicitly so the unit checker can reject order-of-magnitude
// errors of the kind reviewers found in V2 (3.3 W vs 3300 W, 0.14 uW vs
// 137.85 uW, 12 um^2 x 1.7e9 Zones vs a 3 cm^2 die).
// =============================================================================

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

namespace fea {

// ---- physical constants (SI unless noted) ----
constexpr double kHbar = 1.054571817e-34;  // J s
constexpr double kB = 1.380649e-23;        // J / K
constexpr double kQ = 1.602176634e-19;     // C
constexpr double kEV = 1.602176634e-19;    // J per eV
constexpr double kMebitEV = 1e-3;          // eV per meV
constexpr double kPS = 1e-12;              // s per ps
constexpr double kUS = 1e-6;               // s per us
constexpr double kFEMTO = 1e-15;           // J per fJ
constexpr double kMICRO = 1e-6;            // m per um
constexpr double kCM2_PER_UM2 = 1e-8;      // cm^2 per um^2
constexpr double kUM2_PER_CM2 = 1e8;       // um^2 per cm^2

// ---- device physics ----
struct Device {
    double t_hop_eV = 0.020;        // DBW nearest-neighbour hopping, eV
    double t_cluster_eV = 0.015;    // cluster central hopping, eV
    double a_lattice_m = 0.384e-9;  // Si(100) surface lattice constant, m
    double Ec_eV = 0.65;            // charging energy, eV
    double temperature_K = 300.0;   // operating temperature, K
    double phonon_attempt_Hz = 1.59e12; // bulk Si optical phonon, Hz
    double gamma_lead_meV = 22.5;   // one-lead broadening from self-energy, meV
    double gamma_two_lead_meV = 45.0; // two-lead cross configuration, meV
    double gamma_v1_hardcoded_meV = 8.0; // v1 hardcoded value, kept only for comparison
    double v_gate_V = 0.75;         // CMOS rail, V
};

// ---- architecture ----
struct Arch {
    // OUR reference design die. Everything the V3 design sizes follows this.
    double die_area_cm2 = 0.5;      // cm^2 (design point, edge 0.7071 cm)
    // V2's own stated die. V2's claims (1.7e9 Zones, 1.13e14 Blocks, 14.1 TB,
    // 12 um^2 per Zone) were written for THIS area, so every gate that audits
    // V2 must compare against it. Mixing V2's zone count with our die yields
    // ratios true of neither, which is why M2 reads this instead of
    // die_area_cm2.
    double v2_reference_die_cm2 = 3.0;
    // V2 line 116: practical density is stated after 2x routing overhead.
    double routing_area_overhead = 2.0;
    int blocks_per_zone = 65536;    // 256 x 256 data array
    int zone_side_blocks = 256;
    int word_bits = 64;             // V2 line 138: 64-bit Word, 64 blocks along a DBW
    double block_pitch_nm = 1.15;   // centre-to-centre spacing, nm
    // V2 states "A Fusion Block stores a bit", and its own arithmetic agrees:
    // 1.13e14 Blocks x 1 bit = 14.1 TB. The 16-bit figure came from the older
    // architecture paper, which used a different 16-atom Block. This is the
    // field M2 actually reads for density; a duplicate `int bits_per_block = 1`
    // was removed by peer review 4 as dead.
    double block_bits = 1.0;        // bits per Fusion Block, per V2's capacity math
    double zone_addressed_mm = 0.1; // intra-zone crossbar span, mm
    double segment_um = 1.0;        // FIRE segment length, um
};

// ---- CMOS control plane (V2 assumptions under review) ----
struct Control {
    int decoder_transistors = 3000;
    double decoder_area_um2 = 12.0;     // um^2 per Zone decoder
    double decoder_event_J = 15e-15;    // J per decode event
    double pll_group_K = 256;           // Zones per PLL
    double pll_power_W = 0.5e-3;        // W per PLL instance
    int sequencer_transistors = 500;
    int sense_transistors_per_word = 200;
    int words_per_zone = 1024;
    double sense_zone_W = 0.1e-6;       // W per Zone sequencer + sensing
    double data_plane_mW_per_cm2 = 26.47; // mW/cm^2 data plane
    double data_plane_area_cm2 = 0.5;     // data plane occupies our design die
};

// ---- timing ----
struct Timing {
    double t_arm_ps = 33.0;
    double t_fire_ps = 42.9;
    double t_confirm_ps = 33.0;
    double v_signal_frac_c = 0.1;
};

// ---- IO and pathways ----
struct Io {
    double v_bias_V = 0.010;  // V
    double i_pathway_A = 1e-6; // A placeholder, declared not measured
    int n_path = 1;
};

// ---- total-power terms that V2 omitted ----
// Refresh only lives here as a field. The four V2 omissions became functions
// below, because three of them depend on the die and a struct field would have
// to be a literal that goes stale whenever the design area changes.
struct PowerGaps {
    // Energy to rewrite one bit during refresh. DECLARED from the order of a
    // single CMOS-equivalent write, NOT a measured FEA value.
    double refresh_energy_per_bit_J = 1e-18;
};

struct Params {
    Device device;
    Arch arch;
    Control control;
    Timing timing;
    Io io;
    PowerGaps gaps;
};

inline const Params& params() {
    static const Params p;
    return p;
}

// ---- unit-checked helpers ----

// Check a derived value against an expected magnitude and unit name.
// `relative_tolerance` guards only against arithmetic slip, not against the
// modelling assumption behind the value.
inline double check_unit(const std::string& name, double value, const std::string& unit,
                         double expected_magnitude, double relative_tolerance = 0.05) {
    if (!(value == value)) throw std::runtime_error(name + " is NaN (unit " + unit + ")");
    const double scale = std::abs(expected_magnitude) > 0 ? std::abs(value / expected_magnitude) : 0.0;
    if (std::abs(scale - 1.0) > relative_tolerance) {
        throw std::runtime_error(
            "UNIT CHECK FAILED: " + name + " = " + std::to_string(value) + " " + unit +
            ", expected order " + std::to_string(expected_magnitude) + " " + unit +
            " (ratio " + std::to_string(scale) + "). This is the class of error reviewers found in V2.");
    }
    return value;
}

inline void require(bool ok, const std::string& what) {
    if (!ok) throw std::runtime_error("ASSERTION FAILED: " + what);
}

// joules per decode event from fJ
inline double decode_event_J() { return params().control.decoder_event_J; }

// Zones on the die, from die area and per-Zone raw footprint.
inline double raw_block_pitch_m() { return params().arch.block_pitch_nm * 1e-9; }

inline double raw_block_area_cm2() {
    const double m = raw_block_pitch_m();
    return m * m * 1e4; // m^2 to cm^2
}

// Kramers retention from the stated charging energy, attempt frequency, and
// temperature. Shared by M4 and M13 so refresh cannot diverge from retention.
inline double kramers_tau_s(double Ec_eV, double attempt_Hz, double temperature_K) {
    if (temperature_K <= 0.0 || attempt_Hz <= 0.0) throw std::runtime_error("bad retention input");
    const double kT_eV = kB * temperature_K / kEV;
    return 1.0 / (attempt_Hz * std::exp(-(Ec_eV / kT_eV)));
}

// Refresh interval divisor. V2 states its own policy explicitly at two places:
//   line 373: "periodic refresh at tau_ret/2 = 26.1 ms intervals"
//   line 722: "20 epochs of tau_ret/2 = 26.1 ms"
// So V2 refreshes at HALF the retention time, not a third. This module follows
// V2. An earlier revision used /3, which disagreed with the paper under review
// by 1.5x. M13 sweeps the divisor as a sensitivity.
inline double refresh_interval_divisor() { return 2.0; }

inline double refresh_interval_s() {
    return kramers_tau_s(params().device.Ec_eV,
                         params().device.phonon_attempt_Hz,
                         params().device.temperature_K) / refresh_interval_divisor();
}

// V2-audit claims. These belong to V2's paper and to V2's 3 cm^2 die. They are
// used ONLY where a module reproduces V2's own arithmetic, which the design
// contract requires (M1 hand calcs, M2 floorplan, M10 aggregate, M12 atom
// count). Nothing about our design is derived from them.
inline double zone_count_stated() { return 1.7e9; }  // V2 line 130
inline double v2_stated_blocks() { return 1.13e14; } // V2 line 131, a CLAIM

// Data Blocks in one Zone: 1024 Words x 64 bits x 1 Block/bit (V2 lines 126-127).
inline double zone_data_blocks() {
    return static_cast<double>(params().control.words_per_zone) *
           static_cast<double>(params().arch.word_bits);
}

// FZC Blocks resident in one Zone. ONE definition for M2, M15 and fzc-floorplan.
// 535 is the DECLARED default allocation: state 336 (7 groups x 8 bits x 2 rails
// x 3 replicas) + commands 130 + ports 16 + pathways 4 + spares 49. Confirmed by
// spec/REVIEWER-VALIDATION-MATRIX.md ("a declared ledger now gives FZC 535
// Blocks"), spec/REVIEWER-POINT-BY-POINT.md ("FZC budget at 535 Blocks per
// Zone") and FEA_fzc_floorplan_v3 SCENARIO 2. It deliberately EXCEEDS the
// provisional 512 target, which DECISIONS.md records as provisional.
// The 350 figure is NOT the design value: it is the conditional 4-bits-per-group
// variant that only exists to show 512 depends on state width.
inline int zone_fzc_blocks() { return 535; }

// Zone total: the 256 x 256 data array plus the FZC strip. This preserves V2's
// "Zone = 1024 Words" and "A Zone addresses 1,024 Words" (line 190) while FZC
// takes over the per-Word CMOS control V2 placed at line 208, outside the data.
inline double zone_total_blocks() { return zone_data_blocks() + zone_fzc_blocks(); }

// OUR design's Zone count, derived from OUR die. Not V2's 1.7e9, which belongs
// to V2's 3 cm^2 die. At 0.5 cm^2 with a 65886-Block Zone this gives 2.8691e8.
inline double design_zone_count() {
    const auto& a = params().arch;
    return a.die_area_cm2 /
           (zone_total_blocks() * raw_block_area_cm2() * a.routing_area_overhead);
}

inline double payload_bits() { return design_zone_count() * zone_data_blocks(); }

// FZC Blocks are made of the SAME Fusion Blocks, so their control state leaks
// charge exactly like data. payload_bits() deliberately excludes them, because
// control is not storage and capacity must not count them, which means the
// refresh contract has to add them back explicitly. Without this the controller
// is the one part of a Zone nobody maintains, and it would silently lose its
// Zone ID, routing state and refresh pointer inside a few tau.
inline double fzc_refresh_bits() { return zone_fzc_blocks(); } // 535, one bit each
inline double refresh_bits_per_zone() { return zone_data_blocks() + fzc_refresh_bits(); }

// What a full maintenance pass must rewrite: data plus FZC. Capacity metrics
// keep using payload_bits(); anything about refresh uses this.
inline double refresh_bits() { return design_zone_count() * refresh_bits_per_zone(); }

// Grouped Word by Word so the controller never disappears mid-refresh. A Word is
// 64 Blocks, so 535 needs ceil(535/64) = 9 groups. Refreshing them one at a time
// keeps the longest staleness to one rotation, not to one interval.
inline double fzc_refresh_groups() { return std::ceil(static_cast<double>(zone_fzc_blocks()) / 64.0); }

// ---- peer recovery: FZC-v0 "Declaring an FZC dead" and invariant 12 ----
// Retention tau on its own, so the recovery gate compares against the same
// Kramers number M4 and M13 use. refresh_interval_s() above is tau/2.
inline double retention_tau_s() {
    return kramers_tau_s(params().device.Ec_eV,
                         params().device.phonon_attempt_Hz,
                         params().device.temperature_K);
}

// The recovery peer set is N = 5 (N, S, E, W plus one regional peer) and the
// rule is 2/3 of it rounded up, which resolves to 4. Declared by FZC-v0, so it
// is a decision rather than a measurement.
inline double recovery_peer_count() { return 5.0; }
inline double recovery_quorum() { return std::ceil((2.0 / 3.0) * recovery_peer_count()); }

// A peer reports only after m consecutive silent heartbeats. DECLARED: no
// measurement exists, and FZC-v0 lists m alongside H in Open parameters.
inline double heartbeat_miss_threshold_m() { return 3.0; }

// The controller-state ledger that zone_fzc_blocks() documents as
// 336 = 7 groups x 8 bits x 2 rails x 3 replicas, split into its parts so
// recovery can say what a neighbour actually has to transmit.
inline double fzc_state_groups() { return 7.0; }
inline double fzc_state_bits_per_group() { return 8.0; }
inline double fzc_state_rails() { return 2.0; }
inline double fzc_state_replicas() { return 3.0; }

// A neighbour rebuilds a dead controller from ONE replica, so the 3x
// replication is insurance and not traffic. This is the seed.
inline double fzc_seed_bits() {
    return fzc_state_groups() * fzc_state_bits_per_group() * fzc_state_rails();
}
inline double fzc_ledger_bits() { return fzc_seed_bits() * fzc_state_replicas(); }

// Zone id + epoch + sequence + status. DECLARED framing, not a packet spec.
inline double heartbeat_packet_bits() { return 64.0; }

// One recovery report over Slingshot. DECLARED and UNSOURCED. M17 sweeps it
// across six decades and reports how far the gate moves, so that an input with
// no measurement behind it does not get to carry the result.
inline double slingshot_message_latency_s() { return 1.0e-6; }
inline bool slingshot_latency_sourced() { return false; }

// ---- steady-state thermal (M19) ----
// Forward declarations: these power terms are defined further down, and the
// thermal helpers below compose them, so they must be visible here.
inline double restoration_power_W();
inline double boundary_ring_W();
inline double clock_distribution_W();
inline double external_io_W();

// Where each power term actually DISSIPATES. Split by physical location rather
// than by the floor/declared split used for cost, because heat cares about
// placement: the array terms spread over the die, the boundary terms
// concentrate in the perimeter ring where pads and PHY live.
// In-fabric, distributed across the die. These three reproduce M1's printed
// FLOOR of 0.023384 W, which M19 cross-checks as a one-definition gate.
inline double refresh_power_W() {
    return refresh_bits() / refresh_interval_s() * params().gaps.refresh_energy_per_bit_J;
}
inline double array_power_W() {
    const auto& c = params().control;
    return c.data_plane_mW_per_cm2 * 1e-3 * c.data_plane_area_cm2 +
           restoration_power_W() + refresh_power_W();
}
// Boundary and periphery: ring, clock/bias, and the on-die share of I/O.
// io_ondie_fraction is DECLARED because a DIMM-class pJ/bit figure includes the
// package, the trace and the receiver on the other end, none of which heats our
// die. Swept by M19 rather than assumed.
inline double io_ondie_fraction() { return 0.25; }
inline double boundary_power_W() {
    return boundary_ring_W() + clock_distribution_W() +
           io_ondie_fraction() * external_io_W();
}

// Bulk-silicon thermal conductivity at 300 K. Standard reference value, cited.
inline double si_thermal_conductivity_W_per_mK() { return 148.0; }
// Die thickness. DECLARED: standard 300 mm wafer stock, no measurement of this
// device exists, and M19 sweeps it.
inline double die_thickness_m() { return 775e-6; }
// Package and spreader are not modelled at all; the back face is treated as an
// ideal sink at T0, which UNDERSTATES delta T and is flagged as such.
inline bool thermal_geometry_sourced() { return false; }

// ---- recorded outputs of M1, held here as single definitions ----
// M1 assembles these from a different term set than the modules that cross-check
// them, which is what makes the check worth having. They are recorded rather
// than recomputed because re-deriving them here would make the check tautuous.
// If M1's floor moves, these must move with it; a gate in M19 detects the drift.
inline double m1_printed_floor_W() { return 0.023384; }        // M1 SCENARIO 7, FLOOR (4 of 8)
inline double m1_printed_refresh_W() { return 7.1766e-4; }     // M1 SCENARIO 7, refresh term
inline double per_zone_cmos_control_floor_W() { return 3728.326; } // M1, reviewers' figures + data plane

// ---- rescue path: where the recovery receiver's area lives (M18) ----
// The FZC ledger that zone_fzc_blocks() documents as
// 336 state + 130 command + 16 port + 4 pathway + 49 spare = 535, exposed as
// parts so M18 can gate that the documented split still sums, and so the port
// Blocks a rescue receiver would OCCUPY are counted rather than assumed.
inline double fzc_command_blocks() { return 130.0; }
inline double fzc_port_blocks() { return 16.0; }
inline double fzc_pathway_blocks() { return 4.0; }
inline double fzc_spare_blocks() { return 49.0; }
inline double fzc_ledger_sum() {
    return fzc_ledger_bits() + fzc_command_blocks() + fzc_port_blocks() +
           fzc_pathway_blocks() + fzc_spare_blocks();
}

// Per-transistor area implied by V2's OWN declared decoder: 12 um^2 for 3000
// transistors. Derived once here so no module copies the pair. It is V2's
// declared number, so it is a claim rather than a PDK figure.
inline double tx_area_um2() {
    return params().control.decoder_area_um2 /
           static_cast<double>(params().control.decoder_transistors);
}

// Rescue tag detector width, in transistors. DECLARED: no measurement and no
// PDK, and M18 sweeps it rather than endorsing it.
inline double rescue_tag_transistors() { return 20.0; }

// Wire pitch for a dedicated boundary access tree. DECLARED, and swept across
// four decades, because no public 2 nm PDK exists to source it.
inline double rescue_wire_pitch_um() { return 0.1; } // 100 nm
inline bool rescue_geometry_sourced() { return false; }

// V2's single-FIRE absorption probability. ONE definition for M7, M9 and M11,
// which previously each carried their own copy of 0.46.
inline double p_abs_single() { return 0.46; }

// SECDED syndrome + correct pass on the critical path. ONE definition for M9
// and M11, which previously each carried their own copy of 5.0 ps.
inline double t_secded_ps() { return 5.0; }

// Boundary ring. ONE declared input: a fixed PHY THICKNESS, because the
// boundary pads do not shrink when the die shrinks. The FRACTION it costs is
// therefore DERIVED and rises as the die gets smaller: 5.02% at V2's 3 cm^2 but
// 12.06% at our 0.5 cm^2 design point. Budgeting a flat 5% at the design die
// understated the cost by 2.4x, so every consumer derives from this width.
// (An earlier revision had M15 back-solve 0.02166 cm to land on a declared 5%,
// which made M15 and M2 agree by construction rather than by choice.)
inline double boundary_ring_width_cm() { return 0.022; }

inline double boundary_ring_area_cm2(double die_area_cm2) {
    const double side = std::sqrt(die_area_cm2);
    const double w = boundary_ring_width_cm();
    return 4.0 * w * side - 4.0 * w * w;
}

inline double boundary_ring_fraction_at(double die_area_cm2) {
    return boundary_ring_area_cm2(die_area_cm2) / die_area_cm2;
}

// ---- the four terms V2 omitted: DECLARED, each still unsourced ----
// M1 scenario 8 gives each one a stated basis and sweeps it into an interval.
// These are the declared MID-POINT that M1 and M11 both read, so there is one
// declaration rather than two. Derived where they depend on the die, so they
// cannot go stale when the design area changes. None is a measurement.
inline double boundary_ring_density_W_per_cm2() { return 10.0; }   // declared, W/cm^2
inline double boundary_ring_clock_share() { return 0.375; }        // cited 30-45% mid
inline double fabric_clock_share() { return 0.10; }                // declared
inline double external_io_bandwidth_GBps() { return 136.0; }       // LPDDR5X 8-channel, cited
inline double external_io_pJ_per_bit() { return 12.0; }            // DIMM-class, cited
inline double fabric_bandwidth_ceiling_GBps() { return 1064.0; }   // M10 finding
inline double pdn_loss_fraction() { return 0.07; }                 // declared

// Boundary ring power: derived ring AREA times declared density.
inline double boundary_ring_W() {
    return boundary_ring_area_cm2(params().arch.die_area_cm2) * boundary_ring_density_W_per_cm2();
}

// Clock and bias: cited share of boundary power plus a declared fabric share.
inline double clock_distribution_W() {
    const auto& c = params().control;
    const double data_plane = c.data_plane_mW_per_cm2 * 1e-3 * c.data_plane_area_cm2;
    return boundary_ring_W() * boundary_ring_clock_share() + data_plane * fabric_clock_share();
}

// External I/O: attachment bandwidth times cited pJ/bit, capped by M10's ceiling.
inline double external_io_W() {
    const double bw =
        external_io_bandwidth_GBps() < fabric_bandwidth_ceiling_GBps()
            ? external_io_bandwidth_GBps()
            : fabric_bandwidth_ceiling_GBps();
    return bw * 1e9 * 8.0 * external_io_pJ_per_bit() * 1e-12;
}

// Sourcing flags. FALSE means declared but not sourced, and only a real source
// flips these. Openness is keyed on the flag, never on the value being zero, so
// declaring a number can never make the budget look complete by accident.
inline bool boundary_ring_sourced() { return false; }      // needs ring RTL
inline bool external_io_sourced() { return false; }        // needs a package model
inline bool clock_distribution_sourced() { return false; } // needs ring RTL
inline bool pdn_loss_sourced() { return false; }           // needs a floorplan

inline int open_power_term_count() {
    return (boundary_ring_sourced() ? 0 : 1) + (external_io_sourced() ? 0 : 1) +
           (clock_distribution_sourced() ? 0 : 1) + (pdn_loss_sourced() ? 0 : 1);
}

// Reference cells for M11's comparison column. Sourced, not invented.
// Chatterjee et al., "Architecting an Energy-Efficient DRAM System For GPUs",
// HPCA 2017, local PDF PDFs/Chatterjee2017-energy-efficient-DRAM-GPU.pdf,
// extracted with `pdftotext -q` to .read/chatterjee.txt (1767 lines):
//   line 1200: "estimated to be 112 fJ/bit or 1.8nJ for a 2KB row"
//   line 1202: "Previous work ... report even higher values for row-energy
//              (5-6 nJ per 2KB)"  -> 305-366 fJ/bit over a 16384-bit row
//   line 841-843: "for HBM, the column-energy ... can vary between 1.5 pJ/bit
//              when there is no toggling and 5.7 pJ/bit ... 100%"
// Line numbers depend on the pdftotext extraction, so each is paired with its
// verbatim quote: a different extraction shifts lines, not wording.
// These are DRAM-die figures under the paper's own boundary. They are NOT yet
// comparable with FEA's data-plane-only number, which is why M11 still refuses.
inline double dram_row_energy_J_per_bit() { return 112e-15; }        // Chatterjee HPCA 2017
inline double dram_row_energy_prior_low() { return 5e-9 / 16384.0; }  // 5 nJ / 2KB row
inline double dram_row_energy_prior_high() { return 6e-9 / 16384.0; } // 6 nJ / 2KB row
// CITED reference cells from open-access sources. Each carries the boundary it
// was measured under, because FEA's own metric is whole-die retention energy
// and these are mostly cell-write or macro figures. The metric mismatch is the
// reason M11 keeps refusing rather than a gap we can paper over.
// CIM cell write energy, arXiv 2406.08413 Table II (CC-BY), Verbatim row:
// "Write energy: SRAM <0.1 nJ; ReRAM 2 nJ; PCM 6 nJ; FeFET 0.1 J; MRAM <1 nJ".
// Boundary: single cell write, not a macro and not retention.
inline double cim_cell_write_energy_J_per_bit() { return 2e-9; } // ReRAM
inline double pcm_cell_write_energy_J_per_bit() { return 6e-9; }
// RC cell write, Frontiers in Neuroscience 2015 (open access):
// "energies to fully write a resistive memory cell as low as 6 fJ have been
// demonstrated (Cheng et al., 2010)". Boundary: single crossbar cell.
inline double rc_cell_write_energy_J_per_bit() { return 6e-15; }
// CIM macro storage density, Frontiers in Science 2025 (open access), citing
// Spetalnick ISSCC 2022: "A 40nm 64KB ... 2.37 Mb/mm2 RRAM ... macro".
// Boundary: whole macro including periphery, which matches FEA's whole-die basis.
inline double cim_macro_density_bits_per_mm2() { return 2.37e6; }
// HBM3 die area, Research Square preprint 2026, NOT peer reviewed and the die
// figure is secondhand: "HBM3 DRAM die area is 107 mm2 ... 16 Gb die densities".
// Boundary: whole DRAM die, which matches FEA's whole-die basis. Grade C/D.
inline double hbm_die_area_mm2() { return 107.0; }
inline double hbm_die_bits() { return 16e9; }
// HBM2 ACCESS energy, O'Connor et al. MICRO 2017:
// "with the average of 1.21 pJ/bit of activation energy, each HBM2 access
// incurs 3.92 pJ/bit". This is an ACCESS metric, NOT a stored-energy metric, so
// it is recorded but deliberately not used to fill M11's stored column.
inline double hbm_access_energy_J_per_bit() { return 3.97e-12; }

inline double hbm_column_energy_low_J_per_bit() { return 1.5e-12; }
inline double hbm_column_energy_high_J_per_bit() { return 5.7e-12; }

// ---- restoration model (declared, not measured) ----
// Shared by M1, M6 and M13 so the endpoint count cannot drift between modules.
inline double restoration_decay_um() { return 1.0; }
inline double restoration_input_threshold() { return 0.10; }
inline double restoration_endpoint_power_W() { return 1e-9; }
inline double restoration_endpoint_area_um2() { return 0.05; }

inline double restoration_max_spacing_m() {
    return restoration_decay_um() * 1e-6 * std::log(1.0 / restoration_input_threshold());
}

inline long long restoration_endpoint_count() {
    const double edge = std::sqrt(params().arch.die_area_cm2 * 1e-4);
    const double spacing = restoration_max_spacing_m();
    if (spacing <= 0.0) return 0;
    const int per_row = static_cast<int>(std::ceil(edge / spacing));
    return static_cast<long long>(per_row) * per_row;
}

inline double restoration_power_W() {
    return static_cast<double>(restoration_endpoint_count()) * restoration_endpoint_power_W();
}

inline double restoration_area_cm2() {
    return static_cast<double>(restoration_endpoint_count()) * restoration_endpoint_area_um2() *
           kCM2_PER_UM2;
}

} // namespace fea
