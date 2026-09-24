// =============================================================================
// FEA_bandwidth_v3.cpp -- M10 per-Zone rate to chip aggregate
//
// Reviewer 6 point 5: 133 GOPS/Zone for 64-bit ops already equals about
// 1.06 TB/s, nearly identical to the reported chip-level value, so the
// aggregation across Zones must be clarified. This module reproduces V2's
// per-Zone figures, converts them to bytes per second, compares that against
// V2's stated chip aggregate, and then applies the Zone multiplier that V2
// appears to have omitted. A previous revision of this module also charged M13's
// refresh traffic against that aggregate; that comparison was wrong and is
// withdrawn in SCENARIO 4, because refresh writes stay inside a Zone.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

namespace bandwidth {

using fea::params;
using fea::require;
using fea::p_abs_single;
using fea::zone_count_stated;

static constexpr double kBytesPerWord = 8.0; // 64-bit operations

// V2's own SIM 12 contention-model figures.
static double gops_row_broadcast() { return 147.0e9; }
static double gops_random() { return 133.0e9; }
static double gops_strided() { return 9.0e9; }

// V2's stated aggregate, in bytes per second.
static double v2_chip_Bps_random() { return 1.1e12; }   // "~1.1 TB/s"
static double v2_chip_Bps_strided() { return 75.0e9; }  // "~75 GB/s"

// Zone count. V2's headline uses 1.7e9; its bandwidth prose says ~1e9.
static double zone_count_v2_prose() { return 1.0e9; }
static double zone_count_headline() { return zone_count_stated(); }

static void scenario_reproduce_per_zone() {
    std::cout << "\n[SCENARIO 1] reproduce V2's per-Zone contention figures\n";
    std::cout << std::fixed << std::setprecision(0);
    std::cout << "  row-broadcast / sequential : " << (gops_row_broadcast() / 1e9) << " GOPS/Zone\n";
    std::cout << "  random access              : " << (gops_random() / 1e9) << " GOPS/Zone\n";
    std::cout << "  worst-case strided         : " << (gops_strided() / 1e9) << " GOPS/Zone\n";
    require(gops_random() > 0.0 && gops_strided() < gops_random(),
            "strided access must be slower than random access");
    const double random_Bps = gops_random() * kBytesPerWord;
    const double strided_Bps = gops_strided() * kBytesPerWord;
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  random -> bytes/s per Zone : " << random_Bps << " B/s ("
              << std::fixed << std::setprecision(3) << (random_Bps / 1e12) << " TB/s)\n";
    std::cout << "  strided -> bytes/s per Zone: " << std::scientific << strided_Bps
              << " B/s (" << std::fixed << std::setprecision(3) << (strided_Bps / 1e9)
              << " GB/s)\n";
    std::cout << "  these are PER ZONE figures. They have not been aggregated yet.\n";
}

static void scenario_v2_chip_equals_one_zone() {
    std::cout << "\n[SCENARIO 2] V2's chip aggregate is numerically one Zone (Reviewer 6 #5)\n";
    const double per_zone_random = gops_random() * kBytesPerWord;
    const double per_zone_strided = gops_strided() * kBytesPerWord;
    const double v2_random = v2_chip_Bps_random();
    const double v2_strided = v2_chip_Bps_strided();

    const double ratio_random = v2_random / per_zone_random;
    const double ratio_strided = v2_strided / per_zone_strided;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  random:  one Zone = " << per_zone_random / 1e12 << " TB/s, V2 chip = "
              << v2_random / 1e12 << " TB/s, ratio " << ratio_random << "x\n";
    std::cout << "  strided: one Zone = " << per_zone_strided / 1e9 << " GB/s, V2 chip = "
              << v2_strided / 1e9 << " GB/s, ratio " << ratio_strided << "x\n\n";

    require(std::abs(ratio_random - 1.0) < 0.15,
            "V2's chip random bandwidth must equal about one Zone's bandwidth");
    require(std::abs(ratio_strided - 1.0) < 0.15,
            "V2's chip strided bandwidth must equal about one Zone's bandwidth");
    std::cout << "  both V2 chip figures land within 15% of a SINGLE Zone. The multiplication\n";
    std::cout << "  by 10^9 Zones was never performed. Reviewer 6 #5 is exactly right.\n";
    std::cout << "  so V2 has two internally consistent but mutually exclusive readings:\n";
    std::cout << "    (a) the chip really delivers 1.1 TB/s, and 10^9 Zones are unreachable, or\n";
    std::cout << "    (b) every Zone delivers, and the aggregate is 10^9 times larger.\n";
}

static void scenario_true_aggregate() {
    std::cout << "\n[SCENARIO 3] applying the Zone multiplier V2 omitted\n";
    const double per_zone = gops_random() * kBytesPerWord;
    const double aggregate_prose = per_zone * zone_count_v2_prose();
    const double aggregate_headline = per_zone * zone_count_headline();
    const double v2_claim = v2_chip_Bps_random();

    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  per-Zone random                 : " << per_zone << " B/s\n";
    std::cout << "  x 1e9 Zones (V2 prose)          : " << aggregate_prose << " B/s ("
              << std::fixed << std::setprecision(0) << (aggregate_prose / 1e18)
              << " EB/s)\n";
    std::cout << "  x 1.7e9 Zones (V2 headline)     : " << std::scientific << aggregate_headline
              << " B/s\n";
    std::cout << "  V2 stated chip aggregate        : " << v2_claim << " B/s\n";
    std::cout << "  understatement factor            : " << std::fixed << std::setprecision(1)
              << (aggregate_prose / v2_claim) << "x\n\n";

    require(aggregate_prose > 1.0e6 * v2_claim,
            "the true aggregate must exceed V2's stated chip figure by many orders of magnitude");
    require(zone_count_v2_prose() > 1.0e8, "the Zone multiplier must be a large number");

    std::cout << "  a true 10^9-Zone aggregate of " << std::fixed << std::setprecision(0)
              << (aggregate_prose / 1e18) << " EB/s cannot be routed on any die.\n";
    std::cout << "  so reading (b) is physically impossible, which forces reading (a):\n";
    std::cout << "  V2's usable bandwidth is that of one Zone, and the other 10^9 Zones do not\n";
    std::cout << "  contribute. Either way the 1.1 TB/s figure does not describe a 10^9-Zone chip.\n";
}

// CORRECTION. An earlier revision of this scenario compared refresh traffic
// against V2's advertised port and called it a 737x overrun. That comparison is
// wrong: a refresh write is a local FIRE onto a Block inside its own Word, so it
// never crosses an inter-Zone link or the chip port. Comparing internal state
// maintenance with external bandwidth conflates two different things. The honest
// figure is a local duty cycle, which turns out to be small. Refresh is an
// energy and actuator-throughput problem, not a bandwidth problem.
static void scenario_refresh_is_local() {
    std::cout << "\n[SCENARIO 4] refresh is LOCAL, so it is not port bandwidth\n";
    const auto& c = params().control;
    const auto& d = params().device;
    const auto& a = params().arch;

    // A Zone holds 1024 Words by definition, whatever the Zone count is.
    // payload/zones was used here before, which mixes our design payload with
    // V2's stated Zone count and silently under-reports the per-Zone size.
    const double bits_per_zone = fea::zone_data_blocks();
    const double bits_per_word = 64.0; // a Word is 64 Blocks, one bit each
    const double fires_per_zone = bits_per_zone / bits_per_word;

    // Expected FIRE latency, derived the same way M9 derives it.
    const double vg = 2.0 * d.t_hop_eV * fea::kEV * d.a_lattice_m / fea::kHbar;
    const double t_fire_ps = (a.segment_um * 1e-6) / vg / fea::kPS;
    const double expected_fire_ps = t_fire_ps / p_abs_single();

    const double interval = fea::refresh_interval_s();
    const double local_time_s = fires_per_zone * expected_fire_ps * fea::kPS;
    const double duty = local_time_s / interval;
    const double refresh_Bps = fea::refresh_bits() / 8.0 / interval;
    const double bit_writes_per_s = fea::refresh_bits() / interval;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  payload per Zone             : " << std::setprecision(0) << bits_per_zone
              << " bits = " << c.words_per_zone << " Words\n";
    std::cout << "  FIREs to refresh one Zone    : " << std::setprecision(0) << fires_per_zone
              << " (one per Word, 64 bits in parallel)\n";
    std::cout << "  expected FIRE latency (M9)   : " << std::setprecision(3) << expected_fire_ps
              << " ps\n";
    std::cout << "  sequential full-Zone refresh : " << std::setprecision(1)
              << (local_time_s * 1e9) << " ns\n";
    std::cout << "  refresh interval             : " << std::setprecision(4) << (interval * 1e3)
              << " ms\n";
    std::cout << "  LOCAL duty cycle             : " << std::scientific << std::setprecision(3)
              << duty << "  (" << std::fixed << std::setprecision(4) << (duty * 100.0) << " %)\n";
    std::cout << "  V2 claimed refresh overhead  : 1.1e-4  (0.011 %)\n";
    std::cout << "  V2 is " << std::fixed << std::setprecision(1) << (1.1e-4 / duty)
              << "x MORE conservative than this local estimate.\n\n";

    require(local_time_s < interval, "a full local Zone refresh must fit inside one interval");
    require(duty < 1.1e-4,
            "local refresh duty cycle must sit below V2's own stated 1.1e-4 overhead");
    require(fires_per_zone > 1000.0, "a Zone needs on the order of 1000 local FIREs to refresh");

    std::cout << "  WITHDRAWN: an earlier revision of this module called refresh "
              << std::setprecision(0) << (refresh_Bps / 1e12) << " TB/s against a "
              << (v2_chip_Bps_random() / 1e12) << " TB/s port, a "
              << (refresh_Bps / v2_chip_Bps_random()) << "x overrun.\n";
    std::cout << "  That mixed internal state maintenance with external bandwidth. It is wrong\n";
    std::cout << "  and must not be quoted in the manuscript or the point-by-point response.\n\n";

    std::cout << "  what refresh really costs is ACTIVITY and ENERGY, not port bandwidth:\n";
    std::cout << "    bit-writes per second       : " << std::scientific << std::setprecision(3)
              << bit_writes_per_s << " /s\n";
    std::cout << "    at declared 1e-18 J/bit     : " << (bit_writes_per_s * 1e-18) << " W\n";
    std::cout << "    at 1e-16 J/bit              : " << (bit_writes_per_s * 1e-16) << " W\n";
    std::cout << "  every Zone refreshes in parallel, so no single path carries the sum. This is\n";
    std::cout << "  exactly the local scheduling FZC already owns in FZC-v0.\n";
    std::cout << "  note: FZC's REFRESH opcode scrubs protected state from redundant copies,\n";
    std::cout << "  which is ECC correction, not retention refresh. Two different operations.\n";
    std::cout << "  label: duty derived from M9 latency and V2 Kramers cadence, energy per bit DECLARED.\n";
}

static void scenario_serviceable_zones() {
    std::cout << "\n[SCENARIO 5] how many Zones the claimed bandwidth can actually serve\n";
    const double per_zone = gops_random() * kBytesPerWord;
    const double claimed = v2_chip_Bps_random();
    const double serviceable = claimed / per_zone;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  claimed chip bandwidth          : " << claimed / 1e12 << " TB/s\n";
    std::cout << "  per-Zone demand                 : " << per_zone / 1e12 << " TB/s\n";
    std::cout << "  Zones simultaneously serviceable: " << serviceable << "\n";
    std::cout << "  Zones on the die                : " << std::setprecision(0)
              << (zone_count_headline() / 1e6) << " million\n";
    const double fraction = serviceable / zone_count_headline();
    std::cout << "  fraction reachable              : " << std::scientific << std::setprecision(3)
              << fraction << "  (" << fraction * 100.0 << " %)\n\n";

    require(serviceable < 2.0, "the claimed bandwidth must service about one Zone");
    require(fraction < 1.0e-6, "the reachable fraction of the die must be vanishingly small");
    std::cout << std::fixed << std::setprecision(0);
    std::cout << "  V2's bandwidth serves about " << serviceable
              << " Zone(s) of " << (zone_count_headline() / 1e6) << " million.\n";
    std::cout << "  finding: a 14.1 TB array behind a 1.1 TB/s port behaves like memory with a\n";
    std::cout << "  single-tenant crossbar. Capacity does not create bandwidth. Any bandwidth\n";
    std::cout << "  claim must state how many Zones are concurrently reachable through shared\n";
    std::cout << "  routing, which is the physical link model V3 still does not have.\n";
    std::cout << "  label: per-Zone figures from V2, aggregation arithmetic derived, link model open.\n";
}

} // namespace bandwidth

int main() {
    using namespace bandwidth;
    try {
        std::cout << "FEA V3 M10 bandwidth aggregation\n";
        std::cout << "All figures derived from V2's own per-Zone contention numbers.\n";
        scenario_reproduce_per_zone();
        scenario_v2_chip_equals_one_zone();
        scenario_true_aggregate();
        scenario_refresh_is_local();
        scenario_serviceable_zones();
        std::cout << "\nPASS: V2's chip bandwidth identified as one Zone, multiplier omitted, refresh shown local.\n";
        std::cout << "LABEL: derived aggregation, V2 per-Zone inputs, physical link model open.\n";
        std::cout << "NEXT EVIDENCE GATE: a shared-routing link capacity model so concurrent Zone count is derived, not assumed.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
