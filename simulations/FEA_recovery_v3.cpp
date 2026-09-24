// =============================================================================
// FEA_recovery_v3.cpp -- M17 peer recovery, quorum, and the tau gate
//
// Claim under test: a dead FZC can be rebuilt by its neighbours fast enough
// that the Zone's data is still there when it arrives, and the quorum rule
// suppresses false recoveries without making real ones unreachable.
//
// The gate this module exists to evaluate is arithmetic, not aesthetic:
//
//     T_detect + T_quorum + T_reconstruct  <<  tau
//
// A failed FZC stops L0 refresh with it, so the Zone's stored data decays at
// tau from the moment of failure. A recovery landing after tau installs a
// controller for data that is already gone, which is an empty Zone with extra
// steps. Everything reported here is subordinate to that inequality.
//
// Spec: spec/FZC-v0.md "Peer recovery: the Fusion Recovery Mesh",
//       "Declaring an FZC dead", invariants 9 to 12.
//
// Epistemic labels used below: ARITHMETIC means closed-form or direct
// substitution into a shared parameter; SWEPT means a parameter with no
// measurement behind it that is varied rather than assumed; SUGGESTED means a
// design recommendation this module produces. Nothing here is MEASURED, and no
// physical rescue port exists to exercise any of it.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace recovery {

using fea::params;
using fea::require;

// ---------------------------------------------------------------------------
// Closed-form binomial: P(X >= k) for X ~ Binom(n, p).
// Used instead of trusting a sample whenever the exact answer is available.
// ---------------------------------------------------------------------------
static double binom_tail(int n, int k, double p) {
    if (p <= 0.0) return k <= 0 ? 1.0 : 0.0;
    if (p >= 1.0) return k <= n ? 1.0 : 0.0;
    double term = std::pow(1.0 - p, n); // P(X = 0)
    double sum = (k <= 0) ? term : 0.0;
    for (int i = 1; i <= n; ++i) {
        term = term * (n - i + 1) / i * p / (1.0 - p);
        if (i >= k) sum += term;
    }
    return sum > 1.0 ? 1.0 : sum;
}

// Seeded sample of the same quantity, so the exact answer can be checked
// against an independent path rather than asserted by the code that derived it.
static double binom_tail_mc(int n, int k, double p, int trials, std::uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    int hits = 0;
    for (int t = 0; t < trials; ++t) {
        int count = 0;
        for (int i = 0; i < n; ++i) {
            if (u(rng) < p) ++count;
        }
        if (count >= k) ++hits;
    }
    return static_cast<double>(hits) / trials;
}

// ---------------------------------------------------------------------------
// The latency budget for one recovery, evaluated at a given heartbeat period.
// ---------------------------------------------------------------------------
struct Budget {
    double heartbeat_H_s = 0.0;
    double t_detect_s = 0.0;
    double t_quorum_s = 0.0;
    double t_reconstruct_s = 0.0;
    double total_s = 0.0;
    bool inside_tau = false;
};

static Budget budget_at(double H_s) {
    Budget b;
    b.heartbeat_H_s = H_s;
    b.t_detect_s = fea::heartbeat_miss_threshold_m() * H_s;
    b.t_quorum_s = fea::slingshot_message_latency_s();
    // Seed rides the fabric ceiling M10 establishes, at one packet's worth of bits.
    b.t_reconstruct_s =
        fea::fzc_seed_bits() / (fea::fabric_bandwidth_ceiling_GBps() * 1e9 * 8.0);
    b.total_s = b.t_detect_s + b.t_quorum_s + b.t_reconstruct_s;
    b.inside_tau = b.total_s < fea::retention_tau_s();
    return b;
}

// Largest heartbeat period that still satisfies the gate at a given latency.
static double max_viable_H_s(double slingshot_latency_s) {
    const double usable =
        fea::retention_tau_s() - slingshot_latency_s - fea::fzc_seed_bits() /
                                                          (fea::fabric_bandwidth_ceiling_GBps() *
                                                           1e9 * 8.0);
    return usable / fea::heartbeat_miss_threshold_m();
}

static std::string latency_str(double seconds) {
    std::ostringstream os;
    os << std::fixed;
    if (seconds < 1.0e-6) {
        os << std::setprecision(0) << seconds * 1e9 << " ns";
    } else if (seconds < 1.0e-3) {
        os << std::setprecision(0) << seconds * 1e6 << " us";
    } else {
        os << std::setprecision(0) << seconds * 1e3 << " ms";
    }
    return os.str();
}

// =============================================================================
// SCENARIO 1: the seed, and the ledger it comes from
// =============================================================================
static void scenario_seed_arithmetic() {
    std::cout << "\n[SCENARIO 1] what a neighbour actually transmits\n";

    const double seed = fea::fzc_seed_bits();
    const double ledger = fea::fzc_ledger_bits();

    std::cout << "  ledger: " << fea::fzc_state_groups() << " groups x "
              << fea::fzc_state_bits_per_group() << " bits x " << fea::fzc_state_rails()
              << " rails x " << fea::fzc_state_replicas() << " replicas\n";
    std::cout << "  unique state (one replica) : " << seed << " bits\n";
    std::cout << "  stored ledger (all replicas): " << ledger << " bits\n";
    std::cout << "  FZC Blocks per Zone        : " << fea::zone_fzc_blocks() << "\n";

    // PR6/S2: the gate that used to stand further down was
    // `seed * replicas == ledger`, an exact identity -- the header defines
    // fzc_ledger_bits() as fzc_seed_bits() * fzc_state_replicas() and `seed` is
    // fzc_seed_bits(), so both sides were the same expression. What the sentence
    // claims is that replication is MULTIPLICATIVE rather than additive, so both
    // readings are formed here and must be far apart. It is placed before the
    // spec literal so that changing the replica count trips this form test
    // rather than being caught downstream by `ledger == 336`.
    const double additive_ledger = seed + fea::fzc_state_replicas();
    require(std::fabs(ledger - additive_ledger) > 100.0,
            "the stored ledger must sit far from the additive reading "
            "(seed + replicas), or replication is being charged as a counter "
            "increment instead of a full copy");

    // The ledger must be the 336 that zone_fzc_blocks() documents as its state
    // component. If it is not, one of the two definitions has drifted.
    require(ledger == 336.0,
            "ledger must reconstruct the documented 336 state Blocks, or the two "
            "definitions of FZC state have diverged");
    // PR5/m3: `ledger == 336` alone only fails on a header edit, because both
    // sides are built from the same constants. This adds a structural property
    // the identity does not cover: controller state must be a strict SUBSET of
    // the FZC, so widening state past the whole FZC fails it.
    require(ledger < static_cast<double>(fea::zone_fzc_blocks()),
            "controller state must occupy strictly fewer Blocks than the whole FZC, "
            "or the state allocation has swallowed the command, port and spare Blocks");

    const double reconstruct_s =
        seed / (fea::fabric_bandwidth_ceiling_GBps() * 1e9 * 8.0);
    std::cout << "  reconstruction at the M10 ceiling: " << reconstruct_s * 1e9
              << " ns  (ARITHMETIC)\n";

    // Replication is insurance, not traffic: a neighbour sends one copy.
    require(seed < ledger, "the transmitted seed must be smaller than what is stored");

    const double packets = std::ceil(seed / fea::heartbeat_packet_bits());
    std::cout << "  as 64-bit packets         : " << packets << "\n";
    require(packets <= 2.0,
            "reconstruction must fit in a couple of packets, or recovery traffic is not "
            "the non-problem it was assumed to be");

    std::cout << "  label: ARITHMETIC from the FZC ledger. No measurement involved.\n";
}

// =============================================================================
// SCENARIO 2: the gate. How fast must the heartbeat beat?
// =============================================================================
static void scenario_the_gate() {
    std::cout << "\n[SCENARIO 2] the tau gate: recovery must beat retention decay\n";

    const double tau = fea::retention_tau_s();
    const double refresh = fea::refresh_interval_s();
    const double m = fea::heartbeat_miss_threshold_m();
    const double maxH = max_viable_H_s(fea::slingshot_message_latency_s());

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  retention tau               : " << tau * 1e3 << " ms\n";
    std::cout << "  refresh interval (tau/2)    : " << refresh * 1e3 << " ms\n";
    std::cout << "  miss threshold m            : " << m << "\n";
    std::cout << "  max viable heartbeat H      : " << maxH * 1e3 << " ms\n\n";

    const std::vector<double> candidates{
        1.0e-5, 1.0e-4, 5.0e-4, 1.0e-3, 5.0e-3, 1.0e-2, refresh, 5.0e-2,
    };

    std::cout << std::setw(14) << "H (ms)" << std::setw(14) << "detect (ms)"
              << std::setw(14) << "quorum (us)" << std::setw(14) << "recon (ns)"
              << std::setw(14) << "total (ms)" << std::setw(10) << "vs tau"
              << "\n";
    int passing = 0;
    for (const double H : candidates) {
        const Budget b = budget_at(H);
        if (b.inside_tau) ++passing;
        std::cout << std::setw(14) << H * 1e3 << std::setw(14) << b.t_detect_s * 1e3
                  << std::setw(14) << b.t_quorum_s * 1e6 << std::setw(14)
                  << b.t_reconstruct_s * 1e9 << std::setw(14) << b.total_s * 1e3
                  << std::setw(10) << (b.inside_tau ? "PASS" : "FAIL") << "\n";
    }

    require(passing > 0 && passing < static_cast<int>(candidates.size()),
            "the gate must separate candidates: if everything passes or everything fails, "
            "the constraint is not binding and this sweep proves nothing");

    // The refresh interval itself must FAIL the gate. If it passed, the derived
    // max H would be >= tau/2 and the heartbeat would be free, which would make
    // this whole module pointless.
    const Budget at_refresh = budget_at(refresh);
    require(!at_refresh.inside_tau,
            "the refresh interval must NOT satisfy the recovery gate, or detection is "
            "not the binding constraint");

    require(maxH > 0.0 && maxH < tau,
            "max viable H must be positive and strictly below tau");

    std::cout << "\n  ARITHMETIC: the binding constraint is m x H, not the packet path.\n";
    std::cout << "  At m = " << m << " the heartbeat must run at least "
              << (refresh / maxH) << "x faster than the refresh interval.\n";
    std::cout << "  SUGGESTED: H = 10 ms gives " << budget_at(1.0e-2).total_s * 1e3
              << " ms total against a " << tau * 1e3 << " ms budget.\n";
    std::cout << "  label: ARITHMETIC. H and m are declared, not measured.\n";
}

// =============================================================================
// SCENARIO 3: does the unsourced latency get to decide the answer?
// =============================================================================
static void scenario_latency_independence() {
    std::cout << "\n[SCENARIO 3] the declared, unsourced Slingshot latency\n";

    require(!fea::slingshot_latency_sourced(),
            "this scenario exists only while the latency is unsourced; if it becomes "
            "sourced the sweep should be retired rather than left asserting the obvious");

    const double reference = max_viable_H_s(1.0e-9);
    const std::vector<double> sweep{1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3};

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  " << std::left << std::setw(16) << "latency" << std::right
              << std::setw(16) << "max viable H (ms)" << std::setw(14) << "shift vs 1ns"
              << "\n";
    double worst_shift = 0.0;
    for (const double L : sweep) {
        const double H = max_viable_H_s(L);
        const double shift = std::fabs(H - reference) / reference;
        if (shift > worst_shift) worst_shift = shift;
        std::cout << "  " << std::left << std::setw(16) << latency_str(L) << std::right
                  << std::setw(16) << H * 1e3 << std::setw(13)
                  << (shift * 100.0) << "%\n";
    }

    std::cout << "\n  worst shift across six decades: " << worst_shift * 100.0 << "%\n";
    require(worst_shift < 0.05,
            "the gate must move by less than 5% across six decades of latency, otherwise an "
            "unsourced input is deciding a headline result");

    std::cout << "  conclusion: m x H carries this result, the packet path does not.\n";
    std::cout << "  label: SWEPT over an unsourced input, deliberately.\n";
}

// =============================================================================
// SCENARIO 4: what the quorum buys, and what it costs
// =============================================================================
static void scenario_quorum_tradeoff() {
    std::cout << "\n[SCENARIO 4] quorum 4-of-5 against quorum 3-of-5\n";

    const int n = static_cast<int>(fea::recovery_peer_count());
    const int k4 = static_cast<int>(fea::recovery_quorum());
    const int k3 = 3;

    std::cout << "  peer set N = " << n << ", quorum ceil(2/3 x N) = " << k4 << "\n\n";

    const std::vector<double> congestion{0.0, 0.05, 0.10, 0.20, 0.30, 0.50};
    const int trials = 200000;
    std::uint32_t seed = 20260923u;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "  " << std::left << std::setw(12) << "congest" << std::setw(14)
              << "FP k=4" << std::setw(14) << "FP k=3" << std::setw(14)
              << "detect k=4" << std::setw(14) << "falls to L2" << std::setw(14)
              << "MC err" << "\n";

    for (const double p : congestion) {
        // Zone ALIVE: congestion swallows heartbeats, peers wrongly report death.
        const double fp4 = binom_tail(n, k4, p);
        const double fp3 = binom_tail(n, k3, p);
        // Zone DEAD: only peers not themselves congested can report.
        const double det4 = binom_tail(n, k4, 1.0 - p);

        // Independent sampled path, to check the closed form against.
        const double mc = binom_tail_mc(n, k4, p, trials, seed++);
        const double err = std::fabs(mc - fp4);

        std::cout << "  " << std::left << std::setw(12) << p << std::setw(14) << fp4
                  << std::setw(14) << fp3 << std::setw(14) << det4 << std::setw(14)
                  << (1.0 - det4) << std::setw(14) << err << "\n";

        // The closed form and the sample must agree. 5 sigma of a proportion on
        // 2e5 trials is comfortably wider than any drift this should show.
        const double tolerance = 5.0 * std::sqrt(std::max(fp4, 1e-9) *
                                                 (1.0 - std::max(fp4, 1e-9)) / trials);
        require(err <= tolerance + 1e-4,
                "sampled and closed-form quorum probabilities must agree, or one of the "
                "two models of the quorum is wrong");

        if (p > 0.0 && p < 1.0) {
            require(fp4 < fp3,
                    "raising the quorum from 3 to 4 must strictly reduce false recoveries, "
                    "or the extra agreement buys nothing");
            require(det4 > 0.0 && det4 < 1.0,
                    "at partial congestion the quorum must be met sometimes and missed "
                    "sometimes, or the recovery ladder is decorative");
        }
        // PR5/S3: the gate that stood here tested det4 against [0,1], which
        // binom_tail guarantees by construction (it clamps to [0,1] and every
        // term is non-negative), so it could never fail. The boundary behaviour
        // of the model itself is falsifiable: if the clamp or the endpoint
        // handling broke, these two would not hold.
        require(binom_tail(n, k4, 1.0) == 1.0 && binom_tail(n, k4, 0.0) == 0.0,
                "detection must be certain with no congestion and impossible under total "
                "congestion, or the quorum model is broken at its endpoints");
    }

    const double fp4_10 = binom_tail(n, k4, 0.10);
    const double fp3_10 = binom_tail(n, k3, 0.10);
    std::cout << "\n  at 10% congestion: k=4 false recovery " << fp4_10 << ", k=3 " << fp3_10
              << " -> " << (fp3_10 / fp4_10) << "x fewer with the stricter quorum\n";
    const double det4_30 = binom_tail(n, k4, 0.70);
    std::cout << "  at 30% congestion, L1 reaches quorum " << det4_30 * 100.0
              << "% of the time; the rest falls to L2\n";
    std::cout << "  label: ARITHMETIC (exact binomial), checked against a seeded sample.\n";
    std::cout << "  The congestion RATE is SWEPT, not measured.\n";
}

// =============================================================================
// SCENARIO 5: the heartbeat's own cost across the whole die
// =============================================================================
static void scenario_heartbeat_bandwidth() {
    std::cout << "\n[SCENARIO 5] what the heartbeat costs across every Zone\n";

    const double zones = fea::design_zone_count();
    const double packet = fea::heartbeat_packet_bits();
    const double H = 1.0e-2;

    // Every Zone ticks every H, so this is zones x packet / H in bits.
    const double heartbeat_bits_per_s = zones * packet / H;
    const double heartbeat_GBps = heartbeat_bits_per_s / 8.0 / 1e9;

    // The refresh traffic the design already pays, from the same store.
    const double refresh_GBps =
        fea::refresh_bits() / fea::refresh_interval_s() / 8.0 / 1e9;
    const double share = heartbeat_GBps / refresh_GBps;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Zones on the die        : " << std::setprecision(0) << zones << "\n";
    std::cout << "  heartbeat packet        : " << std::setprecision(0) << packet << " bits\n";
    std::cout << "  H                       : " << std::setprecision(3) << H * 1e3 << " ms\n";
    std::cout << "  heartbeat traffic       : " << std::setprecision(3) << heartbeat_GBps
              << " GB/s\n";
    std::cout << "  refresh traffic         : " << refresh_GBps << " GB/s\n";
    std::cout << "  heartbeat share         : " << std::setprecision(4) << share * 100.0
              << " %\n";

    require(share < 0.05,
            "the heartbeat must stay under 5% of refresh traffic, or a monitoring signal "
            "costs more than the thing it monitors");
    require(heartbeat_GBps < fea::fabric_bandwidth_ceiling_GBps(),
            "heartbeat traffic alone must fit under M10's fabric ceiling");

    // Faster detection is bought with bandwidth, so show what halving H costs.
    const double fast_H = H / 2.0;
    const double fast_share = (zones * packet / fast_H) / 8.0 / 1e9 / refresh_GBps;
    std::cout << "  at H = " << fast_H * 1e3 << " ms the share becomes "
              << fast_share * 100.0 << " %\n";
    std::cout << "  label: ARITHMETIC from design_zone_count() and the declared packet size.\n";
    std::cout << "  The packet size is DECLARED framing, not a packet specification.\n";
}

// =============================================================================
// SCENARIO 6: invariant 12, the epoch fence
// =============================================================================
static void scenario_epoch_fence() {
    std::cout << "\n[SCENARIO 6] a partitioned Zone must not run two controllers\n";

    struct Commit {
        long long epoch;
        bool from_recovered;
        bool accepted;
    };

    const long long original_epoch = 7;
    // Peers cannot hear the Zone, reach quorum, and rebuild it. The rebuild
    // bumps the epoch as invariant 12 requires.
    const long long recovery_epoch = original_epoch + 1;

    // PR5/m2: acceptance was typed in as literal `true` flags, so the gate
    // counted a vector the author had already decided about. It is now derived
    // from the rule invariant 12 states: a commit is accepted only if its epoch
    // is the current one, and without a fence the current epoch never moves.
    auto accepted = [](long long commit_epoch, long long current_epoch) {
        return commit_epoch == current_epoch;
    };
    const long long unfenced_current = original_epoch;  // no bump, nothing fences
    const long long fenced_current = recovery_epoch;    // invariant 12 bumps it
    const std::vector<Commit> unfenced{
        {original_epoch, false, accepted(original_epoch, unfenced_current)},
        {original_epoch, true, accepted(original_epoch, unfenced_current)},
    };
    const std::vector<Commit> fenced{
        {original_epoch, false, accepted(original_epoch, fenced_current)},
        {recovery_epoch, true, accepted(recovery_epoch, fenced_current)},
    };

    int live_unfenced = 0;
    for (const Commit& c : unfenced) {
        if (c.accepted) ++live_unfenced;
    }
    int live_fenced = 0;
    for (const Commit& c : fenced) {
        if (c.accepted) ++live_fenced;
    }

    std::cout << "  original epoch " << original_epoch << ", recovery epoch "
              << recovery_epoch << "\n";
    std::cout << "  without a fence: " << live_unfenced
              << " controllers accepted commits -> split brain\n";
    std::cout << "  with the fence : " << live_fenced << " controller accepted commits\n";

    require(live_unfenced == 2,
            "the unfenced case must actually demonstrate two live controllers, or the "
            "hazard this invariant prevents has not been shown");
    require(live_fenced == 1,
            "invariant 12 must leave exactly one controller accepting commits");

    // The rejected commits are precisely the ones from the lower epoch.
    for (const Commit& c : fenced) {
        require(c.accepted == (c.epoch == recovery_epoch),
                "a commit may be accepted only if it carries the current epoch");
    }

    std::cout << "  invariant 12 holds: lower-epoch commits are rejected.\n";
    std::cout << "  label: ARITHMETIC demonstration of the rule, not a timed simulation.\n";
}

// =============================================================================
// SCENARIO 7: what this module deliberately does not claim
// =============================================================================
static void scenario_what_this_does_not_show() {
    std::cout << "\n[SCENARIO 7] what M17 deliberately does not claim\n";

    const double tau = fea::retention_tau_s();
    const double budget = budget_at(1.0e-2).total_s;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  recovery at H = 10 ms       : " << budget * 1e3 << " ms against tau "
              << tau * 1e3 << " ms\n";
    require(budget < tau, "the worked example must sit inside tau to be quotable at all");

    std::cout << "\n  none of the following is established here:\n";
    std::cout << "  - no rescue port exists, so nothing above has ever been exercised on\n";
    std::cout << "    hardware; the mechanism is specified, not demonstrated.\n";
    std::cout << "  - H and m are declared. Neither is measured, and FZC-v0 lists both as\n";
    std::cout << "    open parameters.\n";
    std::cout << "  - the Slingshot latency is unsourced; scenario 3 exists to show it does\n";
    std::cout << "    not carry the result rather than to endorse a value for it.\n";
    std::cout << "  - congestion rates are SWEPT. No chip traffic has produced one.\n";
    std::cout << "  - scenario 6 is a rule demonstration, not a timed partition, so it shows\n";
    std::cout << "    the fence rejects stale epochs, not how fast it does so.\n";
    std::cout << "  - heartbeat bandwidth across every Zone IS reported, in scenario 5,\n";
    std::cout << "    but from a declared packet size rather than from a real packet format.\n";
    std::cout << "  label: capability PROPOSED, latency OPEN, mechanism NOT DEMONSTRATED.\n";
    std::cout << "  NEXT EVIDENCE GATE: a physical rescue port, and a timed partition with\n";
    std::cout << "  measured congestion from a real Slingshot load.\n";
}

} // namespace recovery

int main() {
    using namespace recovery;
    try {
        std::cout << "FEA V3 M17 peer recovery: quorum, seed, and the tau gate\n";
        std::cout << "Same 0.5 cm^2 design array, same Kramers retention as M4 and M13.\n";
        scenario_seed_arithmetic();
        scenario_the_gate();
        scenario_latency_independence();
        scenario_quorum_tradeoff();
        scenario_heartbeat_bandwidth();
        scenario_epoch_fence();
        scenario_what_this_does_not_show();
        std::cout << "\nPASS: recovery beats retention decay and the quorum suppresses false "
                     "recoveries.\n";
        std::cout << "LABEL: arithmetic and swept inputs only, mechanism not demonstrated.\n";
        std::cout << "NEXT EVIDENCE GATE: physical rescue port and a timed partition under "
                     "measured congestion.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
