// =============================================================================
// FEA_fzc_opensystem_v3.cpp -- open-system capture and restoration accounting
//
// Scope: named probability channels for one FZC selector stage. It checks that
// reflection, transmission, temporary actuator occupation, stored occupation,
// and reservoir exchange are separately accounted, and that a clock/bias energy
// source is required before any restoration or cascade is allowed.
//
// These are accounting gates, not calibrated Si dangling-bond rates, energies,
// times, or device-physics predictions.
// =============================================================================

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace opensystem {

struct Channels {
    double reflected = 0.0;
    double transmitted = 0.0;
    double actuator = 0.0;
    double stored = 0.0;
    double reservoir = 0.0;

    double sum() const { return reflected + transmitted + actuator + stored + reservoir; }
};

struct Stage {
    Channels channels{};
    double clock_energy = 0.0;
    double output_margin = 0.0;
    // NOTE: no label or intent flag. A previous revision carried a
    // `labeled_capture` boolean that the harness set by hand, which let the
    // capture claim pass because the test asked for it rather than because the
    // channels earned it. Capture is a pure function of the channels below.
};

static void require(bool ok, const std::string& what) {
    if (!ok) throw std::runtime_error("ASSERTION FAILED: " + what);
}

static bool conserved(const Channels& c, double tol = 1e-12) {
    return std::abs(c.sum() - 1.0) <= tol;
}

static bool capture_allowed(const Stage& stage) {
    // Capture needs a retained stored channel AND an explicit reservoir that
    // received the relaxation energy, and the five channels must conserve.
    // No caller flag decides this: the channels alone do.
    return stage.channels.stored > 0.0 && stage.channels.reservoir > 0.0 &&
           conserved(stage.channels);
}

static bool restoration_allowed(const Stage& stage, double required_margin) {
    return capture_allowed(stage) && stage.clock_energy > 0.0 && stage.output_margin >= required_margin;
}

static void conservation_gate() {
    std::cout << "\n[SCENARIO 1] named channels must conserve probability\n";
    const Channels ok{0.20, 0.30, 0.10, 0.25, 0.15};
    const Channels leak{0.20, 0.30, 0.10, 0.25, 0.00};
    require(conserved(ok), "explicit five-channel model must sum to 1");
    require(!conserved(leak), "missing reservoir must fail conservation when probability was removed");
    std::cout << "  five named channels sum to 1; dropped reservoir fails\n";
}

static void transmission_is_not_capture() {
    std::cout << "\n[SCENARIO 2] transmission or unlabeled loss is not capture\n";
    Stage passed_through;
    passed_through.channels = Channels{0.40, 0.60, 0.0, 0.0, 0.0};
    Stage norm_loss;
    norm_loss.channels = Channels{0.20, 0.50, 0.0, 0.30, 0.0};
    require(conserved(passed_through.channels), "pure transport must still conserve probability");
    require(!capture_allowed(passed_through), "transmission without stored occupation is not capture");
    require(conserved(norm_loss.channels), "a stored label without reservoir exchange can still sum to 1");
    require(!capture_allowed(norm_loss), "stored label without reservoir exchange is not capture");
    std::cout << "  transport-only and reservoir-free loss both rejected as capture\n";
}

static void explicit_relaxation_accounts_capture() {
    std::cout << "\n[SCENARIO 3] explicit relaxation can pass the accounting gate only\n";
    Stage relaxed;
    relaxed.channels = Channels{0.15, 0.10, 0.05, 0.40, 0.30};
    require(capture_allowed(relaxed), "stored occupation plus reservoir exchange passes the accounting gate");
    require(!restoration_allowed(relaxed, 1.0), "capture accounting without clock energy is not restoration");
    std::cout << "  stored=0.40 reservoir=0.30 passes capture accounting; clock energy still missing\n";
}

static void clocked_cascade() {
    std::cout << "\n[SCENARIO 4] cascade requires clock energy and downstream margin\n";
    Stage a;
    a.channels = Channels{0.10, 0.05, 0.05, 0.50, 0.30};
    a.clock_energy = 1.0;
    a.output_margin = 1.5;
    const double downstream_required = 1.2;
    require(restoration_allowed(a, downstream_required), "stage A must clear the downstream input requirement");
    Stage b;
    b.channels = Channels{0.05, 0.05, 0.05, 0.55, 0.30};
    b.clock_energy = 1.0;
    b.output_margin = 0.4;
    require(restoration_allowed(b, 0.2), "stage B must pass its own lower restoration gate");
    require(!restoration_allowed(b, downstream_required), "stage B output must fail the tighter downstream requirement");

    Stage unpowered = a;
    unpowered.clock_energy = 0.0;
    require(!restoration_allowed(unpowered, downstream_required), "removing clock energy must block cascade");
    std::cout << "  clock energy=1.0 and output margin=1.5 allow a two-stage accounting cascade\n";
    std::cout << "  clock energy=0 blocks the same cascade\n";
}

static void channels_alone_decide_capture() {
    std::cout << "\n[SCENARIO 5] capture is a pure function of the channels, not of intent\n";
    // Truth table. Each row breaks exactly one of the three conditions capture
    // requires, so each row must fail for its own reason rather than because a
    // flag was left unset.
    struct Row { const char* name; Channels c; bool expect_capture; };
    const Row rows[] = {
        {"all conditions met",        {0.15, 0.10, 0.05, 0.40, 0.30}, true},
        {"no stored occupation",      {0.40, 0.60, 0.00, 0.00, 0.00}, false},
        {"stored but no reservoir",   {0.20, 0.50, 0.00, 0.30, 0.00}, false},
        {"channels do not conserve",  {0.20, 0.50, 0.00, 0.50, 0.30}, false},
        {"reservoir but no stored",   {0.40, 0.30, 0.00, 0.00, 0.30}, false},
    };
    // PR5/m6: the `matched` counter was removed with its gate; the per-row
    // require above throws on any mismatch, so counting matches was dead.
    const int row_count = static_cast<int>(sizeof(rows) / sizeof(rows[0]));
    for (int i = 0; i < row_count; ++i) {
        Stage s;
        s.channels = rows[i].c;
        const bool got = capture_allowed(s);
        const bool cons = conserved(s.channels);
        std::cout << "  " << std::left << std::setw(26) << rows[i].name << std::right
                  << " conserved=" << (cons ? "yes" : "no ")
                  << " stored=" << std::fixed << std::setprecision(2) << s.channels.stored
                  << " reservoir=" << s.channels.reservoir
                  << " -> capture=" << (got ? "YES" : "no ")
                  << (got == rows[i].expect_capture ? "   ok" : "   WRONG") << "\n";
        require(got == rows[i].expect_capture,
                "capture verdict must follow the channels, not a label");
    }
    // PR5/m6: the summary gate that stood here (matched == row_count) could only
    // fail if the per-row require above had already thrown for the same row. It
    // is replaced by a coverage property the loop never checks: the truth table
    // must exercise both outcomes, or matching every row proves nothing.
    bool saw_capture = false, saw_transmit = false;
    for (const auto& row : rows) {
        if (row.expect_capture) saw_capture = true;
        else saw_transmit = true;
    }
    require(saw_capture && saw_transmit,
            "the channel truth table must contain both capture and no-capture rows, "
            "or matching all of them is vacuous");
    // Independence, tested by MUTATION rather than self-comparison. The gate
    // that used to sit here built two Stages with IDENTICAL channels and asserted
    // f(x) == f(x), which a pure function always satisfies. Now each capture
    // condition is broken in turn and the verdict must flip, plus one mutation
    // that must NOT flip. Both directions have to hold, so neither direction
    // can pass for free.
    Stage a; a.channels = Channels{0.15, 0.10, 0.05, 0.40, 0.30};
    const bool before = capture_allowed(a);
    require(before, "the full-condition Stage must be capture");
    int flips = 0;
    Stage b = a; b.channels.reservoir = 0.0;   // drop a required condition
    if (capture_allowed(b) != before) ++flips;
    Stage c = a; c.channels.stored = 0.0;      // drop the other required condition
    if (capture_allowed(c) != before) ++flips;
    Stage d = a; d.channels.stored += 0.2;     // break conservation
    if (capture_allowed(d) != before) ++flips;
    // Must NOT flip: redistributes the non-capture channels while keeping
    // stored, reservoir and the sum of 1.0 all intact.
    Stage e = a; e.channels.reflected = 0.20; e.channels.transmitted = 0.05;
    const bool held = (capture_allowed(e) == before);
    require(flips == 3, "breaking each capture condition in turn must flip the verdict");
    require(held, "redistributing non-capture channels must NOT flip the verdict");
    std::cout << "  3 of 3 condition mutations flipped the verdict; 1 non-condition "
              << "mutation did not.\n";
    std::cout << "  no intent flag remains: the channels alone decide capture.\n";
    std::cout << "  label: capture decided by stored + reservoir + conservation only.\n";
}

} // namespace opensystem

int main() {
    using namespace opensystem;
    try {
        std::cout << "FEA V3 FZC open-system accounting simulator\n";
        std::cout << "Channels and clock energy are accounting parameters, not device predictions.\n";
        conservation_gate();
        transmission_is_not_capture();
        explicit_relaxation_accounts_capture();
        clocked_cascade();
        channels_alone_decide_capture();
        std::cout << "\nPASS: capture and restoration accounting gates held by channels alone.\n";
        std::cout << "NEXT EVIDENCE GATE: replace these probabilities with a stated geometry, Hamiltonian or rate model, clock waveform, and sensor path.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
