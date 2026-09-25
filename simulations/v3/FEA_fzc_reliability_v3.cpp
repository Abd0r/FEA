// =============================================================================
// FEA_fzc_reliability_v3.cpp -- protected FZC state and recovery simulator
//
// Scope: logical controller protection semantics with injected abstract faults.
// It does not assign a physical bit error rate, retention time, ECC area, energy,
// or capture mechanism. Those values must come from lower-level device models.
// =============================================================================

#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace reliability {

enum class Phase { Idle, Armed, Actuating, Confirming, Committed, Recovering, Halted };
struct State {
    Phase phase = Phase::Idle;
    uint64_t transaction = 0;
    uint32_t sequence = 0;
    uint32_t retry = 0;
    bool valid = true;
    bool operator==(const State& o) const {
        return phase == o.phase && transaction == o.transaction && sequence == o.sequence && retry == o.retry && valid == o.valid;
    }
};

struct Vote {
    bool recoverable = false;
    bool uncorrectable = false;
    State value{};
    int disagreeing_replicas = 0;
};

class TripleState {
public:
    explicit TripleState(State initial = {}) { replicas_.fill(initial); }
    void write_all(const State& next) { replicas_.fill(next); }
    void corrupt(int replica, const State& value) {
        if (replica < 0 || replica >= 3) throw std::out_of_range("replica");
        replicas_[static_cast<std::size_t>(replica)] = value;
    }
    Vote vote() const {
        for (int a = 0; a < 3; ++a) {
            int matches = 0;
            for (int b = 0; b < 3; ++b) if (replicas_[a] == replicas_[b]) ++matches;
            if (matches >= 2) {
                int disagree = 3 - matches;
                return Vote{true, false, replicas_[a], disagree};
            }
        }
        return Vote{false, true, {}, 3};
    }
    void scrub() {
        Vote v = vote();
        if (!v.recoverable) throw std::runtime_error("cannot scrub without a valid majority source");
        replicas_.fill(v.value);
    }
    const State& replica(int i) const { return replicas_.at(static_cast<std::size_t>(i)); }
private:
    std::array<State, 3> replicas_;
};

enum class TransactionOutcome { None, Committed, RetryRequired, Halted };

class Controller {
public:
    bool arm(uint64_t tx, uint32_t sequence) {
        if (halted_ || transaction_outcome_[tx] == TransactionOutcome::Committed) return false;
        State next{Phase::Armed, tx, sequence, 0, true};
        state_.write_all(next);
        return true;
    }
    bool start_fire(uint64_t tx) {
        auto v = assess();
        if (!v.recoverable || v.value.phase != Phase::Armed || v.value.transaction != tx) return false;
        State next = v.value; next.phase = Phase::Actuating;
        state_.write_all(next);
        return true;
    }
    TransactionOutcome confirm(uint64_t tx, bool sensor_confirmed) {
        Vote v = assess();
        if (!v.recoverable || v.value.transaction != tx || v.value.phase != Phase::Actuating) return halt(tx, "invalid-confirm-state");
        State checking = v.value; checking.phase = Phase::Confirming;
        state_.write_all(checking);
        if (sensor_confirmed) {
            State committed = checking; committed.phase = Phase::Committed;
            state_.write_all(committed);
            transaction_outcome_[tx] = TransactionOutcome::Committed;
            ++commits_;
            State idle{}; state_.write_all(idle);
            return TransactionOutcome::Committed;
        }
        return retry_or_halt(tx, checking);
    }
    void inject_state_fault(int replica, const State& state) { state_.corrupt(replica, state); }
    bool scrub() {
        Vote v = assess();
        if (!v.recoverable) return false;
        state_.scrub();
        ++scrubs_;
        return true;
    }
    Vote assess() {
        Vote v = state_.vote();
        if (v.uncorrectable) {
            halted_ = true;
            ++halts_;
            return v;
        }
        return v;
    }
    bool halted() const { return halted_; }
    int commits() const { return commits_; }
    int scrubs() const { return scrubs_; }
    int halts() const { return halts_; }
    TransactionOutcome outcome(uint64_t tx) const {
        auto it = transaction_outcome_.find(tx);
        return it == transaction_outcome_.end() ? TransactionOutcome::None : it->second;
    }
private:
    TransactionOutcome retry_or_halt(uint64_t tx, State current) {
        uint32_t& count = retry_count_[tx];
        if (count >= retry_budget_) return halt(tx, "retry-budget-exhausted");
        ++count;
        current.retry = count;
        current.phase = Phase::Recovering;
        state_.write_all(current);
        transaction_outcome_[tx] = TransactionOutcome::RetryRequired;
        ++retries_;
        State idle{}; state_.write_all(idle);
        return TransactionOutcome::RetryRequired;
    }
    TransactionOutcome halt(uint64_t tx, const std::string&) {
        halted_ = true;
        State stop{}; stop.phase = Phase::Halted; stop.transaction = tx;
        state_.write_all(stop);
        transaction_outcome_[tx] = TransactionOutcome::Halted;
        ++halts_;
        return TransactionOutcome::Halted;
    }
    TripleState state_;
    std::map<uint64_t, TransactionOutcome> transaction_outcome_;
    std::map<uint64_t, uint32_t> retry_count_;
    bool halted_ = false;
    int commits_ = 0, retries_ = 0, scrubs_ = 0, halts_ = 0;
    const uint32_t retry_budget_ = 2;
};

static void require(bool ok, const std::string& what) {
    if (!ok) throw std::runtime_error("ASSERTION FAILED: " + what);
}

static void one_fault_recovery() {
    std::cout << "\n[SCENARIO 1] one divergent protected-state replica is scrubbed\n";
    Controller c;
    require(c.arm(100, 7), "arm");
    require(c.start_fire(100), "fire");
    c.inject_state_fault(0, State{Phase::Halted, 100, 7, 0, false});
    Vote before = c.assess();
    require(before.recoverable && before.disagreeing_replicas == 1, "two valid replicas must outvote one fault");
    require(c.scrub(), "valid majority must support scrub");
    require(c.confirm(100, true) == TransactionOutcome::Committed, "operation may commit after protected recovery");
    require(c.commits() == 1 && !c.halted(), "one fault must not halt a protected controller");
    std::cout << "  majority recovery -> scrub -> confirmed commit\n";
}

static void two_fault_halt() {
    std::cout << "\n[SCENARIO 2] no majority source causes safe halt\n";
    Controller c;
    require(c.arm(200, 1), "arm");
    c.inject_state_fault(0, State{Phase::Halted, 200, 1, 0, false});
    c.inject_state_fault(1, State{Phase::Recovering, 200, 1, 2, false});
    Vote v = c.assess();
    require(v.uncorrectable, "three different replicas are uncorrectable");
    require(c.halted(), "uncorrectable controller state must halt");
    require(!c.scrub(), "no majority source means no fabricated scrub");
    std::cout << "  three-way disagreement -> halt, no invented restoration\n";
}

static void retry_and_exactly_once() {
    std::cout << "\n[SCENARIO 3] bounded retry and exactly-once transaction commit\n";
    Controller c;
    require(c.arm(300, 9), "first arm");
    require(c.start_fire(300), "first fire");
    require(c.confirm(300, false) == TransactionOutcome::RetryRequired, "ambiguous sensor must request retry");
    require(c.arm(300, 10), "retry arm with new sequence");
    require(c.start_fire(300), "retry fire");
    require(c.confirm(300, true) == TransactionOutcome::Committed, "retry may commit");
    require(!c.arm(300, 11), "committed transaction may not execute again");
    require(c.commits() == 1, "exactly one commit for transaction");
    std::cout << "  ambiguous result -> retry -> one committed transaction\n";
}

static void seeded_abstract_fault_campaign() {
    std::cout << "\n[SCENARIO 4] seeded campaign distinguishes one fault from two faults\n";
    std::mt19937 rng(20260921);
    std::uniform_int_distribution<int> which(0, 2);
    int corrected = 0;
    int halted = 0;
    for (int i = 0; i < 100; ++i) {
        Controller c;
        require(c.arm(1000 + i, 1), "campaign arm");
        require(c.start_fire(1000 + i), "campaign fire");
        if (i % 5 == 0) {
            c.inject_state_fault(0, State{Phase::Halted, 1, 1, 1, false});
            c.inject_state_fault(1, State{Phase::Recovering, 2, 2, 2, false});
            Vote v = c.assess();
            if (v.uncorrectable && c.halted()) ++halted;
        } else {
            c.inject_state_fault(which(rng), State{Phase::Recovering, 999999, 99, 99, false});
            if (c.scrub() && c.confirm(1000 + i, true) == TransactionOutcome::Committed) ++corrected;
        }
    }
    require(corrected == 80, "single-replica faults must still be corrected");
    require(halted == 20, "two-replica faults must halt instead of being invented as a majority");
    std::cout << "  corrected=" << corrected << " halted=" << halted << "\n";
}

} // namespace reliability

int main() {
    using namespace reliability;
    try {
        std::cout << "FEA V3 FZC protected-state simulator\n";
        std::cout << "Faults are abstract logical injections, not physical error-rate predictions.\n";
        one_fault_recovery();
        two_fault_halt();
        retry_and_exactly_once();
        seeded_abstract_fault_campaign();
        std::cout << "\nPASS: protected-state recovery invariants held.\n";
        std::cout << "NEXT EVIDENCE GATE: map these logical faults to physical retention, capture, sensing, and correlated-error models.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
