// =============================================================================
// FEA_fzc_v3.cpp -- FZC and unified Slingshot architectural simulator
//
// Scope: controller and network semantics only. This program deliberately does
// NOT model DB transport, irreversible capture, retention, a physical actuator,
// sensing physics, link attenuation, or chip power. Those are explicit inputs to
// later evidence gates, not invented constants.
//
// It verifies v0 architecture invariants through deterministic scenarios:
// boot, typed packet delivery, remote FIRE request, confirm-and-commit, retry,
// duplicate suppression, controller recovery, and safe halt.
// =============================================================================

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <deque>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fzc {

struct Coord {
    int x = 0;
    int y = 0;
    bool operator<(const Coord& o) const { return x != o.x ? x < o.x : y < o.y; }
    bool operator==(const Coord& o) const { return x == o.x && y == o.y; }
    bool operator!=(const Coord& o) const { return !(*this == o); }
};

static std::string text(Coord c) {
    return "(" + std::to_string(c.x) + "," + std::to_string(c.y) + ")";
}

enum class PacketClass { Boot, Control, Status, Ack, Nack, Refresh, Recovery, Diagnostic };
enum class Opcode { Boot, Arm, Fire, Confirm, Refresh, Recover, Halt };
enum class ControllerState { Reset, Booting, Idle, Validate, Armed, Actuating, Confirming, Recovering, Halted };
enum class Result { None, Accepted, Confirmed, Rejected, Ambiguous, Timeout, Fault };

static const char* name(PacketClass v) {
    switch (v) {
        case PacketClass::Boot: return "BOOT";
        case PacketClass::Control: return "CONTROL";
        case PacketClass::Status: return "STATUS";
        case PacketClass::Ack: return "ACK";
        case PacketClass::Nack: return "NACK";
        case PacketClass::Refresh: return "REFRESH";
        case PacketClass::Recovery: return "RECOVERY";
        case PacketClass::Diagnostic: return "DIAGNOSTIC";
    }
    return "?";
}
static const char* name(Opcode v) {
    switch (v) {
        case Opcode::Boot: return "BOOT";
        case Opcode::Arm: return "ARM";
        case Opcode::Fire: return "FIRE";
        case Opcode::Confirm: return "CONFIRM";
        case Opcode::Refresh: return "REFRESH";
        case Opcode::Recover: return "RECOVER";
        case Opcode::Halt: return "HALT";
    }
    return "?";
}
static const char* name(ControllerState v) {
    switch (v) {
        case ControllerState::Reset: return "RESET";
        case ControllerState::Booting: return "BOOTING";
        case ControllerState::Idle: return "IDLE";
        case ControllerState::Validate: return "VALIDATE";
        case ControllerState::Armed: return "ARMED";
        case ControllerState::Actuating: return "ACTUATING";
        case ControllerState::Confirming: return "CONFIRMING";
        case ControllerState::Recovering: return "RECOVERING";
        case ControllerState::Halted: return "HALTED";
    }
    return "?";
}
static const char* name(Result v) {
    switch (v) {
        case Result::None: return "NONE";
        case Result::Accepted: return "ACCEPTED";
        case Result::Confirmed: return "CONFIRMED";
        case Result::Rejected: return "REJECTED";
        case Result::Ambiguous: return "AMBIGUOUS";
        case Result::Timeout: return "TIMEOUT";
        case Result::Fault: return "FAULT";
    }
    return "?";
}

struct Packet {
    PacketClass type = PacketClass::Control;
    Opcode opcode = Opcode::Arm;
    Coord source{};
    Coord destination{};
    uint64_t transaction = 0;
    uint32_t attempt = 0;
    uint32_t sequence = 0;
    bool idempotent = true;
    bool protected_valid = true; // Abstract integrity result. Physical code/ECC comes later.
    std::string payload;
};

struct Trace {
    std::vector<std::string> lines;
    void add(const std::string& s) { lines.push_back(s); }
    bool contains(const std::string& token) const {
        return std::any_of(lines.begin(), lines.end(), [&](const std::string& x) { return x.find(token) != std::string::npos; });
    }
};

// Explicit abstract boundary. A later device model must replace this with a
// physical FZC-pattern -> actuator -> sensor model.
class Actuator {
public:
    virtual ~Actuator() = default;
    virtual Result request(Opcode op, const std::string& payload) = 0;
    virtual Result confirm() = 0;
    virtual void reset() = 0;
};

class ScriptedActuator final : public Actuator {
public:
    explicit ScriptedActuator(std::deque<Result> script) : script_(std::move(script)) {}
    Result request(Opcode op, const std::string&) override {
        if (op != Opcode::Fire) return Result::Rejected;
        ++requests_;
        return Result::Accepted;
    }
    Result confirm() override {
        if (script_.empty()) return Result::Fault;
        return pop();
    }
    void reset() override { ++resets_; }
    int requests() const { return requests_; }
    int resets() const { return resets_; }
private:
    Result pop() { Result r = script_.front(); script_.pop_front(); return r; }
    std::deque<Result> script_;
    int requests_ = 0;
    int resets_ = 0;
};

class Zone {
public:
    Zone(Coord id, Actuator& actuator, Trace& trace) : id_(id), actuator_(actuator), trace_(trace) {}
    Coord id() const { return id_; }
    ControllerState state() const { return state_; }
    bool booted() const { return booted_; }
    int commits() const { return commits_; }
    int halts() const { return halts_; }

    std::vector<Packet> receive(const Packet& p) {
        std::vector<Packet> out;
        trace_.add(text(id_) + " RECEIVE " + name(p.type) + " " + name(p.opcode) + " tx=" + std::to_string(p.transaction));
        if (state_ == ControllerState::Halted && p.opcode != Opcode::Boot) {
            out.push_back(reply(p, PacketClass::Nack, Opcode::Halt, "zone-halted"));
            return out;
        }
        if (!p.protected_valid) {
            enter_recovery("integrity-failure");
            out.push_back(reply(p, PacketClass::Nack, Opcode::Recover, "invalid-packet"));
            return out;
        }
        if (p.destination != id_) throw std::logic_error("wrong endpoint delivery");
        if (seen_.count({p.source, p.transaction, p.sequence})) {
            trace_.add(text(id_) + " DUPLICATE-SUPPRESSED tx=" + std::to_string(p.transaction));
            out.push_back(reply(p, PacketClass::Ack, p.opcode, "duplicate-already-processed"));
            return out;
        }
        seen_.insert({p.source, p.transaction, p.sequence});
        switch (p.opcode) {
            case Opcode::Boot: return boot(p);
            case Opcode::Arm: return arm(p);
            case Opcode::Fire: return fire(p);
            case Opcode::Refresh: return refresh(p);
            case Opcode::Recover: return recover(p);
            case Opcode::Halt: return halt(p, "remote-halt");
            case Opcode::Confirm: return {reply(p, PacketClass::Nack, Opcode::Confirm, "confirm-is-local-only")};
        }
        return out;
    }

private:
    struct SeenKey {
        Coord source;
        uint64_t transaction;
        uint32_t sequence;
        bool operator<(const SeenKey& o) const {
            if (source < o.source) return true;
            if (o.source < source) return false;
            return transaction != o.transaction ? transaction < o.transaction : sequence < o.sequence;
        }
    };

    Packet reply(const Packet& p, PacketClass type, Opcode opcode, std::string payload) const {
        return Packet{type, opcode, id_, p.source, p.transaction, p.attempt, p.sequence, true, true, std::move(payload)};
    }
    void transition(ControllerState target) {
        trace_.add(text(id_) + " " + name(state_) + " -> " + name(target));
        state_ = target;
    }
    void enter_recovery(const std::string& why) {
        transition(ControllerState::Recovering);
        trace_.add(text(id_) + " RECOVERY " + why);
        actuator_.reset();
    }
    std::vector<Packet> boot(const Packet& p) {
        if (state_ != ControllerState::Reset && state_ != ControllerState::Halted) return {reply(p, PacketClass::Nack, Opcode::Boot, "already-booted")};
        transition(ControllerState::Booting);
        booted_ = true;
        transition(ControllerState::Idle);
        return {reply(p, PacketClass::Ack, Opcode::Boot, "boot-committed")};
    }
    std::vector<Packet> arm(const Packet& p) {
        if (!booted_ || state_ != ControllerState::Idle) return {reply(p, PacketClass::Nack, Opcode::Arm, "not-idle-or-unbooted")};
        transition(ControllerState::Validate);
        armed_transaction_ = p.transaction;
        transition(ControllerState::Armed);
        return {reply(p, PacketClass::Ack, Opcode::Arm, "resources-reserved")};
    }
    std::vector<Packet> fire(const Packet& p) {
        if (!booted_ || state_ != ControllerState::Armed || armed_transaction_ != p.transaction)
            return {reply(p, PacketClass::Nack, Opcode::Fire, "not-armed-for-transaction")};
        transition(ControllerState::Actuating);
        Result started = actuator_.request(Opcode::Fire, p.payload);
        if (started != Result::Accepted) {
            enter_recovery("actuator-request-" + std::string(name(started)));
            transition(ControllerState::Idle);
            return {reply(p, PacketClass::Nack, Opcode::Fire, "actuator-not-accepted")};
        }
        transition(ControllerState::Confirming);
        Result observed = actuator_.confirm();
        trace_.add(text(id_) + " CONFIRM " + name(observed));
        if (observed == Result::Confirmed) {
            ++commits_;
            transition(ControllerState::Idle);
            return {reply(p, PacketClass::Status, Opcode::Confirm, "commit-confirmed")};
        }
        enter_recovery("confirm-" + std::string(name(observed)));
        if (p.attempt < max_retries_) {
            transition(ControllerState::Idle);
            Packet retry = p;
            retry.type = PacketClass::Nack;
            retry.opcode = Opcode::Recover;
            retry.payload = "retry-required";
            return {reply(p, PacketClass::Nack, Opcode::Recover, "retry-required")};
        }
        return halt(p, "retry-budget-exhausted");
    }
    std::vector<Packet> refresh(const Packet& p) {
        if (!booted_) return {reply(p, PacketClass::Nack, Opcode::Refresh, "unbooted")};
        if (p.payload != "valid-redundant-source") return {reply(p, PacketClass::Nack, Opcode::Refresh, "no-valid-source")};
        trace_.add(text(id_) + " REFRESH-COMMIT protected-state-from-redundant-source");
        return {reply(p, PacketClass::Ack, Opcode::Refresh, "scrub-committed")};
    }
    std::vector<Packet> recover(const Packet& p) {
        if (!booted_) return {reply(p, PacketClass::Nack, Opcode::Recover, "unbooted")};
        enter_recovery("external-request");
        transition(ControllerState::Idle);
        return {reply(p, PacketClass::Ack, Opcode::Recover, "recovery-complete")};
    }
    std::vector<Packet> halt(const Packet& p, const std::string& why) {
        actuator_.reset();
        transition(ControllerState::Halted);
        ++halts_;
        trace_.add(text(id_) + " HALT " + why);
        return {reply(p, PacketClass::Nack, Opcode::Halt, why)};
    }

    Coord id_;
    Actuator& actuator_;
    Trace& trace_;
    ControllerState state_ = ControllerState::Reset;
    bool booted_ = false;
    uint64_t armed_transaction_ = 0;
    int commits_ = 0;
    int halts_ = 0;
    const uint32_t max_retries_ = 2;
    std::set<SeenKey> seen_;
};

class Mesh {
public:
    void add(Zone& zone) { zones_.emplace(zone.id(), &zone); }
    std::vector<Packet> deliver(Packet packet) {
        std::vector<Packet> terminal;
        std::deque<Packet> q;
        q.push_back(packet);
        while (!q.empty()) {
            Packet p = q.front(); q.pop_front();
            ++hop_count_;
            auto it = zones_.find(p.destination);
            if (it == zones_.end()) throw std::logic_error("destination missing");
            std::vector<Packet> responses = it->second->receive(p);
            for (Packet& response : responses) {
                if (response.destination == p.source && (response.type == PacketClass::Ack || response.type == PacketClass::Nack || response.type == PacketClass::Status)) {
                    terminal.push_back(response);
                } else {
                    q.push_back(response);
                }
            }
        }
        return terminal;
    }
    int hop_count() const { return hop_count_; }
private:
    std::map<Coord, Zone*> zones_;
    int hop_count_ = 0;
};

static Packet packet(PacketClass kind, Opcode op, Coord from, Coord to, uint64_t tx, uint32_t seq, std::string body = "") {
    return Packet{kind, op, from, to, tx, 0, seq, true, true, std::move(body)};
}

static void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error("ASSERTION FAILED: " + message);
}

static void print_terminal(const std::vector<Packet>& r) {
    for (const Packet& p : r) {
        std::cout << "  reply " << name(p.type) << " " << name(p.opcode)
                  << " " << text(p.source) << " -> " << text(p.destination)
                  << " tx=" << p.transaction << " " << p.payload << "\n";
    }
}

static void scenario_happy_path() {
    std::cout << "\n[SCENARIO 1] boot -> ARM -> FIRE -> CONFIRM -> commit\n";
    Trace trace;
    ScriptedActuator actuator({Result::Confirmed});
    Zone source({0, 0}, actuator, trace), destination({1, 0}, actuator, trace);
    Mesh mesh; mesh.add(source); mesh.add(destination);
    print_terminal(mesh.deliver(packet(PacketClass::Boot, Opcode::Boot, {0, 0}, {1, 0}, 1, 1)));
    print_terminal(mesh.deliver(packet(PacketClass::Control, Opcode::Arm, {0, 0}, {1, 0}, 2, 1)));
    print_terminal(mesh.deliver(packet(PacketClass::Control, Opcode::Fire, {0, 0}, {1, 0}, 2, 2, "target-word=7")));
    require(destination.commits() == 1, "confirmed FIRE must commit exactly once");
    require(destination.state() == ControllerState::Idle, "confirmed operation returns to IDLE");
    require(trace.contains("CONFIRM CONFIRMED"), "confirmation trace required");
}

static void scenario_retry_and_duplicate() {
    std::cout << "\n[SCENARIO 2] ambiguous confirm -> recovery -> retry, duplicate suppression\n";
    Trace trace;
    ScriptedActuator actuator({Result::Ambiguous, Result::Confirmed});
    Zone source({0, 0}, actuator, trace), destination({1, 0}, actuator, trace);
    Mesh mesh; mesh.add(source); mesh.add(destination);
    mesh.deliver(packet(PacketClass::Boot, Opcode::Boot, {0, 0}, {1, 0}, 10, 1));
    mesh.deliver(packet(PacketClass::Control, Opcode::Arm, {0, 0}, {1, 0}, 11, 1));
    auto first = mesh.deliver(packet(PacketClass::Control, Opcode::Fire, {0, 0}, {1, 0}, 11, 2));
    require(first.size() == 1 && first[0].type == PacketClass::Nack, "ambiguous confirm must NACK");
    mesh.deliver(packet(PacketClass::Control, Opcode::Arm, {0, 0}, {1, 0}, 12, 1));
    auto second = mesh.deliver(packet(PacketClass::Control, Opcode::Fire, {0, 0}, {1, 0}, 12, 2));
    require(second.size() == 1 && second[0].type == PacketClass::Status, "retry confirmation must commit");
    auto duplicate = mesh.deliver(packet(PacketClass::Control, Opcode::Fire, {0, 0}, {1, 0}, 12, 2));
    require(destination.commits() == 1, "duplicate may not create a second commit");
    require(duplicate[0].payload == "duplicate-already-processed", "duplicate response must be explicit");
    require(actuator.resets() >= 1, "ambiguous action must reset actuator before retry");
}

static void scenario_refresh_and_safe_halt() {
    std::cout << "\n[SCENARIO 3] refresh source rule and safe halt\n";
    Trace trace;
    ScriptedActuator actuator({Result::Fault});
    Zone source({0, 0}, actuator, trace), destination({1, 0}, actuator, trace);
    Mesh mesh; mesh.add(source); mesh.add(destination);
    mesh.deliver(packet(PacketClass::Boot, Opcode::Boot, {0, 0}, {1, 0}, 20, 1));
    auto bad_refresh = mesh.deliver(packet(PacketClass::Refresh, Opcode::Refresh, {0, 0}, {1, 0}, 21, 1, "lost-state"));
    require(bad_refresh[0].type == PacketClass::Nack, "refresh without redundant source must fail");
    auto good_refresh = mesh.deliver(packet(PacketClass::Refresh, Opcode::Refresh, {0, 0}, {1, 0}, 22, 1, "valid-redundant-source"));
    require(good_refresh[0].type == PacketClass::Ack, "refresh with valid source must commit");
    mesh.deliver(packet(PacketClass::Control, Opcode::Arm, {0, 0}, {1, 0}, 23, 1));
    auto failure = mesh.deliver(packet(PacketClass::Control, Opcode::Fire, {0, 0}, {1, 0}, 23, 2));
    require(failure[0].type == PacketClass::Nack, "faulted actuator must NACK");
    require(destination.state() == ControllerState::Idle, "first failure remains recoverable");
}

static void scenario_retry_exhausted() {
    std::cout << "\n[SCENARIO 4] retry budget exhaustion halts\n";
    Trace trace;
    ScriptedActuator actuator({Result::Fault});
    Zone source({0, 0}, actuator, trace), destination({1, 0}, actuator, trace);
    Mesh mesh;
    mesh.add(source);
    mesh.add(destination);
    mesh.deliver(packet(PacketClass::Boot, Opcode::Boot, {0, 0}, {1, 0}, 30, 1));
    mesh.deliver(packet(PacketClass::Control, Opcode::Arm, {0, 0}, {1, 0}, 31, 1));
    Packet exhausted{PacketClass::Control, Opcode::Fire, {0, 0}, {1, 0}, 31, 2, 2, true, true, "target-word=7"};
    auto halted = mesh.deliver(exhausted);
    require(halted.size() == 1 && halted[0].payload == "retry-budget-exhausted", "attempt at the retry ceiling must halt");
    require(destination.state() == ControllerState::Halted && destination.halts() == 1, "exhausted retries must halt the zone");
    std::cout << "  attempt=2 fault -> halt\n";
}

} // namespace fzc

int main() {
    using namespace fzc;
    try {
        std::cout << "FEA V3 FZC/Slingshot architectural simulator\n";
        std::cout << "Scope: protocol semantics only. Device actuation is abstract.\n";
        scenario_happy_path();
        scenario_retry_and_duplicate();
        scenario_refresh_and_safe_halt();
        scenario_retry_exhausted();
        std::cout << "\nPASS: all v0 architectural invariants held in deterministic scenarios.\n";
        std::cout << "NEXT EVIDENCE GATE: replace ScriptedActuator with a validated physical model.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
