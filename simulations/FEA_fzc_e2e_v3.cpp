// =============================================================================
// FEA_fzc_e2e_v3.cpp -- one FZC transaction computed across every gate
//
// Each gate below is calculated inside this file:
//   dual-rail recognition score
//   screened geometry bias
//   clocked storage integral
//   sensor reading = stored + declared offset
//   three-replica majority vote
//   hop walk over an explicit link list
//   duplicate check against a commit log
//
// Declared parameters are not silicon measurements.
// =============================================================================

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace e2e {

struct Vec2 { double x = 0.0; double y = 0.0; };
static const Vec2 kPlusA{0.5, 0.5};
static const Vec2 kPlusB{-0.5, -0.5};
static const Vec2 kMinusA{-0.5, 0.5};
static const Vec2 kMinusB{0.5, -0.5};

using Rails = std::array<int, 6>;
struct Replica {
    int phase = 1;
    int tx = 7;
    bool operator==(const Replica& o) const { return phase == o.phase && tx == o.tx; }
};
struct Link { bool up = true; };
struct Outcome { double lead = 1.0; double actuator = 0.0; double stored = 0.0; double reservoir_energy = 0.0; };

struct Request {
    Rails rails{1, 0, 0, 1, 0, 1}; // target bits 0,1,1
    double separation = 1.0;
    double lambda = 100.0;
    double clock = 1.2;
    double sensor_offset = 0.0;
    std::array<Replica, 3> replicas{};
    std::vector<Link> links{Link{true}, Link{true}};
    uint64_t tx = 7;
};

struct Report {
    double score = 0.0;
    double bias = 0.0;
    double stored = 0.0;
    double reading = 0.0;
    int hops = 0;
    bool committed = false;
    std::string stop = "none";
};

class CommitLog {
public:
    bool seen(uint64_t tx) const { return ids_.count(tx) != 0; }
    void add(uint64_t tx) { ids_.insert(tx); }
private:
    std::set<uint64_t> ids_;
};

static void require(bool ok, const std::string& what) {
    if (!ok) throw std::runtime_error("ASSERTION FAILED: " + what);
}

static double recognize(const Rails& rails) {
    const int target[3] = {0, 1, 1};
    double score = 0.0;
    for (int bit = 0; bit < 3; ++bit) {
        const int selected = target[bit];
        score += rails[2 * bit + selected] - rails[2 * bit + (1 - selected)];
    }
    return score;
}

static double distance(Vec2 a, Vec2 b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

static double bias_for(double separation, double lambda) {
    const Vec2 command{-separation, 0.5};
    auto couple = [&](Vec2 site) { return std::exp(-distance(command, site) / lambda) / distance(command, site); };
    return (couple(kMinusA) + couple(kMinusB)) - (couple(kPlusA) + couple(kPlusB));
}

static double clock_at(double amplitude, double time) {
    if (amplitude <= 0.0 || time < 0.0 || time > 1.0) return 0.0;
    return amplitude * std::sin(3.14159265358979323846 * time);
}

static Outcome integrate(double bias, double amplitude) {
    Outcome state;
    const double dt = 0.001;
    for (int step = 0; step < 1500; ++step) {
        const double time = step * dt;
        const double clock = clock_at(amplitude, time);
        const double barrier = std::max(0.0, 1.5 - bias - clock);
        const double enter = dt * 2.0 * std::exp(-barrier) * state.lead;
        const double leak = dt * 0.04 * state.lead;
        const bool falling = clock_at(amplitude, time) > clock_at(amplitude, time + dt);
        const double relax = dt * ((bias > 0.5 && falling) ? 8.0 : 0.0) * std::max(state.actuator, 0.0);
        const double escape = dt * 0.05 * std::exp(-std::abs(bias)) * std::max(state.actuator, 0.0);
        state.lead -= enter + leak;
        state.actuator += enter - relax - escape;
        state.stored += relax;
        state.reservoir_energy += relax + escape;
    }
    return state;
}

static bool majority(const std::array<Replica, 3>& replicas, Replica& agreed) {
    for (int a = 0; a < 3; ++a) {
        int matches = 0;
        for (int b = 0; b < 3; ++b) if (replicas[a] == replicas[b]) ++matches;
        if (matches >= 2) {
            agreed = replicas[a];
            return true;
        }
    }
    return false;
}

static int walk(const std::vector<Link>& links) {
    int hops = 0;
    for (const Link& link : links) {
        if (!link.up) return -1;
        ++hops;
    }
    return hops;
}

static Report run(const Request& request, CommitLog& log) {
    Report report;
    report.score = recognize(request.rails);
    if (report.score <= 1.0) {
        report.stop = "bad-code";
        return report;
    }
    report.bias = bias_for(request.separation, request.lambda);
    if (report.bias <= 0.5) {
        report.stop = "weak-bias";
        return report;
    }
    const Outcome physical = integrate(report.bias, request.clock);
    report.stored = physical.stored;
    if (physical.stored <= 0.5 || physical.reservoir_energy <= 0.0) {
        report.stop = request.clock <= 0.0 ? "no-clock" : "no-capture";
        return report;
    }
    report.reading = physical.stored + request.sensor_offset;
    if (report.reading < 0.5) {
        report.stop = "sensor-reject";
        return report;
    }
    Replica agreed;
    if (!majority(request.replicas, agreed) || agreed.tx != static_cast<int>(request.tx)) {
        report.stop = "no-majority";
        return report;
    }
    report.hops = walk(request.links);
    if (report.hops < 0) {
        report.stop = "link-down";
        return report;
    }
    if (log.seen(request.tx)) {
        report.stop = "duplicate";
        return report;
    }
    log.add(request.tx);
    report.committed = true;
    report.stop = "none";
    return report;
}

static void success_path() {
    std::cout << "\n[SCENARIO 1] one valid transaction computes every gate and commits once\n";
    CommitLog log;
    const Report done = run(Request{}, log);
    require(done.stop == "none", "clean path must pass every computed gate");
    require(done.score == 3.0 && done.bias > 0.5 && done.stored > 0.5, "score, bias, and storage must be computed");
    require(done.reading >= 0.5 && done.hops == 2 && done.committed, "sensor, vote, and link walk must allow one commit");
    std::cout << "  score=" << done.score << " bias=" << done.bias << " stored=" << done.stored
              << " reading=" << done.reading << " hops=" << done.hops << "\n";
}

static void failures_do_not_commit() {
    std::cout << "\n[SCENARIO 2] a failed computed gate stops the path before commit\n";
    CommitLog log;
    Request bad;
    bad.rails = {1, 0, 1, 0, 0, 1};
    Request weak;
    weak.lambda = 0.35;
    Request dark;
    dark.clock = 0.0;
    Request noisy;
    noisy.sensor_offset = -1.0;
    Request split;
    split.replicas = {Replica{1, 7}, Replica{2, 7}, Replica{3, 8}};
    Request down;
    down.links = {Link{true}, Link{false}};
    const Request inputs[] = {bad, weak, dark, noisy, split, down};
    const char* expected[] = {"bad-code", "weak-bias", "no-clock", "sensor-reject", "no-majority", "link-down"};
    for (int i = 0; i < 6; ++i) {
        const Report report = run(inputs[i], log);
        require(report.stop == expected[i], "transaction stopped at the wrong computed gate");
        require(!report.committed, "stopped transaction committed");
        std::cout << "  " << report.stop << " -> no commit\n";
    }
    Request first;
    first.tx = 42;
    first.replicas = {Replica{1, 42}, Replica{1, 42}, Replica{1, 42}};
    Request second = first;
    require(run(first, log).committed, "first transaction id must commit");
    const Report duplicate = run(second, log);
    require(duplicate.stop == "duplicate" && !duplicate.committed, "commit log must reject the same transaction id");
    std::cout << "  duplicate -> no commit\n";
}

static void majority_recovers() {
    std::cout << "\n[SCENARIO 3] one bad replica still allows one commit\n";
    CommitLog log;
    Request one;
    one.tx = 99;
    one.replicas = {Replica{9, 1}, Replica{1, 99}, Replica{1, 99}};
    const Report done = run(one, log);
    require(done.committed && done.hops == 2, "two matching replicas must outvote one and still deliver");
    std::cout << "  commits=1 hops=" << done.hops << "\n";
}

} // namespace e2e

int main() {
    using namespace e2e;
    try {
        std::cout << "FEA V3 FZC end-to-end composed simulation\n";
        std::cout << "Recognition, storage, sensing, vote, and delivery are computed. Parameters are declared.\n";
        success_path();
        failures_do_not_commit();
        majority_recovers();
        std::cout << "\nPASS: the path committed only when every computed gate passed.\n";
        std::cout << "NEXT EVIDENCE GATE: replace these declared parameters with measured device rates.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
