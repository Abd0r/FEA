// =============================================================================
// FEA_slingshot_v3.cpp -- unified Slingshot architectural network simulator
//
// Scope: finite-buffer, typed-packet, coordinate-routing semantics in discrete
// arbitration rounds. A round is NOT a physical time unit. This program makes
// no claim about DBW delay, link energy, capture, restoration, or chip bandwidth.
// Those require a physical-link model.
// =============================================================================

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace sling {

struct Coord {
    int x = 0, y = 0;
    bool operator<(const Coord& o) const { return x != o.x ? x < o.x : y < o.y; }
    bool operator==(const Coord& o) const { return x == o.x && y == o.y; }
};
static std::string text(Coord c) { return "(" + std::to_string(c.x) + "," + std::to_string(c.y) + ")"; }

enum class Class { Recovery, Ack, Refresh, Boot, Control, Status, Data, Diagnostic };
enum class Terminal { Active, Delivered, RejectedIngress, RejectedLink, Expired };

static int priority(Class c) { return static_cast<int>(c); }

struct Packet {
    uint64_t id = 0;
    Class type = Class::Data;
    Coord source{}, destination{};
    uint32_t sequence = 0;
    int hop_limit = 32;
    int hops = 0;
    int injected_round = 0;
    int delivered_round = -1;
    Terminal terminal = Terminal::Active;
    std::string reason;
};

struct Link {
    Coord a{}, b{};
    bool live = true;
    bool operator<(const Link& o) const {
        return a < o.a ? true : o.a < a ? false : b < o.b;
    }
};

class Network {
public:
    Network(int width, int height, std::size_t ingress_capacity)
        : width_(width), height_(height), ingress_capacity_(ingress_capacity) {
        if (width < 1 || height < 1 || ingress_capacity < 1) throw std::invalid_argument("invalid network geometry");
        for (int y = 0; y < height_; ++y)
            for (int x = 0; x < width_; ++x) queues_[{x, y}] = {};
    }

    uint64_t inject(Class type, Coord source, Coord destination, uint32_t sequence = 0, int hop_limit = 32) {
        ensure(source); ensure(destination);
        Packet p{next_id_++, type, source, destination, sequence, hop_limit, 0, round_, -1, Terminal::Active, ""};
        packets_[p.id] = p;
        auto& q = queues_[source];
        if (q.size() >= ingress_capacity_) {
            packets_[p.id].terminal = Terminal::RejectedIngress;
            packets_[p.id].reason = "origin-ingress-full";
            ++rejected_ingress_;
        } else {
            q.push_back(p.id);
        }
        return p.id;
    }

    void set_link(Coord from, Coord to, bool live) {
        ensure(from); ensure(to);
        if (manhattan(from, to) != 1) throw std::invalid_argument("link endpoints must be neighbours");
        links_[{from, to}].live = live;
    }

    void tick() {
        ++round_;
        std::map<Coord, std::vector<uint64_t>> arrivals;
        std::set<Link> used;
        for (auto& entry : queues_) {
            const Coord at = entry.first;
            auto& q = entry.second;
            if (q.empty()) continue;
            const std::size_t idx = select(q);
            const uint64_t id = q[idx];
            Packet& p = packets_.at(id);
            if (p.destination == at) {
                q.erase(q.begin() + static_cast<std::ptrdiff_t>(idx));
                deliver(p);
                continue;
            }
            if (p.hops >= p.hop_limit) {
                q.erase(q.begin() + static_cast<std::ptrdiff_t>(idx));
                p.terminal = Terminal::Expired;
                p.reason = "hop-limit";
                ++expired_;
                continue;
            }
            const Coord next = route(at, p.destination);
            Link link{at, next, true};
            auto found = links_.find(link);
            const bool live = found == links_.end() ? true : found->second.live;
            if (!live) {
                q.erase(q.begin() + static_cast<std::ptrdiff_t>(idx));
                p.terminal = Terminal::RejectedLink;
                p.reason = "link-down " + text(at) + "->" + text(next);
                ++rejected_link_;
                continue;
            }
            if (used.count(link) || queues_.at(next).size() + arrivals[next].size() >= ingress_capacity_) {
                ++blocked_;
                continue;
            }
            q.erase(q.begin() + static_cast<std::ptrdiff_t>(idx));
            used.insert(link);
            ++p.hops;
            arrivals[next].push_back(id);
        }
        for (auto& entry : arrivals)
            for (uint64_t id : entry.second) queues_.at(entry.first).push_back(id);
        // Delivery consumes a local arbitration round, making endpoint handling explicit.
        for (auto& entry : queues_) {
            auto& q = entry.second;
            for (std::size_t i = 0; i < q.size();) {
                Packet& p = packets_.at(q[i]);
                if (p.destination == entry.first) {
                    q.erase(q.begin() + static_cast<std::ptrdiff_t>(i));
                    deliver(p);
                } else {
                    ++i;
                }
            }
        }
        assert_conservation();
    }

    void run_until_terminal(uint64_t id, int maximum_rounds) {
        for (int i = 0; i < maximum_rounds && packets_.at(id).terminal == Terminal::Active; ++i) tick();
        if (packets_.at(id).terminal == Terminal::Active) throw std::runtime_error("packet remained active past test limit");
    }

    const Packet& packet(uint64_t id) const { return packets_.at(id); }
    int round() const { return round_; }
    int blocked() const { return blocked_; }
    int rejected_ingress() const { return rejected_ingress_; }
    int rejected_link() const { return rejected_link_; }

private:
    void ensure(Coord c) const {
        if (c.x < 0 || c.x >= width_ || c.y < 0 || c.y >= height_) throw std::out_of_range("coordinate outside mesh");
    }
    static int manhattan(Coord a, Coord b) { return std::abs(a.x - b.x) + std::abs(a.y - b.y); }
    Coord route(Coord at, Coord destination) const {
        if (at.x != destination.x) return {at.x + (destination.x > at.x ? 1 : -1), at.y};
        return {at.x, at.y + (destination.y > at.y ? 1 : -1)};
    }
    std::size_t select(const std::deque<uint64_t>& q) const {
        std::size_t selected = 0;
        for (std::size_t i = 1; i < q.size(); ++i) {
            const Packet& a = packets_.at(q[i]);
            const Packet& b = packets_.at(q[selected]);
            if (priority(a.type) < priority(b.type) || (priority(a.type) == priority(b.type) && a.id < b.id)) selected = i;
        }
        return selected;
    }
    void deliver(Packet& p) {
        if (p.terminal != Terminal::Active) throw std::logic_error("only active packet can deliver");
        p.terminal = Terminal::Delivered;
        p.delivered_round = round_;
        ++delivered_;
    }
    void assert_conservation() const {
        std::set<uint64_t> buffered;
        for (const auto& entry : queues_) for (uint64_t id : entry.second) {
            if (!buffered.insert(id).second) throw std::logic_error("packet appears in two buffers");
            if (packets_.at(id).terminal != Terminal::Active) throw std::logic_error("terminal packet remains buffered");
        }
        int terminal_count = 0;
        for (const auto& entry : packets_) {
            const Packet& p = entry.second;
            if (p.terminal == Terminal::Active) {
                if (!buffered.count(p.id)) throw std::logic_error("active packet lost outside buffers");
            } else ++terminal_count;
        }
        if (terminal_count != delivered_ + rejected_ingress_ + rejected_link_ + expired_)
            throw std::logic_error("terminal packet accounting mismatch");
    }

    int width_, height_;
    std::size_t ingress_capacity_;
    int round_ = 0;
    uint64_t next_id_ = 1;
    int delivered_ = 0, rejected_ingress_ = 0, rejected_link_ = 0, expired_ = 0, blocked_ = 0;
    std::map<Coord, std::deque<uint64_t>> queues_;
    std::map<Link, Link> links_;
    std::map<uint64_t, Packet> packets_;
};

static void require(bool result, const std::string& message) {
    if (!result) throw std::runtime_error("ASSERTION FAILED: " + message);
}

static void scenario_priority() {
    std::cout << "\n[SCENARIO 1] priority: RECOVERY before DATA at one finite port\n";
    Network n(3, 1, 4);
    const uint64_t data = n.inject(Class::Data, {0, 0}, {2, 0});
    const uint64_t recovery = n.inject(Class::Recovery, {0, 0}, {2, 0});
    n.run_until_terminal(data, 20);
    n.run_until_terminal(recovery, 20);
    require(n.packet(recovery).terminal == Terminal::Delivered, "recovery must deliver");
    require(n.packet(data).terminal == Terminal::Delivered, "data must deliver");
    require(n.packet(recovery).delivered_round < n.packet(data).delivered_round, "recovery must win arbitration");
    std::cout << "  recovery delivered in round " << n.packet(recovery).delivered_round
              << "; data delivered in round " << n.packet(data).delivered_round << "\n";
}

static void scenario_backpressure() {
    std::cout << "\n[SCENARIO 2] finite ingress buffers and backpressure\n";
    Network n(3, 1, 1);
    const uint64_t first = n.inject(Class::Control, {0, 0}, {2, 0});
    const uint64_t second = n.inject(Class::Data, {0, 0}, {2, 0});
    require(n.packet(second).terminal == Terminal::RejectedIngress, "full ingress must reject explicitly");
    n.run_until_terminal(first, 20);
    require(n.packet(first).terminal == Terminal::Delivered, "accepted packet must deliver");
    std::cout << "  explicit ingress rejection count=" << n.rejected_ingress() << "\n";
}

static void scenario_fault_and_hops() {
    std::cout << "\n[SCENARIO 3] failed link is explicit; hop count is conserved\n";
    Network n(3, 3, 3);
    n.set_link({1, 0}, {2, 0}, false);
    const uint64_t fault = n.inject(Class::Control, {0, 0}, {2, 0});
    n.run_until_terminal(fault, 20);
    require(n.packet(fault).terminal == Terminal::RejectedLink, "failed link must not become silent loss");
    require(n.packet(fault).reason.find("link-down") != std::string::npos, "failure reason must name link");

    Network clean(3, 3, 3);
    const uint64_t path = clean.inject(Class::Data, {0, 0}, {2, 2});
    clean.run_until_terminal(path, 30);
    require(clean.packet(path).terminal == Terminal::Delivered, "clean path must deliver");
    require(clean.packet(path).hops == 4, "XY route must use Manhattan hop count");
    std::cout << "  failed-link rejection and 4-hop delivery verified\n";
}

static void scenario_refresh_priority() {
    std::cout << "\n[SCENARIO 4] refresh deadline traffic outranks ordinary control/data\n";
    Network n(3, 1, 4);
    const uint64_t control = n.inject(Class::Control, {0, 0}, {2, 0});
    const uint64_t refresh = n.inject(Class::Refresh, {0, 0}, {2, 0});
    n.run_until_terminal(control, 20);
    n.run_until_terminal(refresh, 20);
    require(n.packet(refresh).delivered_round < n.packet(control).delivered_round, "refresh must outrank ordinary control");
    std::cout << "  refresh delivered in round " << n.packet(refresh).delivered_round
              << "; control delivered in round " << n.packet(control).delivered_round << "\n";
}

} // namespace sling

int main() {
    using namespace sling;
    try {
        std::cout << "FEA V3 unified Slingshot architectural simulator\n";
        std::cout << "Unit: arbitration rounds. No physical time, energy, or bandwidth claim.\n";
        scenario_priority();
        scenario_backpressure();
        scenario_fault_and_hops();
        scenario_refresh_priority();
        std::cout << "\nPASS: finite-buffer routing invariants held in all v0 scenarios.\n";
        std::cout << "NEXT EVIDENCE GATE: physical link, restoration, and packet-energy model.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
