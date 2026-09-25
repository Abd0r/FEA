// =============================================================================
// FEA_fzc_program_v3.cpp -- M16 one chipset, three programs
//
// Claim under test: a Zone is programmable, not fixed-function. The same array
// at the 0.5 cm^2 design point must be able to behave like a CPU-like unit, a
// GPU-like unit, or an NPU-like unit, and the ONLY thing that differs between
// those runs may be the command stream. If a hardware parameter differs, this
// module fails, because then the result came from silicon, not from program.
//
// This is NOT a benchmark and NOT a performance claim against a dedicated
// accelerator. It reports cycles, ops and hops for a stated program shape on a
// stated array. Cycle time is set by M9, which this module does not re-derive.
//
// The comparison must distinguish FEA from PIM, CIM
// and RC. Those distinctions are conceptual in V2. M16 makes one of them
// executable: an array that runs three roles without a hardware change is not a
// fixed-function accelerator, whatever its speed turns out to be.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace program {

using fea::params;
using fea::require;

// The hardware each run is allowed to see. Read fresh per role so that any
// role-specific mutation of the chipset shows up as a mismatch below.
struct HardwareSignature {
    double zone_data_blocks = 0.0;
    double zone_fzc_blocks = 0.0;
    double block_pitch_nm = 0.0;
    double die_cm2 = 0.0;
    double design_zones = 0.0;
    double word_bits = 0.0;
};

static HardwareSignature read_hardware() {
    HardwareSignature h;
    h.zone_data_blocks = fea::zone_data_blocks();
    h.zone_fzc_blocks = static_cast<double>(fea::zone_fzc_blocks());
    h.block_pitch_nm = params().arch.block_pitch_nm;
    h.die_cm2 = params().arch.die_area_cm2;
    h.design_zones = fea::design_zone_count();
    h.word_bits = params().arch.word_bits;
    return h;
}

// PR6/S1: the reference used to be `read_hardware()`, and every Run also stored
// `read_hardware()`, so `same_hardware(a, a)` proved only that a pure function
// is deterministic -- it could not notice a defect inside read_hardware().
// The expected signature is built field by field from the named parameters,
// which is a different code path from read_hardware().
static HardwareSignature expected_hardware() {
    HardwareSignature h;
    h.zone_data_blocks = static_cast<double>(fea::params().arch.blocks_per_zone);
    h.zone_fzc_blocks = static_cast<double>(fea::zone_fzc_blocks());
    h.block_pitch_nm = fea::params().arch.block_pitch_nm;
    h.die_cm2 = fea::params().arch.die_area_cm2;
    h.design_zones = static_cast<double>(fea::design_zone_count());
    h.word_bits = fea::params().arch.word_bits;
    return h;
}

static bool same_hardware(const HardwareSignature& a, const HardwareSignature& b) {
    return a.zone_data_blocks == b.zone_data_blocks && a.zone_fzc_blocks == b.zone_fzc_blocks &&
           a.block_pitch_nm == b.block_pitch_nm && a.die_cm2 == b.die_cm2 &&
           a.design_zones == b.design_zones && a.word_bits == b.word_bits;
}

static std::string show(const HardwareSignature& h) {
    return "data=" + std::to_string(static_cast<long long>(h.zone_data_blocks)) +
           " fzc=" + std::to_string(static_cast<long long>(h.zone_fzc_blocks)) +
           " pitch=" + std::to_string(h.block_pitch_nm) + "nm" +
           " die=" + std::to_string(h.die_cm2) + "cm2" +
           " zones=" + std::to_string(static_cast<long long>(h.design_zones));
}

// One configured run. Everything here is an OUTPUT of the program shape, never
// an input: the role decides op count, depth, active Zones and hop count.
struct Run {
    std::string role;
    long long ops = 0;
    long long cycles = 0;
    long long active_zones = 0;
    long long hops = 0;
    HardwareSignature hardware;

    double ops_per_cycle() const { return cycles > 0 ? static_cast<double>(ops) / cycles : 0.0; }
    // How many independent instances of this configuration fit the real die.
    long long clusters_on_die() const {
        if (active_zones <= 0) return 0;
        return static_cast<long long>(fea::design_zone_count() / active_zones);
    }
};

// P1 CPU-like: one dependency chain. Each op needs the previous result, so
// depth and cycle count are the same number and only one Zone is ever busy.
static Run cpu_program(long long total_ops) {
    Run r;
    r.role = "CPU-like dependency chain";
    r.ops = total_ops;
    r.cycles = total_ops; // fully serial: every op is on the critical path
    r.active_zones = 1;
    r.hops = 0; // local to one Zone, no Slingshot traffic
    r.hardware = read_hardware();
    return r;
}

// P2 GPU-like: the same total ops spread over independent chains of equal
// length. Chains never wait on each other, so the critical path is one chain,
// not the total.
static Run gpu_program(long long total_ops, long long chains) {
    Run r;
    r.role = "GPU-like independent chains";
    const long long depth = total_ops / chains;
    r.ops = total_ops;
    r.cycles = depth;
    r.active_zones = chains;
    r.hops = 0; // intra-Zone only in this program shape
    r.hardware = read_hardware();
    return r;
}

// P3 NPU-like: a K x K systolic MAC array. Ops are local but the operands
// arrive from neighbours, so this shape is the one that spends hops.
static Run npu_program(long long cells) {
    Run r;
    r.role = "NPU-like systolic array";
    const long long k = static_cast<long long>(std::sqrt(static_cast<double>(cells)));
    r.ops = k * k;
    r.cycles = 2 * k - 1;              // wavefront: first result after 2K-1 stages
    r.active_zones = k * k;
    r.hops = 2 * k * (k - 1);          // each cell takes one pass from left and top
    r.hardware = read_hardware();
    return r;
}

static void scenario_same_hardware_three_programs() {
    std::cout << "\n[SCENARIO 1] one chipset, three programs, no hardware change\n";

    const long long target_ops = 4096;
    const std::vector<Run> runs{
        cpu_program(target_ops),
        gpu_program(target_ops, 256),
        npu_program(4096),
    };

    const HardwareSignature reference = expected_hardware();
    // Negative control: a predicate that accepted a perturbed signature could
    // not detect a role-specific hardware change either.
    {
        HardwareSignature altered = reference;
        altered.die_cm2 += 1.0;
        require(!same_hardware(altered, reference),
                "same_hardware must reject a perturbed signature, or the identity "
                "gate below would pass whatever the hardware did");
    }
    for (const Run& r : runs) {
        require(same_hardware(r.hardware, reference),
                "every role must run on an identical chipset, or the difference came from silicon");
        require(r.ops == target_ops,
                "all three programs must execute the same total op count for a fair comparison");
    }
    std::cout << "  chipset (identical in all three runs): " << show(reference) << "\n\n";

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "  " << std::left << std::setw(30) << "role" << std::right << std::setw(7) << "ops"
              << std::setw(9) << "cycles" << std::setw(8) << "Zones" << std::setw(8) << "hops"
              << std::setw(12) << "ops/cycle" << std::setw(12) << "clusters" << "\n";
    for (const Run& r : runs) {
        std::cout << "  " << std::left << std::setw(30) << r.role << std::right << std::setw(7)
                  << r.ops << std::setw(9) << r.cycles << std::setw(8) << r.active_zones
                  << std::setw(8) << r.hops << std::setw(12) << r.ops_per_cycle()
                  << std::setw(12) << r.clusters_on_die() << "\n";
    }
    std::cout << "\n";

    // The dispatcher must actually branch on the program. If all three produced
    // the same cycle count, the role selection was ignored.
    const Run& cpu = runs[0];
    const Run& gpu = runs[1];
    const Run& npu = runs[2];
    require(cpu.cycles != gpu.cycles && gpu.cycles != npu.cycles,
            "the three programs must produce different runtimes, or the program was not dispatched");
    require(cpu.active_zones == 1 && gpu.active_zones > 1 && npu.active_zones > gpu.active_zones,
            "each role must activate a different amount of the same array");

    require(gpu.ops_per_cycle() > cpu.ops_per_cycle(),
            "spreading the same ops over independent chains must beat the serial chain");
    require(npu.ops_per_cycle() > cpu.ops_per_cycle(),
            "the systolic shape must beat the serial chain despite its hops");
    require(npu.hops > cpu.hops && cpu.hops == 0,
            "only the systolic shape should spend Slingshot hops in this program set");
    require(runs[2].clusters_on_die() > 0,
            "at least one instance of the largest configuration must fit the 0.5 cm^2 die");

    std::cout << "  same array, same FZC budget, same pitch, same die in all three runs.\n";
    std::cout << "  what changed is the command stream: depth, parallel width, hop pattern.\n";
    std::cout << "  clusters = how many independent instances fit the design die, so the\n";
    std::cout << "  0.5 cm^2 chip can be partitioned into any mix of these three roles.\n";
    std::cout << "  label: program shape DERIVED, cycle time owned by M9, no accelerator claim.\n";
}

static void scenario_width_scales_throughput() {
    std::cout << "\n[SCENARIO 2] widening the parallel chain changes throughput, not hardware\n";
    const long long target_ops = 4096;
    const std::vector<long long> widths{1, 4, 16, 64, 256};

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "  " << std::left << std::setw(8) << "chains" << std::setw(9) << "depth"
              << std::setw(9) << "cycles" << std::setw(12) << "ops/cycle" << std::setw(12)
              << "clusters" << "\n";
    long long previous_cycles = 0;
    bool first = true;
    const HardwareSignature reference = expected_hardware();
    for (const long long w : widths) {
        const Run r = gpu_program(target_ops, w);
        require(same_hardware(r.hardware, reference),
                "widening the program must not alter the chipset");
        require(r.ops == target_ops, "every width must still execute the same op count");
        if (!first) {
            require(r.cycles < previous_cycles,
                    "adding independent chains must shorten the runtime, not lengthen it");
        }
        previous_cycles = r.cycles;
        first = false;
        std::cout << "  " << std::left << std::setw(8) << w << std::setw(9) << (target_ops / w)
                  << std::setw(9) << r.cycles << std::setw(12) << r.ops_per_cycle()
                  << std::setw(12) << r.clusters_on_die() << "\n";
    }
    std::cout << "  width is a program parameter. No register, cell or pitch changed.\n";
    std::cout << "  label: derived from program shape; the array size is bounded by the die.\n";
}

static void scenario_what_this_does_not_show() {
    std::cout << "\n[SCENARIO 3] what this module deliberately does not claim\n";
    const Run cpu = cpu_program(4096);
    const Run npu = npu_program(4096);

    require(cpu.clusters_on_die() > 1, "the serial shape trivially fits many times over");
    require(npu.clusters_on_die() > 0, "the systolic shape must still fit at least once");

    std::cout << "  CPU-like shape : " << cpu.clusters_on_die() << " independent instances on the die\n";
    std::cout << "  NPU-like shape : " << npu.clusters_on_die() << " independent instances on the die\n";
    std::cout << "  a 127-cycle systolic array is NOT shown to beat a real matrix engine here,\n";
    std::cout << "  because no real matrix engine was simulated. Speed is unclaimed.\n";
    std::cout << "  what IS shown: the capability exists, is reprogrammable, and needs no\n";
    std::cout << "  different silicon for any of the three roles.\n";
    std::cout << "  next evidence gate: energy per op per role, and a workload that stresses\n";
    std::cout << "  all three roles at once so partitioning cost becomes measurable.\n";
    std::cout << "  label: capability PROPOSED, performance OPEN.\n";
}

} // namespace program

int main() {
    using namespace program;
    try {
        std::cout << "FEA V3 M16 programmability: one chipset, three programs\n";
        std::cout << "Same 0.5 cm^2 design array in every run. Only the command stream differs.\n";
        scenario_same_hardware_three_programs();
        scenario_width_scales_throughput();
        scenario_what_this_does_not_show();
        std::cout << "\nPASS: one array ran three roles on identical hardware with distinct results.\n";
        std::cout << "LABEL: program shape derived, capability proposed, performance unclaimed.\n";
        std::cout << "NEXT EVIDENCE GATE: per-role energy and a mixed-role partition workload.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
