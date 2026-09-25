// =============================================================================
// FEA_secded_v3.cpp -- M8 SECDED code cost, correction capability, and limits
//
// V2 asserts a standard SECDED Hamming code at the Word
// level with "~2-3x overhead" but never simulates it and never reports its
// area, power, latency, or residual failure. This module derives the code from
// the Hamming bound for the stated Word width, charges the overhead, and checks
// it against V2's claimed factor. It also tests what SECDED cannot do: correct
// a multi-bit burst of the kind a common-mode failure produces.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace secded {

using fea::params;
using fea::require;

struct Code {
    int data_bits = 64;
    int hamming_parity = 0;   // single-error-correcting Hamming parity bits
    int overall_parity = 1;   // SECDED extra overall parity bit
    int total_parity = 0;
    int codeword = 0;
    double overhead_ratio = 0.0;
};

// Smallest r with 2^r >= k + r + 1.
static int hamming_parity_bits(int k) {
    // 1LL not 1: `1 << r` overflows signed 32-bit int at r=31 and is undefined
    // behaviour for r>=31, and the old loop ran to r<=64. Capped at 60, far
    // beyond any Word width this suite models, so the shift always stays in
    // range for a 64-bit long long.
    for (int r = 1; r <= 60; ++r) {
        if ((1LL << r) >= k + r + 1) return r;
    }
    throw std::runtime_error("no Hamming parity solution for this word width");
}

static Code make_code(int data_bits) {
    Code c;
    c.data_bits = data_bits;
    c.hamming_parity = hamming_parity_bits(data_bits);
    c.overall_parity = 1;
    c.total_parity = c.hamming_parity + c.overall_parity;
    c.codeword = c.data_bits + c.total_parity;
    c.overhead_ratio = static_cast<double>(c.codeword) / c.data_bits;
    return c;
}

// ---- real SECDED codec for the derived 71-bit Hamming portion + overall bit ----
// Positions are 1-indexed. Power-of-two positions 1,2,4,...,64 hold parity.
static bool is_parity_position(int pos) {
    return pos > 0 && (pos & (pos - 1)) == 0 && pos <= 64;
}

// Encode 64 data bits into a 72-bit SECDED codeword.
static std::vector<int> encode(uint64_t data) {
    const int ham = 71;
    std::vector<int> cw(ham + 1, 0);
    int bit = 0;
    for (int pos = 1; pos <= ham; ++pos) {
        if (is_parity_position(pos)) continue;
        cw[pos - 1] = static_cast<int>((data >> bit) & 1ULL);
        ++bit;
    }
    for (int p = 1; p <= 64; p <<= 1) {
        int x = 0;
        for (int j = 1; j <= ham; ++j)
            if (((j & p) != 0) && j != p) x ^= cw[j - 1];
        cw[p - 1] = x;
    }
    int all = 0;
    for (int i = 0; i < ham; ++i) all ^= cw[i];
    cw[ham] = all;
    return cw;
}

static uint64_t extract_data(const std::vector<int>& cw) {
    const int ham = 71;
    uint64_t data = 0;
    int bit = 0;
    for (int pos = 1; pos <= ham; ++pos) {
        if (is_parity_position(pos)) continue;
        data |= (static_cast<uint64_t>(cw[pos - 1]) << bit);
        ++bit;
    }
    return data;
}

enum DecodeResult { kClean = 0, kCorrected = 1, kDetectedUncorrectable = 2 };

// Returns whether a correction was attempted. Leaves decoded data in data_out.
static DecodeResult decode(std::vector<int> cw, uint64_t& data_out) {
    const int ham = 71;
    int syndrome = 0;
    for (int p = 1; p <= 64; p <<= 1) {
        int x = 0;
        for (int j = 1; j <= ham; ++j)
            if (((j & p) != 0) && j != p) x ^= cw[j - 1];
        if (cw[p - 1] != x) syndrome |= p;
    }
    int overall = 0;
    for (int i = 0; i < ham; ++i) overall ^= cw[i];
    const bool overall_mismatch = (overall != cw[ham]);

    DecodeResult result = kClean;
    if (syndrome == 0 && !overall_mismatch) {
        result = kClean;
    } else if (syndrome == 0 && overall_mismatch) {
        cw[ham] ^= 1;          // only the overall parity bit was wrong
        result = kCorrected;
    } else if (syndrome != 0 && overall_mismatch) {
        // GUARD: syndrome is built from parity positions {1,2,4,8,16,32,64} OR'd
        // together, so for multi-bit errors it can reach 127 while the codeword
        // only has 71 Hamming positions. Indexing at cw[syndrome-1] with
        // syndrome > 71 is out of bounds (ASan: heap-buffer-overflow). A real
        // SECDED decoder only corrects when syndrome names a valid position.
        if (syndrome >= 1 && syndrome <= ham) {
            cw[syndrome - 1] ^= 1; // single error at the syndromed position
            result = kCorrected;
        } else {
            result = kDetectedUncorrectable; // syndrome names no real position
        }
    } else {
        result = kDetectedUncorrectable; // syndrome set, parity even => even count
    }
    data_out = extract_data(cw);
    return result;
}

static void scenario_derive_the_code() {
    std::cout << "\n[SCENARIO 1] SECDED code derived from the Hamming bound for a 64-bit Word\n";
    const Code c = make_code(params().arch.word_bits);
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  data bits                  : " << c.data_bits << "\n";
    std::cout << "  Hamming parity (2^r>=k+r+1): " << c.hamming_parity << "  (2^" << c.hamming_parity
              << " = " << (1 << c.hamming_parity) << " >= " << (c.data_bits + c.hamming_parity + 1) << ")\n";
    std::cout << "  overall parity (double det): " << c.overall_parity << "\n";
    std::cout << "  total parity bits          : " << c.total_parity << "\n";
    std::cout << "  codeword length            : " << c.codeword << " bits\n";
    std::cout << "  overhead ratio             : " << c.overhead_ratio << "x\n";
    require(c.hamming_parity == 7, "a 64-bit Word needs 7 Hamming parity bits");
    require(c.codeword == 72, "SECDED on 64 bits gives a 72-bit codeword");
    std::cout << "  capability: corrects exactly 1 bit, detects 2, cannot correct 2.\n";
}

static void scenario_v2_overhead_claim() {
    std::cout << "\n[SCENARIO 2] V2's ~2-3x SECDED overhead does not match a 64-bit Word\n";
    const Code c64 = make_code(params().arch.word_bits);
    const Code c8 = make_code(8);
    const Code c4 = make_code(4);
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  Word width   codeword   overhead\n";
    std::cout << "  " << std::setw(6) << c64.data_bits << "       " << std::setw(6) << c64.codeword
              << "      " << c64.overhead_ratio << "x\n";
    std::cout << "  " << std::setw(6) << c8.data_bits << "       " << std::setw(6) << c8.codeword
              << "      " << c8.overhead_ratio << "x\n";
    std::cout << "  " << std::setw(6) << c4.data_bits << "       " << std::setw(6) << c4.codeword
              << "      " << c4.overhead_ratio << "x\n";
    require(c64.overhead_ratio < 2.0,
            "SECDED on a 64-bit Word must come in under V2's claimed 2x lower bound");
    require(c64.overhead_ratio > 1.0, "any code adds overhead above 1x");
    std::cout << "  V2 states ~2-3x overhead at the Word level, but its Words are 64-bit.\n";
    std::cout << "  the derived cost is " << std::setprecision(3) << c64.overhead_ratio
              << "x. A factor of 2 only appears at a " << "4-bit Word.\n";
    std::cout << "  finding: either V2 means a different code, a different Word width, or the\n";
    std::cout << "  factor is wrong: it was asserted, not simulated.\n";
}

static void scenario_charged_into_budget() {
    std::cout << "\n[SCENARIO 3] SECDED overhead must be charged to capacity and area\n";
    const Code c = make_code(params().arch.word_bits);
    const double raw_bits = fea::payload_bits();
    const double protected_bits = raw_bits * c.overhead_ratio;
    const double added_bits = protected_bits - raw_bits;
    const double payload_capacity_TB = raw_bits / 8.0 / 1e12;
    const double physical_capacity_TB = protected_bits / 8.0 / 1e12;
    const double area_growth = c.overhead_ratio;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  payload bits (unchanged)   : " << (raw_bits / 1e12) << "e12\n";
    std::cout << "  physical bits with parity  : " << (protected_bits / 1e12) << "e12\n";
    std::cout << "  parity bits added          : " << (added_bits / 1e12) << "e12\n";
    std::cout << "  payload capacity           : " << payload_capacity_TB << " TB\n";
    std::cout << "  physical capacity          : " << physical_capacity_TB << " TB\n";
    std::cout << "  physical bits occupy       : " << std::setprecision(4) << area_growth
              << "x the cell area\n";
    require(added_bits > 0.0, "parity must add physical bits");
    require(physical_capacity_TB > payload_capacity_TB,
            "parity must raise physical capacity while payload capacity stays fixed");
    require(area_growth < 1.2, "SECDED overhead must stay under 20% for a 64-bit Word");
    std::cout << "  payload does NOT grow: parity buys correction, not storage. Physical\n";
    std::cout << "  occupancy does grow, and that is what must be charged to the floorplan.\n";
    std::cout << "  at V2's practical density the cells already fill ~99% of the die (M2).\n";
    std::cout << "  adding " << std::setprecision(2) << ((area_growth - 1.0) * 100.0)
              << "% more physical bits makes the floorplan worse, not neutral.\n";
    std::cout << "  SECDED area must be added inside the FZC Block budget, not assumed free.\n";
}

static void scenario_power_latency_overhead() {
    std::cout << "\n[SCENARIO 4] SECDED adds parity generation, syndrome, and correction latency\n";
    const Code c = make_code(params().arch.word_bits);
    // Declared, not measured: syndrome XOR cost scales with parity bits.
    const int syndrome_xors = c.total_parity;
    const double per_bit_energy_J = params().control.decoder_event_J / params().control.decoder_transistors;
    const double parity_energy_J = syndrome_xors * per_bit_energy_J;
    const double cycles_added = 1.0; // one syndrome+correct pass, declared

    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  parity bits in syndrome     : " << syndrome_xors << "\n";
    std::cout << "  energy per parity XOR       : " << per_bit_energy_J << " J (declared)\n";
    std::cout << "  syndrome energy per Word    : " << parity_energy_J << " J\n";
    std::cout << std::fixed << std::setprecision(0);
    std::cout << "  added latency per Word      : " << cycles_added << " cycle\n";
    // cycles_added is DECLARED 1.0, so `require(cycles_added > 0.0)` gated a
    // constant and could never fail. That gate is removed; the declared nature
    // of the latency is stated in the output instead. What is computed here is
    // the parity work, which scales with the code parameters derived in SCENARIO 1.
    require(parity_energy_J > 0.0, "parity work must cost energy");
    require(syndrome_xors > 0, "the derived code must have parity bits to compute");
    std::cout << "  label: the per-XOR energy is a declared scale, not a measured CMOS figure.\n";
    std::cout << "  what matters structurally: correction is on the critical path, so it enters\n";
    std::cout << "  M9 cycle time, and it consumes energy, so it enters M1.\n";
}

static void scenario_burst_limit() {
    std::cout << "\n[SCENARIO 5] SECDED exercised for real: single corrected, double detected, bursts unsafe\n";
    const uint64_t original = 0x0123456789ABCDEFULL;
    const int ham = 71;
    const std::vector<int> base = encode(original);
    uint64_t decoded = 0;
    const DecodeResult clean = decode(base, decoded);
    require(clean == kClean && decoded == original,
            "a clean codeword must decode to the original data with no correction");

    // (a) every single-bit position must be corrected back to the original.
    int single_ok = 0;
    const int single_total = ham + 1;
    for (int pos = 0; pos < single_total; ++pos) {
        std::vector<int> cw = base;
        cw[pos] ^= 1;
        const DecodeResult r = decode(cw, decoded);
        if (r == kCorrected && decoded == original) ++single_ok;
    }
    require(single_ok == single_total,
            "every single-bit error in the codeword must be corrected to the original data");

    // (b) every pair must be detected and NOT corrected.
    int pair_detected = 0;
    int pair_wrongly_corrected = 0;
    int pair_total = 0;
    for (int i = 0; i < single_total; ++i) {
        for (int j = i + 1; j < single_total; ++j) {
            ++pair_total;
            std::vector<int> cw = base;
            cw[i] ^= 1;
            cw[j] ^= 1;
            const DecodeResult r = decode(cw, decoded);
            if (r == kDetectedUncorrectable) ++pair_detected;
            if (r == kCorrected) ++pair_wrongly_corrected;
        }
    }
    require(pair_detected == pair_total,
            "every double-bit error must be detected and refused, not corrected");
    require(pair_wrongly_corrected == 0, "no double-bit error may ever be corrected");

    // (c) three-bit bursts: count how often the decoder silently miscorrects.
    int triple_miscorrect = 0;
    int triple_detected = 0;
    int triple_lucky = 0;
    int triple_total = 0;
    for (int i = 0; i < single_total; ++i) {
        for (int j = i + 1; j < single_total; ++j) {
            for (int k = j + 1; k < single_total; ++k) {
                ++triple_total;
                std::vector<int> cw = base;
                cw[i] ^= 1; cw[j] ^= 1; cw[k] ^= 1;
                const DecodeResult r = decode(cw, decoded);
                if (r == kDetectedUncorrectable) ++triple_detected;
                else if (decoded == original) ++triple_lucky;
                else ++triple_miscorrect;
            }
        }
    }

    std::cout << std::fixed << std::setprecision(0);
    std::cout << "  clean codeword             : " << (clean == kClean ? "decoded clean" : "not clean") << "\n";
    std::cout << "  single-bit positions tested: " << single_total
              << ", corrected to original: " << single_ok << "\n";
    std::cout << "  double-bit pairs tested    : " << pair_total
              << ", detected: " << pair_detected
              << ", wrongly corrected: " << pair_wrongly_corrected << "\n";
    std::cout << "  three-bit bursts tested    : " << triple_total << "\n";
    std::cout << "    detected, not corrected  : " << triple_detected << "\n";
    std::cout << "    miscorrected silently    : " << triple_miscorrect << "  <- UNSAFE\n";
    std::cout << "    landed back on original  : " << triple_lucky << "\n";

    require(triple_miscorrect > 0,
            "some three-bit burst must silently miscorrect, which is the hazard being shown");
    require(triple_miscorrect > triple_detected,
            "most three-bit bursts must miscorrect, since SECDED only classifies 1 and 2");

    const double miscorrect_rate = static_cast<double>(triple_miscorrect) / triple_total;
    std::cout << "  miscorrect rate at 3 bits  : " << std::setprecision(4) << miscorrect_rate
              << " of " << triple_total << " combinations\n";
    std::cout << "  finding: SECDED is verified correct for 1 bit and verified refusing for 2,\n";
    std::cout << "  but a 3+ bit common-mode burst is silently miscorrected about "
              << std::setprecision(1) << (miscorrect_rate * 100.0) << "% of the time.\n";
    std::cout << "  M7's 5% common-mode ceiling is a fault class SECDED cannot cover, because\n";
    std::cout << "  shared causes produce bursts rather than isolated flips.\n";
    std::cout << "  label: codec properties computed by encoding and decoding, not asserted.\n";
}

} // namespace secded

int main() {
    using namespace secded;
    try {
        std::cout << "FEA V3 M8 SECDED overhead and limits\n";
        std::cout << "Code derived from the Hamming bound for the stated Word width.\n";
        scenario_derive_the_code();
        scenario_v2_overhead_claim();
        scenario_charged_into_budget();
        scenario_power_latency_overhead();
        scenario_burst_limit();
        std::cout << "\nPASS: SECDED cost derived, V2's factor challenged, burst limit exposed.\n";
        std::cout << "NEXT EVIDENCE GATE: measure the actual fault class (independent vs correlated) so code choice follows evidence.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
