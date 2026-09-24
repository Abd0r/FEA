// =============================================================================
// FEA_fzc_recognition_v3.cpp -- parameterized FZC pattern-recognition model
//
// Scope: a dimensionless electrostatic matched-pattern recognizer for one FZC
// actuator. It tests codeword selectivity, rail corruption, neighbour offset,
// and static disorder. Couplings, threshold, and disorder are abstract model
// parameters, not calibrated Si dangling-bond device values or performance
// predictions. This model does not include transport, capture, relaxation,
// sensor physics, powered gain/restoration, timing, energy, or retention.
// =============================================================================

#include <array>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace recognition {

constexpr int kBits = 3;
constexpr int kRails = 2 * kBits;
using Codeword = std::array<int, kBits>;
using Rails = std::array<int, kRails>;

struct Decision {
    double score = 0.0;
    bool fires = false;
};

struct Campaign {
    int valid = 0;
    int invalid = 0;
    int missed = 0;
    int false_fires = 0;
};

class DualRailRecognizer {
public:
    DualRailRecognizer(Codeword target, double threshold, double neighbour_offset = 0.0)
        : target_(target), threshold_(threshold), neighbour_offset_(neighbour_offset) {}

    Decision decide(const Rails& rails, double disorder = 0.0) const {
        double score = neighbour_offset_ + disorder;
        for (int bit = 0; bit < kBits; ++bit) {
            const int selected = target_[bit] == 0 ? 0 : 1;
            const int other = 1 - selected;
            score += rails[2 * bit + selected] - rails[2 * bit + other];
        }
        return Decision{score, score > threshold_};
    }

    double threshold() const { return threshold_; }
    double neighbour_offset() const { return neighbour_offset_; }

private:
    Codeword target_;
    double threshold_;
    double neighbour_offset_;
};

static Rails encode(const Codeword& code) {
    Rails rails{};
    for (int bit = 0; bit < kBits; ++bit) rails[2 * bit + code[bit]] = 1;
    return rails;
}

static Rails flip_pair(const Rails& rails, int bit) {
    Rails result = rails;
    std::swap(result[2 * bit], result[2 * bit + 1]);
    return result;
}

static void require(bool ok, const std::string& what) {
    if (!ok) throw std::runtime_error("ASSERTION FAILED: " + what);
}

static Campaign campaign(const DualRailRecognizer& recognizer, const Codeword& target, const std::vector<Codeword>& codebook,
                         double disorder_sigma, int trials, uint32_t seed) {
    Campaign result{};
    std::vector<Codeword> invalid;
    for (const Codeword& code : codebook) if (code != target) invalid.push_back(code);
    if (invalid.empty()) throw std::runtime_error("campaign requires at least one invalid codeword");
    std::mt19937 rng(seed);
    std::normal_distribution<double> disorder(0.0, disorder_sigma);
    for (int i = 0; i < trials; ++i) {
        const bool is_valid = (i % 4) == 0;
        const Codeword& code = is_valid ? target : invalid[static_cast<std::size_t>(i % static_cast<int>(invalid.size()))];
        const Decision d = recognizer.decide(encode(code), disorder(rng));
        if (is_valid) {
            ++result.valid;
            if (!d.fires) ++result.missed;
        } else {
            ++result.invalid;
            if (d.fires) ++result.false_fires;
        }
    }
    return result;
}

static void ideal_selectivity() {
    std::cout << "\n[SCENARIO 1] ideal dual-rail codeword separation\n";
    const Codeword target{0, 1, 1};
    const std::vector<Codeword> codebook{{0, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0}};
    const DualRailRecognizer r(target, 1.0);
    const double correct = r.decide(encode(target)).score;
    double nearest_invalid = -1e9;
    for (const Codeword& code : codebook) {
        if (code == target) continue;
        nearest_invalid = std::max(nearest_invalid, r.decide(encode(code)).score);
        require(!r.decide(encode(code)).fires, "invalid codeword must not fire in the ideal model");
    }
    require(r.decide(encode(target)).fires, "target codeword must fire in the ideal model");
    require(correct == 3.0 && nearest_invalid == -1.0, "expected dimensionless matched-pattern separation");
    std::cout << "  target score=" << correct << ", nearest invalid=" << nearest_invalid
              << ", target threshold margin=" << (correct - r.threshold()) << "\n";
}

static void rail_fault_rejection() {
    std::cout << "\n[SCENARIO 2] one corrupted dual-rail pair is rejected\n";
    const Codeword target{0, 1, 1};
    const DualRailRecognizer r(target, 1.5);
    const double ideal = r.decide(encode(target)).score;
    const Decision d = r.decide(flip_pair(encode(target), 1));
    require(r.decide(encode(target)).fires, "uncorrupted target must clear a 1.5 threshold");
    require(d.score == ideal - 2.0, "one rail-pair flip must reduce the matched score by 2");
    require(d.score < r.threshold() && !d.fires, "the reduced score must fall below the threshold");
    std::cout << "  ideal=" << ideal << " corrupted=" << d.score << " threshold=" << r.threshold() << " -> rejected\n";
}

static void neighbour_sensitivity() {
    std::cout << "\n[SCENARIO 3] bounded neighbour offset erodes but does not erase margin\n";
    const Codeword target{0, 1, 1};
    const Codeword invalid{1, 0, 1};
    const DualRailRecognizer r(target, 1.0, 0.25);
    const Decision valid = r.decide(encode(target));
    const Decision wrong = r.decide(encode(invalid));
    require(valid.fires, "target must remain above threshold under declared neighbour offset");
    require(!wrong.fires, "nearest invalid must remain below threshold under declared neighbour offset");
    require(valid.score - r.threshold() == 2.25, "valid margin must include declared neighbour offset");
    std::cout << "  neighbour offset=" << r.neighbour_offset() << ", valid margin=" << (valid.score - r.threshold())
              << ", invalid headroom=" << (r.threshold() - wrong.score) << "\n";
}

static void seeded_disorder_falsification() {
    std::cout << "\n[SCENARIO 4] seeded disorder campaign exposes a selectivity failure regime\n";
    const Codeword target{0, 1, 1};
    const std::vector<Codeword> codebook{{0, 0, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0}};
    const DualRailRecognizer r(target, 1.0);
    const Campaign low = campaign(r, target, codebook, 0.10, 10000, 20260922);
    const Campaign high = campaign(r, target, codebook, 1.00, 10000, 20260922);
    require(low.missed == 0 && low.false_fires == 0, "low abstract disorder must preserve ideal separation in this seeded run");
    require(high.missed > 0 && high.false_fires > 0, "high abstract disorder must reveal both miss and false-fire modes");
    std::cout << "  sigma=0.10: missed=" << low.missed << "/" << low.valid << ", false=" << low.false_fires << "/" << low.invalid << "\n";
    std::cout << "  sigma=1.00: missed=" << high.missed << "/" << high.valid << ", false=" << high.false_fires << "/" << high.invalid << "\n";
}

} // namespace recognition

int main() {
    using namespace recognition;
    try {
        std::cout << "FEA V3 FZC dimensionless pattern-recognition simulator\n";
        std::cout << "Model values are abstract selectivity parameters, not device-physics predictions.\n";
        ideal_selectivity();
        rail_fault_rejection();
        neighbour_sensitivity();
        seeded_disorder_falsification();
        std::cout << "\nPASS: abstract recognition invariants held.\n";
        std::cout << "NEXT EVIDENCE GATE: calibrate coupling, disorder, neighbour offsets, clocked gain, and sensing against a physical device model.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
