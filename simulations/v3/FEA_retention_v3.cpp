// =============================================================================
// FEA_retention_v3.cpp -- M4 Kramers retention sensitivity
//
// V2 stated 52.2 ms at 300 K and 2.1 ms at 330 K. The 330 K value was recomputed the
// same formula and got about 5.3 ms at 330 K. This module recomputes retention
// from Ec, the attempt frequency, and temperature must reproduce the reference
// number, and must fail on V2's 2.1 ms. It then sweeps Ec, attempt frequency,
// and temperature so the uncertainty is an output rather than a footnote.
// =============================================================================

#include "fea_params.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace retention {

using fea::check_unit;
using fea::params;
using fea::require;
using fea::kramers_tau_s;

// Kramers comes from fea_params. A previous revision kept a verbatim local copy
// here even though the header exists specifically so M4 and M13 share one
// definition and retention cannot diverge from refresh. Peer review 4 caught it,
// and it was the one remaining hole in the one-definition claim.

static void scenario_reproduce_reference_values() {
    std::cout << "\n[SCENARIO 1] reproduce V2's 300 K figure and the reference's 330 K recomputation\n";
    const auto& d = params().device;
    const double tau_300 = kramers_tau_s(d.Ec_eV, d.phonon_attempt_Hz, 300.0);
    const double tau_330 = kramers_tau_s(d.Ec_eV, d.phonon_attempt_Hz, 330.0);
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Ec = " << d.Ec_eV << " eV, nu0 = " << std::setprecision(2)
              << (d.phonon_attempt_Hz / 1e12) << " THz\n";
    std::cout << std::setprecision(4);
    std::cout << "  tau(300 K) = " << (tau_300 * 1e3) << " ms   (V2 claimed 52.2 ms)\n";
    std::cout << "  tau(330 K) = " << (tau_330 * 1e3) << " ms   (reference value ~5.31 ms, V2 claimed 2.1 ms)\n";

    check_unit("retention at 300 K", tau_300, "s", 0.0522, 0.05);
    check_unit("retention at 330 K", tau_330, "s", 0.00531, 0.05);

    const bool v2_330_ok = std::abs(tau_330 - 0.0021) / 0.0021 < 0.20;
    require(!v2_330_ok, "V2's 2.1 ms at 330 K must not reproduce from V2's own parameters");
    std::cout << "  V2's 330 K figure is off by a factor of " << std::setprecision(1)
              << (tau_330 / 0.0021) << " against its own stated Ec and nu0.\n";
    std::cout << "  reference 5.31 ms reproduces; V2's 2.1 ms does not.\n";
}

static void scenario_attempt_frequency_uncertainty() {
    std::cout << "\n[SCENARIO 2] local-mode attempt frequency is the dominant uncertainty\n";
    const auto& d = params().device;
    const std::vector<double> factors{0.5, 1.0, 2.0, 3.0};
    std::cout << std::fixed << std::setprecision(3);
    double lo = 1e300, hi = 0.0;
    for (const double f : factors) {
        const double tau = kramers_tau_s(d.Ec_eV, d.phonon_attempt_Hz * f, 300.0);
        lo = std::min(lo, tau);
        hi = std::max(hi, tau);
        std::cout << "  nu0 x " << std::setprecision(1) << f << std::setprecision(3)
                  << " -> tau = " << (tau * 1e3) << " ms\n";
    }
    require(hi / lo > 4.0, "a 6x sweep in attempt frequency must spread retention by the same factor");
    std::cout << "  span " << std::setprecision(1) << (lo * 1e3) << " to " << (hi * 1e3)
              << " ms. V2 called this a 2-3x effect; it is multiplicative on tau.\n";
}

static void scenario_charging_energy_sweep() {
    std::cout << "\n[SCENARIO 3] retention versus charging energy at 300 K\n";
    const auto& d = params().device;
    const std::vector<double> ecs{0.30, 0.40, 0.50, 0.65, 0.70};
    std::cout << std::fixed << std::setprecision(1);
    bool saw_catastrophic = false;
    bool saw_design = false;
    for (const double ec : ecs) {
        const double tau = kramers_tau_s(ec, d.phonon_attempt_Hz, 300.0);
        std::cout << "  Ec = " << ec << " eV -> tau = " << std::setprecision(4) << (tau * 1e3)
                  << " ms\n";
        if (tau < 1e-4) saw_catastrophic = true;
        if (std::abs(ec - d.Ec_eV) < 1e-9) saw_design = std::abs(tau - 0.0522) < 0.005;
    }
    require(saw_catastrophic, "a low-Ec point must show catastrophic retention so the sensitivity is visible");
    require(saw_design, "the Ec = 0.65 eV design point must land near V2's 52.2 ms");
    std::cout << "  Ec=0.3 eV collapses retention by orders of magnitude. The design point is a cliff,\n";
    std::cout << "  not a plateau. V2's single-point claim hides that.\n";
}

static void scenario_refresh_overhead() {
    std::cout << "\n[SCENARIO 4] retention sets the refresh contract the system must pay\n";
    const auto& d = params().device;
    const double tau = kramers_tau_s(d.Ec_eV, d.phonon_attempt_Hz, d.temperature_K);
    // Shared refresh contract, not a local copy: fea_params takes V2's tau/2
    // policy at lines 373 and 722, and payload is the derived Zone count.
    const double refresh_interval_s = fea::refresh_interval_s();
    const double refreshes_per_second = 1.0 / refresh_interval_s;
    const double bits = fea::refresh_bits();
    const double bytes_per_s = bits / 8.0 * refreshes_per_second;
    require(refresh_interval_s > 0.0 && refreshes_per_second > 0.0, "refresh quantities must be positive");
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  tau = " << (tau * 1e3) << " ms -> refresh interval " << (refresh_interval_s * 1e3)
              << " ms -> " << std::setprecision(0) << refreshes_per_second << " passes/s\n";
    std::cout << "  moving " << (bits / 1e14) << "e14 bits (data + FZC) every " << std::setprecision(4) << (refresh_interval_s * 1e3)
              << " ms is " << std::setprecision(3) << (bytes_per_s / 1e9) << " GB/s of refresh traffic alone.\n";
    const double array_GBs = bits / 8.0 / 1e9;
    const double required_GBps = bytes_per_s / 1e9;
    require(required_GBps > 0.0, "refresh bandwidth must be positive");
    std::cout << "  array size is " << std::setprecision(1) << array_GBs << " GB.\n";
    std::cout << "  to clear it once per " << std::setprecision(4) << (refresh_interval_s * 1e3)
              << " ms the fabric must sustain " << std::setprecision(0) << required_GBps
              << " GB/s of refresh traffic alone.\n";
    std::cout << "  V2 called refresh overhead < 1e-4. That assumed a pass much cheaper than\n";
    std::cout << "  clearing " << std::setprecision(1) << array_GBs << " GB inside "
              << (refresh_interval_s * 1e3)
              << " ms. This traffic must be charged to M1 power\n";
    std::cout << "  and to M9/M10 bandwidth, and it shares Slingshot with real work.\n";
}

static void scenario_kramers_regime_caveat() {
    std::cout << "\n[SCENARIO 5] the Kramers prefactor is an assumption, not a measurement\n";
    const auto& d = params().device;
    // the attempt frequency comes from bulk Si optical phonons,
    // while the trapped state is a 5-atom surface cluster with local modes.
    const double tau_bulk = kramers_tau_s(d.Ec_eV, d.phonon_attempt_Hz, 300.0);
    // A faster attempt frequency means faster escape, so it shortens retention.
    const double tau_faster = kramers_tau_s(d.Ec_eV, d.phonon_attempt_Hz * 10.0, 300.0);
    const double tau_slower = kramers_tau_s(d.Ec_eV, d.phonon_attempt_Hz * 0.1, 300.0);
    require(tau_faster < tau_bulk && tau_bulk < tau_slower,
            "prefactor uncertainty must bracket the nominal value");
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  nominal nu0 (bulk optical phonon) : " << (tau_bulk * 1e3) << " ms\n";
    std::cout << "  nu0 x 10 (faster local mode)      : " << (tau_faster * 1e3) << " ms\n";
    std::cout << "  nu0 / 10 (slower local mode)      : " << (tau_slower * 1e3) << " ms\n";
    std::cout << "  No local-mode measurement exists for the 5-atom cluster, so retention at 300 K\n";
    std::cout << "  stays an open experimental question. This module bounds it; it does not close it.\n";
    std::cout << "  label: derived formula, estimated inputs, open physical validation.\n";
}

} // namespace retention

int main() {
    using namespace retention;
    try {
        std::cout << "FEA V3 M4 retention sensitivity\n";
        std::cout << "Kramers formula recomputed from Ec, attempt frequency, and temperature.\n";
        scenario_reproduce_reference_values();
        scenario_attempt_frequency_uncertainty();
        scenario_charging_energy_sweep();
        scenario_refresh_overhead();
        scenario_kramers_regime_caveat();
        std::cout << "\nPASS: reference figure reproduced, V2 figure rejected, sensitivity exposed.\n";
        std::cout << "NEXT EVIDENCE GATE: local-mode phonon spectrum or an ab initio/kMC stability study at 300 K.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
