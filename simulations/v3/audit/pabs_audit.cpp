// =============================================================================
// simulations/v3/audit/pabs_audit.cpp
//
// AUDIT TOOL. This is NOT one of the 27 gated V3 targets and is not part of
// `make check`. It reproduces three quantities the manuscript states, using the
// same Crank-Nicolson scheme, units and geometry as the retained reference
// propagator in simulations/v2/FEA_sim_v2.cpp (SIM 4):
//
//   PART 0  the committed on-resonance value: SIM 4 prints "Absorbed (cluster):
//           0.4608" for the 500-site chain, cluster at site 250, packet centre
//           x0 = 100, sigma = 20, 2000 steps. Reproduced here.
//   PART A  packet-width sweep on a longer chain (sigma = 5..80 sites). The
//           absorbed fraction converges to the time-independent single-site
//           scattering coefficient rather than drifting with the pulse, which is
//           why the manuscript states the 0.4608 is not a finite-pulse artifact.
//   PART B  the on/off contrast is convention dependent. The gate factor in
//           SIM 4 is 1/(1 + (E_ctr/GAMMA_DERIVED_meV)^2), i.e. scaled by Gamma:
//           that gives the committed 1,066x. Scaling it by Gamma/2, the
//           half-width convention used for the line shape elsewhere in the
//           manuscript, gives 4,194x. The on-resonance value is unaffected
//           because the gate factor is 1 at E_ctr = 0.
//
// Closed form for a one-site loss channel in a 1D lead on resonance (k = pi/2):
//   eta = Gamma/(2t),  T = 4/(2+eta)^2,  R = eta^2/(2+eta)^2,  A = 1 - T - R
//   at Gamma = 45 meV, t = 20 meV: eta = 1.125, T = 0.4096, R = 0.1296, A = 0.4608
//
// Build and run:  c++ -O2 -std=c++17 -o pabs_audit pabs_audit.cpp && ./pabs_audit
// =============================================================================

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

using cd = std::complex<double>;

namespace {
constexpr double kTtMeV = 20.0;   // lead hopping, matches fea_params t_hop_eV
constexpr double kEta = 1.125;    // Gamma/(2t) = 45/40
constexpr double kDt = 0.1;
const double kKF = M_PI / 2.0;

// One Crank-Nicolson run. Returns the probability lost to the imaginary sink.
double absorbed(int N, int cluster_site, double x0, double sig, int steps,
                double E_ctr_meV, double gate_scale_meV) {
    const double gate =
        1.0 / (1.0 + (E_ctr_meV / gate_scale_meV) * (E_ctr_meV / gate_scale_meV));
    std::vector<cd> psi(N);
    double norm = 0.0;
    for (int i = 0; i < N; ++i) {
        double env = std::exp(-((i - x0) * (i - x0)) / (4.0 * sig * sig));
        psi[i] = cd(env * std::cos(kKF * i), env * std::sin(kKF * i));
        norm += std::norm(psi[i]);
    }
    for (int i = 0; i < N; ++i) psi[i] /= std::sqrt(norm);

    std::vector<cd> V(N, cd(0, 0));
    V[cluster_site] = cd(E_ctr_meV / kTtMeV, -kEta * gate);

    const cd beta(0.0, kDt * 0.5);
    for (int s = 0; s < steps; ++s) {
        std::vector<cd> rhs(N);
        for (int i = 0; i < N; ++i) {
            cd l = (i == 0) ? cd(0, 0) : psi[i - 1];
            cd r = (i == N - 1) ? cd(0, 0) : psi[i + 1];
            rhs[i] = psi[i] + beta * (l + r) - beta * V[i] * psi[i];
        }
        const cd nb = -beta;
        std::vector<cd> cp(N), dp(N), diag(N);
        for (int i = 0; i < N; ++i) diag[i] = cd(1, 0) + beta * V[i];
        cp[0] = nb / diag[0];
        dp[0] = rhs[0] / diag[0];
        for (int i = 1; i < N; ++i) {
            cd den = diag[i] - nb * cp[i - 1];
            cp[i] = nb / den;
            dp[i] = (rhs[i] - nb * dp[i - 1]) / den;
        }
        psi[N - 1] = dp[N - 1];
        for (int i = N - 2; i >= 0; --i) psi[i] = dp[i] - cp[i] * psi[i + 1];
    }
    double tot = 0.0;
    for (int i = 0; i < N; ++i) tot += std::norm(psi[i]);
    return 1.0 - tot;
}
}  // namespace

int main() {
    const double T = 4.0 / ((2.0 + kEta) * (2.0 + kEta));
    const double R = kEta * kEta / ((2.0 + kEta) * (2.0 + kEta));
    const double A = 1.0 - T - R;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "FEA V3 capture audit (not a gated target)\n";
    std::cout << "closed form  eta=" << kEta << "  T=" << T << "  R=" << R
              << "  A=" << A << "\n\n";

    std::cout << "[PART 0] committed SIM 4 geometry (N=500, site 250, x0=100, sigma=20, 2000 steps)\n";
    const double on_committed = absorbed(500, 250, 100.0, 20.0, 2000, 0.0, 45.0);
    std::cout << "   on-resonance absorbed = " << on_committed
              << "   (SIM 4 committed output prints 0.4608)\n\n";

    std::cout << "[PART A] packet-width sweep, longer chain (N=900, site 450, x0=100, 4000 steps)\n";
    for (double sig : {5.0, 10.0, 20.0, 40.0, 80.0}) {
        const double a = absorbed(900, 450, 100.0, sig, 4000, 0.0, 45.0);
        std::cout << "   sigma=" << std::setw(5) << sig << " sites   A=" << a
                  << "   |A-closed|=" << std::fabs(a - A) << "\n";
    }

    std::cout << "\n[PART B] on/off contrast vs the gate-factor scale (committed geometry)\n";
    const double off_gamma = absorbed(500, 250, 100.0, 20.0, 2000, 300.0, 45.0);
    const double off_half = absorbed(500, 250, 100.0, 20.0, 2000, 300.0, 22.5);
    std::cout << "   gate scale Gamma  = 45.0 meV  off=" << off_gamma
              << "  on/off=" << on_committed / off_gamma << "\n";
    std::cout << "   gate scale Gamma/2= 22.5 meV  off=" << off_half
              << "  on/off=" << on_committed / off_half << "\n";
    std::cout << "   the on-resonance value is unaffected: the gate factor is 1 at E_ctr = 0.\n";
    std::cout << "\nAUDIT LABEL: reproduces cited quantities only; no physical claim beyond the model.\n";
    return 0;
}
