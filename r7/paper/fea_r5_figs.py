"""Figures for FEA-R5. Called from fea_r5.py --figs."""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch

BLUE, ORANGE, GREY, RED, GREEN, PURPLE = '#2b6cb0', '#dd6b20', '#718096', '#c53030', '#2f855a', '#6b46c1'
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9.5, 'axes.labelsize': 9,
                     'legend.fontsize': 7.5, 'figure.dpi': 150, 'savefig.bbox': 'tight',
                     'axes.spines.top': False, 'axes.spines.right': False})


def fig1(N, G):
    a, b, c = G['a'], G['b'], G['c']
    fig = plt.figure(figsize=(11, 3.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.35, 1.0], wspace=0.25)
    # (a) cell on the lattice: every DB on a different dimer
    ax = fig.add_subplot(gs[0])
    for i in range(-3, 4):
        for j in range(-2, 2):
            ax.plot(i*a, j*b, 'o', color='#cbd5e0', ms=5, zorder=1)
            ax.plot(i*a, j*b + c, 'o', color='#cbd5e0', ms=5, zorder=1)
            ax.plot([i*a, i*a], [j*b, j*b + c], color='#e2e8f0', lw=1.5, zorder=0)
    sites = [(-a, 0), (a, 0), (0, -(b - c)), (0, b)]
    for x, y in sites:
        ax.plot([0, x], [0, y], color=GREY, lw=1, zorder=2)
        ax.add_patch(Circle((x, y), 0.07, color=BLUE, zorder=3))
    ax.add_patch(Circle((0, 0), 0.085, color=ORANGE, zorder=4))
    ax.text(-a/2, -0.14, '3.84 Å', ha='center', fontsize=7)
    ax.text(0.07, -(b - c)/2, '5.43 Å', fontsize=7)
    ax.text(0.07, b/2, '7.68 Å', fontsize=7)
    for i in range(-3, 4):
        ax.plot(i*a, 2*b - 0.05, 's', color=PURPLE, ms=4)
    ax.text(-1.2, 2*b + 0.07, 'DB wire (row +2)', fontsize=7, color=PURPLE)
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.75, 1.8); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title('(a) Five-DB cell, one DB per dimer', loc='left')
    ax.text(-1.25, -1.72, 'grey pairs: H-terminated dimers\norange: absorber  blue: cell sites', fontsize=6.5, color=GREY)
    # (b) tile
    ax = fig.add_subplot(gs[1])
    nr, nc = 5, 8
    for r in range(nr):
        y = r
        ax.plot([0.6, nc + 0.4], [y + 0.35, y + 0.35], color=PURPLE, lw=1.6)
        for k in range(nc):
            ax.add_patch(Rectangle((k + 0.75, y), 0.5, 0.25, color=ORANGE if (k + r) % 2 else BLUE))
    for k in range(nc):
        ax.plot([k + 1, k + 1], [-0.4, nr - 0.2], color=GREEN, lw=0.8, ls='--')
    for r in range(nr):
        ax.add_patch(Circle((0.1, r + 0.35), 0.17, color=RED))       # source Block
        ax.add_patch(Circle((-0.35, r + 0.35), 0.13, color=PURPLE))  # control Block
        ax.add_patch(Circle((nc + 0.85, r + 0.35), 0.17, color=GREEN))  # sink Block
    ax.text(-0.5, nr - 0.05, 'control  source', fontsize=6.3)
    ax.text(nc + 0.6, nr - 0.05, 'sink', fontsize=6.3)
    ax.text(nc/2 + 0.3, -0.75, 'dashed: buried Si:P column gate lines (hold / write / logic); all Blocks are Fusion Blocks', fontsize=6.3, ha='center', color=GREEN)
    ax.set_xlim(-0.7, nc + 1.4); ax.set_ylim(-1.0, nr + 0.3); ax.axis('off')
    ax.set_title('(b) Tile: one Word per wire, FZC Blocks at the ends', loc='left')
    # (c) energy diagram
    ax = fig.add_subplot(gs[2])
    ax.add_patch(Rectangle((0, -0.04), 1.0, 0.08, color=PURPLE, alpha=0.25))
    ax.text(0.03, 0.05, 'wire band (80 meV)', fontsize=7, color=PURPLE)
    ax.plot([1.4, 2.2], [0.0, 0.0], color=ORANGE, lw=2); ax.text(2.25, 0.03, 'capture level', fontsize=7, va='center')
    ax.plot([1.4, 2.2], [-0.07, -0.07], color=BLUE, lw=2); ax.text(2.25, -0.10, 'storage level', fontsize=7, va='center')
    ax.annotate('', xy=(1.8, -0.068), xytext=(1.8, -0.003), arrowprops=dict(arrowstyle='->', color=RED, lw=1))
    ax.text(0.95, -0.035, 'relaxation', fontsize=6.5, color=RED)
    ax.plot([1.4, 2.2], [-1.0, -1.0], color=BLUE, lw=2, ls='--'); ax.text(2.25, -1.0, 'storage level\nin hold', fontsize=7, va='center')
    ax.annotate('', xy=(0.5, -1.0), xytext=(0.5, -0.045), arrowprops=dict(arrowstyle='<->', color=GREY))
    ax.text(0.55, -0.55, r'$E_B\geq$ 0.98–1.08 eV' + '\n(R1)', fontsize=7)
    ax.set_xlim(0, 3.2); ax.set_ylim(-1.15, 0.15); ax.axis('off')
    ax.set_title('(c) Levels during write and hold', loc='left')
    fig.savefig('figures/r5_fig1_architecture.pdf'); plt.close(fig)


def fig2(N, G):
    scat, kappa = G['scat'], G['kappa']
    fig, ax = plt.subplots(2, 2, figsize=(9, 6.4))
    g = np.logspace(-8, 1.5, 400)*kappa
    Tt, Rr, Aa = scat(0.0, 0.0, g)
    x = g/kappa
    ax[0, 0].loglog(x, Aa, color=ORANGE, label='captured $A$')
    ax[0, 0].loglog(x, Tt, color=BLUE, label='transmitted $T$')
    ax[0, 0].loglog(x, Rr, color=GREY, label='reflected $R$')
    ax[0, 0].axvline(1, color='k', lw=0.6, ls=':'); ax[0, 0].text(1.3, 1e-9, 'write\n(critical)', fontsize=7)
    ax[0, 0].axvline(N['logic_gammaA_over_kappa'], color='k', lw=0.6, ls=':')
    ax[0, 0].text(N['logic_gammaA_over_kappa']*1.4, 1e-9, 'logic', fontsize=7)
    ax[0, 0].set_ylim(1e-14, 2); ax[0, 0].set_xlabel(r'relaxation rate $\gamma/\kappa$')
    ax[0, 0].set_ylabel('probability'); ax[0, 0].legend(frameon=False, loc='lower right')
    ax[0, 0].set_title('(a) Empty cell on resonance', loc='left')
    t = G['t']
    E = np.linspace(-2*t + 1e-6, 2*t - 1e-6, 2001)
    _, _, Ae = scat(E, 0.0, kappa)
    ax[0, 1].plot(E*1e3, Ae, color=ORANGE, label='capture $A(E)$, write mode')
    kT = G['kT']
    w = 1/np.cosh(E/(2*kT))**2; ax[0, 1].fill_between(E*1e3, 0, 0.5*w/w.max(), color=BLUE, alpha=0.15, label='300 K contact (band-limited)')
    wg = np.exp(-E**2/(2*0.001**2)); ax[0, 1].fill_between(E*1e3, 0, 0.5*wg/wg.max(), color=GREEN, alpha=0.3, label='1 meV injector')
    ax[0, 1].set_xlabel(r'$E$ relative to band centre (meV)'); ax[0, 1].set_ylabel('probability')
    ax[0, 1].set_xlim(-40, 40); ax[0, 1].set_ylim(0, 0.62); ax[0, 1].legend(frameon=False, loc='upper right')
    ax[0, 1].set_title('(b) Capture across the wire band', loc='left')
    D = np.logspace(-1.5, 0.3, 200)
    loss = [sum(scat(0.0, d, N['logic_gammaA_over_kappa']*kappa)[1:]) for d in D]
    ax[1, 0].loglog(D, loss, color=PURPLE, label='one hold-mode cell')
    ax[1, 0].loglog(D, 1 - (1 - np.array(loss))**127, color=RED, label='127 hold-mode cells')
    ax[1, 0].axhline(1e-3, color=GREY, ls='--', lw=0.8); ax[1, 0].axvline(1.0, color=GREY, ls=':', lw=0.8)
    ax[1, 0].set_xlabel('detuning of non-addressed cells (eV)'); ax[1, 0].set_ylabel('scattering loss')
    ax[1, 0].legend(frameon=False); ax[1, 0].set_title('(c) Passing the rest of the Word', loc='left')
    W = np.logspace(-11, -6, 200)
    kT = G['kT']
    for td, col in ((0.1e-12, BLUE), (1e-12, PURPLE)):
        ax[1, 1].semilogx(W*1e12, kT*np.log(W/(1e-9*td)), color=col, label=f'write relaxation {td*1e12:g} ps')
    ax[1, 1].axhline(0.0645, color=RED, ls='--', lw=0.8); ax[1, 1].text(12, 0.09, 'one optical phonon', fontsize=7, color=RED)
    ax[1, 1].axvline(100, color='k', lw=0.6, ls=':')
    ax[1, 1].set_xlabel('time a column stays addressed (ps)'); ax[1, 1].set_ylabel('minimum storage depth $\\Delta_s$ (eV)')
    ax[1, 1].set_ylim(0, 1.0); ax[1, 1].legend(frameon=False, loc='lower right')
    ax[1, 1].set_title('(d) Depth that protects stored ones while addressed', loc='left')
    fig.tight_layout(); fig.savefig('figures/r5_fig2_capture.pdf'); plt.close(fig)


def fig3(N, G):
    fig, ax = plt.subplots(figsize=(9.5, 3.4))
    cells = ['$a^+$', '$a^-$', '$b^+$', '$b^-$', '$m^-$', '$m^+$']
    xs = np.arange(len(cells))*1.4 + 1.2
    ax.plot([0, xs[-1] + 1.2], [2.6, 2.6], color=PURPLE, lw=2)
    ax.text(0, 2.75, 'source Block', fontsize=7); ax.text(xs[-1] + 0.3, 2.75, 'sink Block', fontsize=7)
    for x, lab in zip(xs, cells):
        ax.add_patch(Rectangle((x - 0.3, 2.0), 0.6, 0.35, fc='#fefcbf', ec=GREY))
        ax.text(x, 2.17, lab, ha='center', va='center', fontsize=8)
    steps = [
        ('1', r'$a^+,b^+$ logic;  $m^-$ write', r'electron reaches $m^-$ only if $a\wedge b$:  $m^-\leftarrow a\wedge b$'),
        ('2', r'$a^-$ logic;  $m^+$ write', r'$m^+\leftarrow \bar a$'),
        ('3', r'$b^-$ logic;  $m^+$ write', r'$m^+\leftarrow \bar a\vee\bar b = \overline{a\wedge b}$'),
    ]
    for k, (n, setup, res) in enumerate(steps):
        y = 1.35 - 0.55*k
        ax.text(0, y, f'step {n}', fontsize=8, weight='bold')
        ax.text(1.6, y, setup, fontsize=8)
        ax.text(5.0, y, res, fontsize=8)
    ax.text(0, -0.5, r'All other columns stay in hold and are transparent. $m^\pm$ start reset. Result: $m=\mathrm{NAND}(a,b)$ in dual rail,'
                     '\nwritten in place in three write steps (each repeated until verified); no read-out to the tile edge.', fontsize=7.5, color=GREY)
    ax.set_xlim(-0.2, xs[-1] + 1.4); ax.set_ylim(-0.8, 3.0); ax.axis('off')
    fig.savefig('figures/r5_fig3_nand.pdf'); plt.close(fig)


def fig4(N, G):
    tau, kB = G['tau'], G['kB']
    bits = N['design']['logical_bits']
    Eb = np.linspace(0.5, 2.0, 300)
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
    ax[0].semilogy(Eb, [tau(x) for x in Eb], color=PURPLE, label='300 K')
    ax[0].semilogy(Eb, [tau(x, T=350) for x in Eb], color=PURPLE, ls='--', label='350 K')
    ax[0].axvspan(G['Si_gap'], 2.0, color=GREY, alpha=0.12)
    ax[0].text(1.14, 1e-6, 'above the\nSi band gap', fontsize=7, color=GREY)
    for tt, col in ((1, RED), (3, ORANGE), (5, GREEN)):
        e = N['retention'][tt]['Eb'][1]
        ax[0].axvline(e, color=col, lw=0.8)
        ax[0].text(e - 0.01, 1e14 / (10**(2*tt)), f"t={tt}" if tt > 1 else 'SECDED', fontsize=7, color=col, ha='right')
    ax[0].axvline(0.65, color='k', lw=0.6, ls=':'); ax[0].text(0.66, 1e10, '0.65 eV', fontsize=7)
    ax[0].set_xlabel('storage barrier $E_B$ (eV)'); ax[0].set_ylabel('mean hold time (s)')
    ax[0].set_ylim(1e-9, 1e20); ax[0].legend(frameon=False, loc='lower right')
    ax[0].set_title('(a) Required barrier by code strength', loc='left')
    grid = np.linspace(0.85, 1.2, 36)
    need = []
    for x in grid:
        tt, _ = G['min_bch'](x, bits, tmax=30)
        need.append(tt if tt else np.nan)
    ax[1].plot(grid, need, 'o-', color=GREEN, ms=3)
    ax[1].axvline(G['Si_gap'], color=GREY, lw=0.8, ls=':')
    ax[1].set_xlabel('storage barrier $E_B$ (eV)'); ax[1].set_ylabel('BCH errors corrected per 64-bit Word')
    ax[1].set_title('(b) Code needed for <1 failure per array-year', loc='left')
    fig.tight_layout(); fig.savefig('figures/r5_fig4_retention.pdf'); plt.close(fig)


def fig5(N, G):
    a, b = G['a'], G['b']
    p = np.linspace(2, 30, 60)
    fig, ax = plt.subplots(figsize=(5.6, 3.5))
    ax.loglog(p, [G['shift_monopole'](x, x) for x in p], color=RED, label='single-rail, worst case')
    ax.loglog(p, [G['shift_dualrail'](x, x) for x in p], color=BLUE, label='dual-rail, worst case')
    ax.axhline(N['crit_meV']/1e3, color=GREY, ls='--', lw=0.8)
    ax.text(2.1, N['crit_meV']/1e3*1.2, r'criterion $\Gamma/2$ = %.1f meV' % N['crit_meV'], fontsize=7)
    ax.set_xticks([2, 3, 5, 10, 20, 30]); ax.set_xticklabels(['2', '3', '5', '10', '20', '30'])
    ax.set_xlabel('square cell pitch (nm)'); ax.set_ylabel('level shift from neighbours (eV)')
    ax.legend(frameon=False); ax.set_title(r'Neighbour shift, $\varepsilon_r=5.6$, $\lambda_\mathrm{TF}=5$ nm', loc='left')
    fig.savefig('figures/r5_fig5_crosstalk.pdf'); plt.close(fig)


def make_all(N, G):
    fig1(N, G); fig2(N, G); fig3(N, G); fig4(N, G); fig5(N, G)
