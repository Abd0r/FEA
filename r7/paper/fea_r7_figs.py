#!/usr/bin/env python3
"""Figures for FEA-R7. Reads numbers_r7.json (fea_r7.py) and the chip-simulation outputs in ../opt."""
import json, math, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'opt'))
N = json.load(open('numbers_r7.json'))
SIM = json.load(open('../opt/chip_sim_results.json'))
AD = json.load(open('../opt/chip_sim_adaptive_results.json'))
OUT = 'figures'

BLUE, ORANGE, GREY, RED, GREEN, PURPLE, TEAL = '#2b6cb0', '#dd6b20', '#718096', '#c53030', '#2f855a', '#6b46c1', '#2c7a7b'
QC = {'DB0': BLUE, 'DB-': RED, 'DB+': ORANGE}
QL = {'DB0': r'DB$^0$', 'DB-': r'DB$^-$', 'DB+': r'DB$^+$'}
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9.5, 'axes.labelsize': 9, 'legend.fontsize': 7.5,
                     'figure.dpi': 150, 'savefig.bbox': 'tight',
                     'axes.spines.top': False, 'axes.spines.right': False})
def tag(ax, s): ax.text(-0.13, 1.04, s, transform=ax.transAxes, fontweight='bold', fontsize=11)

# ---------------------------------------------------------------- Fig 1  architecture
def fig1():
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4), gridspec_kw=dict(width_ratios=[0.85, 1.35, 1.0]))
    a = ax[0]; d = 1.536; ar, br = 0.384, 0.768
    for y in np.arange(-2.3, 2.31, br):                         # dimer rows
        a.plot([-2.2, 2.2], [y, y], color='#e2e8f0', lw=6, zorder=0)
    a.plot([-2.2, 2.2], [0, 0], color=PURPLE, lw=3.5, zorder=1)
    a.text(-2.05, 0.15, 'wire', color=PURPLE, fontsize=8)
    for s, filled in ((+1, True), (-1, False)):
        pts = [(0, s*d), (-ar, s*d), (ar, s*d), (0, s*(d + br)), (0, s*(d - br))]
        for (x, y) in pts:
            a.add_patch(Circle((x, y), 0.11, color=ORANGE if filled else 'white', ec=ORANGE, lw=1.2, zorder=3))
        a.text(0.62, s*d - 0.1, ('$s$ (rail +)' if s > 0 else r'$\bar s$ (rail $-$)'), fontsize=8)
    a.annotate('', xy=(-1.0, d), xytext=(-1.0, 0), arrowprops=dict(arrowstyle='<->', color=GREY, lw=0.8))
    a.text(-2.05, d/2, '1.54 nm', fontsize=7.5, color=GREY)
    a.set_xlim(-2.1, 2.3); a.set_ylim(-2.6, 2.6); a.set_aspect('equal'); a.axis('off')
    a.set_title('One dual-rail bit'); tag(a, 'a')

    b = ax[1]
    b.add_patch(Rectangle((0, 0), 10, 6, fill=False, ec=GREY, lw=1))
    for i in range(6):
        y = 0.6 + i*0.85
        b.plot([1.6, 9.6], [y, y], color=PURPLE, lw=1.6)
        for x in np.arange(2.2, 9.6, 0.75):
            b.add_patch(Circle((x, y + 0.2), 0.11, color=ORANGE, alpha=0.8)); b.add_patch(Circle((x, y - 0.2), 0.11, ec=ORANGE, fc='white'))
        b.add_patch(Circle((1.35, y), 0.17, color=GREEN)); b.add_patch(Circle((9.85, y), 0.17, color=TEAL))
    b.add_patch(Rectangle((0.2, 0.2), 0.75, 5.5, color=BLUE, alpha=0.25)); b.text(0.33, 2.2, 'FZC', rotation=90, fontsize=8, color=BLUE)
    b.text(1.0, 5.75, 'source', color=GREEN, fontsize=7); b.text(9.2, 5.75, 'sink', color=TEAL, fontsize=7)
    b.text(2.2, -0.6, '64-bit Words, bit-serial, lock-step', fontsize=7.5)
    b.set_xlim(-0.2, 10.5); b.set_ylim(-0.8, 6.4); b.set_aspect('equal'); b.axis('off')
    b.set_title('Zone: Word wires and FZC'); tag(b, 'b')

    c = ax[2]
    c.add_patch(Rectangle((0, 0), 10, 10, fc='#f7fafc', ec=GREY, lw=1.5))
    c.add_patch(Rectangle((0.25, 0.25), 9.5, 9.5, fill=False, ec=GREY, lw=0.6, ls=':'))
    kinds = [('memory + logic', ORANGE), ('in-memory AI', PURPLE), ('probabilistic (outlook)', GREEN)]
    rng = np.random.default_rng(3)
    for i in range(6):
        for j in range(6):
            k = 0 if (i + j) % 3 else (1 if i < 3 else 2)
            if k == 2 and rng.random() < 0.5: k = 0
            c.add_patch(Rectangle((0.7 + i*1.45, 0.7 + j*1.45), 1.25, 1.25, color=kinds[k][1], alpha=0.55))
    for k, (lab, col) in enumerate(kinds):
        c.add_patch(Rectangle((10.5, 7.6 - k*1.6), 0.6, 0.6, color=col, alpha=0.6)); c.text(11.3, 7.75 - k*1.6, lab, fontsize=7.5)
    c.text(10.5, 1.2, 'transistors only\nat the edge\n(power, I/O)', fontsize=7.5, color=GREY)
    c.set_xlim(-0.2, 17); c.set_ylim(-0.3, 10.3); c.set_aspect('equal'); c.axis('off')
    c.set_title('Chip: Zones by role'); tag(c, 'c')
    fig.tight_layout(w_pad=0.5); fig.savefig(f'{OUT}/r7_fig1_architecture.pdf'); plt.close(fig)

# ---------------------------------------------------------------- Fig 2  DFT: pattern stability
def fig2():
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.5))
    for q in ('DB0', 'DB-', 'DB+'):
        p = N['dft'][q]['profile']; x = np.linspace(0, 1, len(p))
        ax[0].plot(x, [r['dE'] for r in p], 'o-', color=QC[q], label=f"{QL[q]}  $E_a$ = {N['dft'][q]['Ea']:.2f} eV", ms=4)
        ax[1].plot(x, [r['dimer'] for r in p], 'o-', color=QC[q], label=QL[q], ms=4)
    ax[0].set_xlabel('reaction coordinate (H on atom $i$ → atom $k$)'); ax[0].set_ylabel('energy (eV)')
    ax[0].legend(frameon=False, loc='lower center', fontsize=6.8); ax[0].set_ylim(-0.45, 1.6); ax[0].set_title('Intra-dimer H hop, CI-NEB'); tag(ax[0], 'a')
    ax[1].set_xlabel('reaction coordinate'); ax[1].set_ylabel('Si–Si dimer length (Å)')
    ax[1].legend(frameon=False); ax[1].set_title(r'DB$^-$ breaks the dimer (bridged H)'); tag(ax[1], 'b')
    T = np.linspace(250, 450, 200); kB = 8.617e-5
    for q in ('DB0', 'DB-', 'DB+'):
        Ea = N['dft'][q]['Ea']; f = 0.5 if q == 'DB-' else 1.0
        lo = 1/(f*10**10.8*np.exp(-Ea/(kB*T)))/3.156e7; hi = 1/(f*10**9.8*np.exp(-Ea/(kB*T)))/3.156e7
        ph = 1/(f*1e13*np.exp(-Ea/(kB*T)))/3.156e7
        ax[2].fill_between(T - 273.15, lo, hi, color=QC[q], alpha=0.25)
        ax[2].plot(T - 273.15, ph, '--', color=QC[q], lw=0.9)
        ax[2].plot(T - 273.15, np.sqrt(lo*hi), color=QC[q], label=QL[q])
    ax[2].set_yscale('log'); ax[2].axvline(27, color=GREY, lw=0.6, ls=':')
    ax[2].set_xlabel('temperature (°C)'); ax[2].set_ylabel('mean time per DB site to one hop (years)')
    ax[2].set_title('Pattern stability (band: $10^{10.3\\pm0.5}$ Hz; dashed: $10^{13}$ Hz)'); ax[2].legend(frameon=False)
    tag(ax[2], 'c')
    fig.tight_layout(); fig.savefig(f'{OUT}/r7_fig2_dft.pdf'); plt.close(fig)

# ---------------------------------------------------------------- Fig 3  retention engineering
def fig3():
    import retention_model as RM
    import fea_chip_sim as S
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.5))
    labels = ['none', '1', '2', '3', '2+3', '1+2+3', '1+2+3+4']
    opts = {1: None, 2: dict(n_s=1e10), 3: dict(d_nm=10.0), 4: {}}
    for li, (lam, col) in enumerate(((1.0, BLUE), (1.5, GREEN))):
        vals = []
        for lab in labels:
            p = dict(dE=0.56, lam=0.0, n_s=1e18, d_nm=1.0, dq=0.0)
            combo = [] if lab == 'none' else [int(x) for x in lab.split('+')]
            for k in combo:
                p.update(dict(lam=lam) if k == 1 else opts[k])
            lv = [x/100 for x in range(26, 87)] if 4 in combo else [0.56]
            vals.append(max(RM.tau(**{**p, 'dE': x}) for x in lv))
        ax[0].bar(np.arange(len(labels)) + (li - 0.5)*0.38, vals, 0.38, color=col, label=f'$\\lambda$ = {lam} eV')
    ax[0].set_yscale('log'); ax[0].set_xticks(range(len(labels))); ax[0].set_xticklabels(labels, rotation=30, fontsize=7.5)
    ax[0].set_xlabel('optimisations applied'); ax[0].set_ylabel('bit lifetime at 300 K (s)')
    ax[0].legend(frameon=False); ax[0].set_title('No single fix suffices; together they do'); tag(ax[0], 'a')
    lams = np.linspace(round(N['lambda_eV']['minus'], 2), 2.2, 60)   # from the DFT value up; the classical form is unreliable below it
    for TC, col in ((27, BLUE), (50, ORANGE), (85, RED)):
        y = [S.tau_bit(TC + 273.15, l, S.best_level(l)) for l in lams]
        ax[1].plot(lams, y, color=col, label=f'{TC} °C')
    if 'lambda_eV' in N:
        v = N['lambda_eV']['minus']
        ax[1].axvline(v, color=RED, ls='--', lw=1.0); ax[1].annotate(f'DFT: {v:.2f} eV', xy=(v, 3e-3), xytext=(v + 0.12, 3e0), fontsize=7, color=RED, arrowprops=dict(arrowstyle='->', color=RED, lw=0.7))
        ax[1].axvspan(1.0, 1.5, color=GREEN, alpha=0.08); ax[1].text(1.03, 1e2, 'needed:\nlaptop | phone', fontsize=6.5, color=GREEN)
    ax[1].set_yscale('log'); ax[1].set_xlabel(r'relaxation energy $\lambda$ (eV)'); ax[1].set_ylabel('bit lifetime (s)')
    ax[1].legend(frameon=False); ax[1].set_title('All optimisations: lifetime vs $\\lambda$'); tag(ax[1], 'b')
    t = N['tilt']; eps = np.linspace(4, 12, 50); COUL = 1.43996
    mod = [COUL/(e*t['r_near_nm'])*math.exp(-t['r_near_nm']/5) - COUL/(e*t['r_far_nm'])*math.exp(-t['r_far_nm']/5) for e in eps]
    ax[2].plot(eps, np.array(mod)*1e3, color=PURPLE, label='screened point charges')
    ax[2].axhline(t['dE_DFT_charge_part']*1e3, color=RED, label='DFT, charge-dependent part')
    ax[2].axhline(t['dE_DFT']*1e3, color=RED, ls='--', lw=0.8, label='DFT, two DB$^-$ (total)')
    ax[2].axhline(t['dE_DFT_neutral']*1e3, color=GREY, ls=':', lw=0.8, label='DFT, two DB$^0$')
    ax[2].axvline(5.6, color=GREY, ls=':', lw=0.8); ax[2].text(5.7, ax[2].get_ylim()[1]*0.9 if ax[2].get_ylim()[1] > 0 else 100, r'$\varepsilon_r$ used', fontsize=6.5)
    ax[2].set_xlabel(r'$\varepsilon_r$ of the point-charge model'); ax[2].set_ylabel('bit-state tilt (meV)')
    ax[2].legend(frameon=False); ax[2].set_title('Neighbour charge: DFT vs crosstalk model'); tag(ax[2], 'c')
    fig.tight_layout(); fig.savefig(f'{OUT}/r7_fig3_retention.pdf'); plt.close(fig)

# ---------------------------------------------------------------- Fig 4  whole-chip simulation
def fig4():
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.6), gridspec_kw=dict(width_ratios=[1.5, 1, 1]))
    a = ax[0]; works = ['sleep', 'office/web', 'video call', 'compile', 'gaming', 'AI inference']
    L = f"{round(N['lambda_eV']['minus'], 3)}"; w = 0.38
    for k, (mode, col) in enumerate((('fixed', GREY), ('adaptive', GREEN))):
        rows = {r['work']: r for r in AD[f'{L}/desktop'][mode]['rows']}
        a.bar(np.arange(len(works)) + (k - 0.5)*w, [rows[x]['P_ref']*1e3 for x in works], w, color=col, label=f'{mode} refresh')
    a.set_xticks(range(len(works))); a.set_xticklabels(works, rotation=20, fontsize=7.5)
    a.set_ylabel('refresh power (mW)'); a.set_title(f'Desktop at the DFT $\\lambda$ = {float(L):.2f} eV')
    a.legend(frameon=False); tag(a, 'a')
    b = ax[1]; machines = ['phone', 'laptop', 'desktop']
    LD = f"{round(N['lambda_eV']['minus'], 3)}"
    for li, (lam, hatch, lab) in enumerate(((LD, 'xx', f'{float(LD):.2f} (DFT)'), ('1.0', '//', '1.0'), ('1.5', '', '1.5'))):
        for k, (mode, col) in enumerate((('fixed', GREY), ('adaptive', GREEN))):
            v = [AD[f'{lam}/{m}'][f'{mode}_max_ambient_gaming_C'] for m in machines]
            b.bar(np.arange(3) + (li*2 + k - 2.5)*0.14, v, 0.14, color=col, hatch=hatch, ec='k', lw=0.3,
                  label=f'{mode}, $\\lambda$={lab}')
    b.axhline(25, color=GREY, ls=':', lw=0.8); b.text(2.45, 26, 'room', fontsize=6.5, color=GREY)
    b.set_xticks(range(3)); b.set_xticklabels(machines); b.set_ylabel('highest ambient while gaming (°C)')
    b.legend(frameon=False, fontsize=6, ncol=3, loc='upper center', bbox_to_anchor=(0.5, -0.12)); b.set_title('Adaptive refresh: thermal margin'); tag(b, 'b')
    c = ax[2]; lamTF = [2.0, 5.0, 10.0]
    per = [N[f'raw_GB_lam{l}']*64/92 for l in lamTF]
    layers = np.array([1, 2, 4, 8, 16])
    for l, v, col in zip(lamTF, per, (GREEN, BLUE, ORANGE)):
        c.plot(layers, v*layers, 'o-', color=col, label=f'$\\lambda_{{TF}}$ = {l:g} nm')
        c.plot(layers, v*layers*6, 's:', color=col, ms=3)
    c.set_xscale('log', base=2); c.set_yscale('log'); c.set_xlabel('stacked FEA layers'); c.set_ylabel('usable working memory (GB)')
    c.legend(frameon=False); c.set_title('Capacity: 0.5 cm$^2$ (solid), 3 cm$^2$ (dotted)'); tag(c, 'c')
    fig.tight_layout(); fig.savefig(f'{OUT}/r7_fig4_chip.pdf'); plt.close(fig)

if __name__ == '__main__':
    for f in (fig1, fig2, fig3, fig4):
        f(); print('done', f.__name__)
