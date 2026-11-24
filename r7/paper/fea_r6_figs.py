"""Figures for FEA-R6."""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

BLUE, ORANGE, GREY, RED, GREEN, PURPLE = '#2b6cb0', '#dd6b20', '#718096', '#c53030', '#2f855a', '#6b46c1'
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9.5, 'axes.labelsize': 9, 'legend.fontsize': 7.5,
                     'figure.dpi': 150, 'savefig.bbox': 'tight',
                     'axes.spines.top': False, 'axes.spines.right': False})


def fig_layout(N, G):
    d, rp, c = G['d_cell'], G['r_port'], G['c']
    px = N['pitch'][0]
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.8), gridspec_kw=dict(width_ratios=[1.1, 1.1, 1.0]))
    x0 = ax[0]
    x0.plot([-3.6, 3.6], [0, 0], color=PURPLE, lw=3)
    x0.text(-3.6, 0.15, 'wire (n-doped dimer line)', fontsize=7, color=PURPLE)
    for k in (0,):
        x = k*px
        for s in (+1, -1):
            x0.add_patch(Circle((x, s*d), 0.22, color=ORANGE if (s > 0) == (k % 2 == 0) else BLUE))
            yd = s*(d + rp)
            x0.plot([x - 0.25, x + 0.25], [yd, yd], color=GREY, lw=1)
            x0.plot([x, x], [s*(d + rp - c), s*(d + rp)], color=GREEN, lw=2.2)
            x0.add_patch(Circle((x, s*(d + rp - c) if (s > 0) == (k % 2 == 0) else s*(d + rp)), 0.09, color=RED, zorder=5))
    x0.text(0.4, d - 0.05, 'Block, rail +', fontsize=7)
    x0.text(0.4, -d - 0.05, 'Block, rail −', fontsize=7)
    x0.text(0.4, d + rp - 0.1, 'H-memory dimer', fontsize=7, color=GREEN)
    x0.annotate('', xy=(-1.5, d), xytext=(-1.5, 0), arrowprops=dict(arrowstyle='<->', color=GREY))
    x0.text(-1.4, d/2, f'{d:.2f} nm', fontsize=7)
    x0.set_xlim(-3.7, 3.7); x0.set_ylim(-d - rp - 0.6, d + rp + 0.7); x0.axis('off'); x0.set_aspect('equal')
    x0.set_title('(a) One bit: mirrored pairs across the wire', loc='left')
    x0.text(-3.6, -d - rp - 0.55, 'red: dangling bond of the memory dimer\n(H on the other atom); next bit 8 nm along the wire', fontsize=6.3, color=GREY)
    # (b) tiers
    x1 = ax[1]
    tiers = [('electron (working)', N['working_tau_s'][0], N['working_tau_s'][1], ORANGE),
             ('H within dimer (memory)', N['tau_dimer_range_300'][0], N['tau_dimer_range_300'][1], GREEN),
             ('H along row (archive)', N['tau_row_range_300'][0], N['tau_row_range_300'][1], BLUE)]
    for i, (lab, lo, hi, col) in enumerate(tiers):
        x1.barh(i, math.log10(hi) - math.log10(lo), left=math.log10(lo), color=col, height=0.5)
        x1.text(math.log10(hi) + 0.3, i, lab, va='center', fontsize=7.5)
    for val, lab in ((math.log10(60), 'min'), (math.log10(86400), 'day'), (math.log10(YEAR := 3.156e7), 'year'), (math.log10(YEAR*1e6), 'Myr')):
        x1.axvline(val, color=GREY, lw=0.5, ls=':'); x1.text(val, 2.55, lab, fontsize=6.5, ha='center', color=GREY)
    x1.set_yticks([]); x1.set_xlim(-5, 24); x1.set_ylim(-0.6, 2.8)
    x1.set_xlabel(r'$\log_{10}$ mean lifetime at 300 K (s)')
    x1.set_title('(b) Three storage tiers at 300 K', loc='left')
    # (c) read contrast and error vs port distance
    x2 = ax[2]
    r = np.linspace(0.6, 2.0, 100)
    cont = [G['contrast'](x)*1e3 for x in r]
    x2.semilogy(r, cont, color=GREEN, label='read contrast (meV)')
    err = [G['read_error'](x) for x in r]
    x2.semilogy(r, err, color=RED, label='single-read error')
    x2.axhline(N['read_err_majority3'], color=RED, ls=':', lw=0.8); x2.text(1.25, N['read_err_majority3']*1.8, 'majority of 3 at 0.77 nm', fontsize=7, color=RED)
    x2.axvline(G['r_port'], color='k', lw=0.6, ls=':')
    x2.set_xlabel('port Block to memory dimer (nm)'); x2.legend(frameon=False, loc='center right')
    x2.set_title('(c) Reading a hydrogen bit', loc='left')
    fig.tight_layout(); fig.savefig('figures/r6_fig1_hybrid.pdf'); plt.close(fig)


def make_all(N, G):
    fig_layout(N, G)
