#!/usr/bin/env python3
"""Write numbers_r7.tex: LaTeX macros for every number quoted in FEA-R7.tex.
Order: fea_r7.py -> ../opt/fea_chip_sim.py, fea_chip_sim_adaptive.py -> this script."""
import json, math, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
N = json.load(open('numbers_r7.json'))
SIM = json.load(open('../opt/chip_sim_results.json'))
AD = json.load(open('../opt/chip_sim_adaptive_results.json'))
M = {}

def sci(x, d=1):
    e = int(math.floor(math.log10(abs(x)))); m = x/10**e
    if round(m, d) >= 10: m /= 10; e += 1
    return f'{m:.{d}f}\\times10^{{{e}}}'
def f(x, d=2): return f'{x:.{d}f}'

# chip
P = N['pitch_lam5.0']
M['Pitch'] = f'{P[0]:.2f}\\times{P[1]:.1f}'
M['RawGB'] = f(N['raw_GB_lam5.0'], 1); M['UsableGB'] = f(N['chip']['usable_GB'], 1)
M['RawGBlo'] = f(N['raw_GB_lam10.0'], 0); M['RawGBhi'] = f(N['raw_GB_lam2.0'], 0)
M['RawBits'] = sci(N['chip']['raw_bits'], 2)
M['ReadSilent'] = sci(N['read_silent_err']); M['ReadDisturb'] = sci(N['read_disturb'])
M['ERefresh'] = f(N['E_refresh_bit_aJ']); M['EGate'] = f(N['E_gate_aJ']); M['ELine'] = f(N['E_line_per_wire_op_aJ'], 0)
M['Kappa'] = f(N['crit_meV'], 1)
# DFT
D = N['dft']
M['EaZero'] = f(D['DB0']['Ea']); M['EaMinus'] = f(D['DB-']['Ea']); M['EaPlus'] = f(D['DB+']['Ea'])
M['BridgeE'] = f(D['DB-']['bridge_dE']); M['DimerEnd'] = f(D['DB-']['dimer_endpoint']); M['DimerBridge'] = f(D['DB-']['bridge_dimer'])
M['DimerPlusSaddle'] = f(min(r['dimer'] for r in D['DB+']['profile']))
M['DzMinus'] = f(D['DB_atom_dz_on_charging']['DB-']); M['DzPlus'] = f(-D['DB_atom_dz_on_charging']['DB+'])
M['SpinMoment'] = f(D['spin_check']['magmom_min'], 2); M['SpinDiff'] = f(abs(D['spin_check']['Ea_spin'] - D['DB0']['Ea'])*1e3, 2)
M['NebStepsMinus'] = str(D['DB-']['steps'])
S = N['stability']
M['LifeZero'] = sci(S['DB0']['Hill_clean_meas']['mean_life_yr']); M['LifeZeroPh'] = sci(S['DB0']['phonon']['mean_life_yr'])
M['LifeMinus'] = sci(S['DB-']['Hill_clean_meas']['mean_life_yr']); M['LifeMinusPh'] = sci(S['DB-']['phonon']['mean_life_yr'])
M['LifePlus'] = f(S['DB+']['Hill_clean_meas']['mean_life_yr'], 1); M['LifePlusPhDays'] = f(S['DB+']['phonon']['mean_life_yr']*365.25, 0)
M['NDB'] = sci(N['n_DB_chip'])
M['DefectsYr'] = sci(N['pattern_defects_per_year_chip']['Hill_clean_meas'])
M['DefectsYrPh'] = sci(N['pattern_defects_per_year_chip']['phonon'])
M['DefFrac'] = sci(S['DB0']['Hill_clean_meas']['per_DB_year']); M['DefFracPh'] = sci(S['DB0']['phonon']['per_DB_year'])
pb = lambda q: 100*(1 - (1 - S['DB0'][q]['per_DB_year'])**(10*10))     # % of bits (10 DBs) with a moved DB in 10 yr
M['BitTenYr'] = f(pb('Hill_clean_meas'), 3); M['BitTenYrPh'] = f(pb('phonon'), 0)
M['FplusMax'] = sci(N['f_plus_max']['Hill_clean_meas'])
A = N['anneal_Tmax_C']
M['TannealA'] = f(A['1e-06']['Hill_clean_meas'], 0); M['TannealB'] = f(A['0.0001']['Hill_clean_meas'], 0)
M['TannealC'] = f(A['0.01']['Hill_clean_meas'], 0); M['TannealAPh'] = f(A['1e-06']['phonon'], 0)
T = N['tilt']
M['TiltDFT'] = f(T['dE_DFT']*1e3, 0); M['TiltModel'] = f(T['dE_model_eps5_6']*1e3, 0); M['TiltModelB'] = f(T['dE_model_eps6_35']*1e3, 0)
M['TiltRatio'] = f(T['bias_ratio_300K'], 0); M['TiltRnear'] = f(T['r_near_nm']*10, 2); M['TiltRfar'] = f(T['r_far_nm']*10, 2)
M['TiltNeutral'] = f(-T['dE_DFT_neutral']*1e3, 0)          # magnitude; neutral pair favours the NEAR configuration
M['NearSTgap'] = f(T['spin']['near_ST_gap_meV'], 0); M['FarSTgap'] = f(abs(T['spin']['far_ST_gap_meV']), 1)
M['TiltCharge'] = f(T['dE_DFT_charge_part']*1e3, 0)
M['TiltAgreeLo'] = f(100*(T['model_over_DFT_eps6_35'] - 1), 0); M['TiltAgreeHi'] = f(100*(T['model_over_DFT_eps5_6'] - 1), 0)
if 'lambda_eV' in N:
    M['LambdaMinus'] = f(N['lambda_eV']['minus']); M['LambdaPlus'] = f(N['lambda_eV']['plus'])
else:
    M['LambdaMinus'] = '\\pending'; M['LambdaPlus'] = '\\pending'

# retention / simulations
import sys; sys.path.insert(0, '../opt')
import fea_chip_sim as Sx
LDFT = round(N['lambda_eV']['minus'], 3)
for lam, tg in ((LDFT, 'Dft'), (1.0, 'One'), (1.5, 'OneHalf')):
    dE = Sx.best_level(lam)
    for TC, tt in ((27, 'RT'), (50, 'Fifty'), (85, 'EightyFive')):
        M[f'Tau{tg}{tt}'] = sci(Sx.tau_bit(TC + 273.15, lam, dE)) if Sx.tau_bit(TC + 273.15, lam, dE) < 0.01 else f(Sx.tau_bit(TC + 273.15, lam, dE), 3)
HRS = dict(Sx.DAY)   # hours per activity in the simulated day
def run(lam, mach, mode):
    return AD[f'{lam}/{mach}'][mode]
for lam, tg in ((str(LDFT), 'Dft'), ('1.0', 'One'), ('1.5', 'OneHalf')):
    for mach in ('phone', 'laptop', 'desktop'):
        mt = mach.capitalize()
        for mode in ('fixed', 'adaptive'):
            md = 'Fix' if mode == 'fixed' else 'Ad'
            M[f'Amb{tg}{mt}{md}'] = f(AD[f'{lam}/{mach}'][f'{mode}_max_ambient_gaming_C'], 0)
            r = run(lam, mach, mode)
            if r['runaway']:
                M[f'Peak{tg}{mt}{md}'] = '\\text{runaway}'; M[f'Day{tg}{mt}{md}'] = '\\text{runaway}'
            else:
                M[f'Peak{tg}{mt}{md}'] = f(max(x['P'] for x in r['rows']), 2)
                M[f'Day{tg}{mt}{md}'] = f(r['Wh'], 1)
                tot = sum(x['P']*HRS[x['work']] for x in r['rows']); ref = sum(x['P_ref']*HRS[x['work']] for x in r['rows'])
                M[f'RefShare{tg}{mt}{md}'] = f(100*ref/tot, 1)
                M[f'Arr{tg}{mt}{md}'] = f(1e3*max(x['P_ref'] + x['P_acc'] + x['P_ops'] for x in r['rows']), 0)
# refresh-energy saving (adaptive vs fixed), lambda = 1.0, laptop
def ref_E(lam, mach, mode):
    r = run(lam, mach, mode); return sum(x['P_ref']*HRS[x['work']] for x in r['rows'])
for mach in ('laptop', 'desktop'):
    M[f'RefSave{mach.capitalize()}'] = f(100*(1 - ref_E('1.0', mach, 'adaptive')/ref_E('1.0', mach, 'fixed')), 0)
M['RefSaveDftDesktop'] = f(100*(1 - ref_E(str(LDFT), 'desktop', 'adaptive')/ref_E(str(LDFT), 'desktop', 'fixed')), 0)
_d = lambda mode: run(str(LDFT), 'desktop', mode)['Wh']
M['DaySaveDftDesktop'] = f(100*(1 - _d('adaptive')/_d('fixed')), 0)
M['LambdaPlusTwoSides'] = f(N['cc']['lam']['DB+ at neutral geometry'], 2) + '--' + f(N['cc']['lam']['neutral, at DB+ geometry'], 2)
M['LambdaMinusA'] = f(N['cc']['lam']['DB- at neutral geometry'], 3); M['LambdaMinusB'] = f(N['cc']['lam']['neutral, at DB- geometry'], 3)
M['TauRatio'] = f(Sx.tau_bit(300, 1.5, Sx.best_level(1.5))/Sx.tau_bit(300, LDFT, Sx.best_level(LDFT)), 0)
M['LambdaShort'] = f(1.5/LDFT, 1)
M['TauDftRTms'] = f(1e3*Sx.tau_bit(300, LDFT, Sx.best_level(LDFT)), 1)
M['EcapDft'] = f(1e3*Sx.ecap(LDFT, Sx.best_level(LDFT)), 0)
M['RefSavePhone'] = f(100*(1 - ref_E('1.5', 'phone', 'adaptive')/ref_E('1.5', 'phone', 'fixed')), 0)
fails = [run(l, m, md)['fails_day'] for l in (str(LDFT), '1.0', '1.5') for m in ('phone', 'laptop', 'desktop') for md in ('fixed', 'adaptive')
         if not run(l, m, md)['runaway']]
M['FailsDay'] = sci(max(fails))
# FEA vs a 100 W budget
M['OpsAtHundredW'] = sci(100/(N['E_gate_aJ']*1e-18))

R5 = json.load(open('numbers_r5.json'))
M['GammaEl'] = f(R5['Gamma_el_meV'], 1); M['Tw'] = f(R5['tw_meV'], 1); M['Amax'] = f(R5['A_max'], 2)
M['AWrite'] = f(R5['write_tau0.1ps']['A'], 3); M['AttemptsW'] = str(R5['attempts_bit_9999'])
M['LogicRelax'] = f(R5['logic_relax_time_s']*1e6, 2); M['DsMin'] = f(R5['Ds_min_W100ps_eV'], 2)
M['AContact'] = f(R5['A_thermal_bandlimited_critical'], 2); M['ASource'] = f(R5['A_injector_sigma1meV'], 2)
R6 = json.load(open('numbers_r6.json'))
M['ControlShift'] = f(R6['control_shift_1.536_eV'], 2); M['PredLeak'] = sci(R6['predicate_leak_1.536'])
vg = R6['wire_len_um']/R6['transit_ps']                       # group velocity of the R6 wire model (um/ps)
wl = 64*N['pitch_lam5.0'][0]*1e-3                             # 64 bits along the wire at the R7 pitch (um)
M['WireLen'] = f(wl, 2); M['Transit'] = f(wl/vg, 0)

with open('numbers_r7.tex', 'w') as fh:
    fh.write('% generated by make_macros_r7.py -- do not edit\n')
    fh.write('\\providecommand{\\pending}{\\textbf{[pending]}}\n')
    for k, v in M.items():
        assert k.isalpha(), k
        fh.write(f'\\newcommand{{\\{k}}}{{{v}}}\n')
print(len(M), 'macros'); print({k: M[k] for k in list(M)[:12]})
