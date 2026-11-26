"""Write reference_r7.json from the current outputs (run once, by the author, after a full build)."""
import json, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
N = json.load(open('paper/numbers_r7.json')); AD = json.load(open('opt/chip_sim_adaptive_results.json'))
ref = {
 'N:dft:DB0:Ea': (N['dft']['DB0']['Ea'], 1e-6), 'N:dft:DB-:Ea': (N['dft']['DB-']['Ea'], 1e-6),
 'N:dft:DB+:Ea': (N['dft']['DB+']['Ea'], 1e-6), 'N:dft:DB-:bridge_dE': (N['dft']['DB-']['bridge_dE'], 1e-6),
 'N:tilt:dE_DFT': (N['tilt']['dE_DFT'], 1e-6), 'N:tilt:dE_model_eps5_6': (N['tilt']['dE_model_eps5_6'], 1e-6),
 'N:raw_GB_lam5.0': (N['raw_GB_lam5.0'], 1e-6), 'N:chip:usable_GB': (N['chip']['usable_GB'], 1e-6),
 'N:E_refresh_bit_aJ': (N['E_refresh_bit_aJ'], 1e-6), 'N:read_disturb': (N['read_disturb'], 1e-6),
 'N:pattern_defects_per_year_chip:Hill_clean_meas': (N['pattern_defects_per_year_chip']['Hill_clean_meas'], 1e-6),
 'N:anneal_Tmax_C:0.0001:Hill_clean_meas': (N['anneal_Tmax_C']['0.0001']['Hill_clean_meas'], 1e-6),
}
ref['N:E_line_per_wire_op_aJ'] = (N['E_line_per_wire_op_aJ'], 1e-6)
ref['N:code:ml_check:s_max_no_fail'] = (N['code']['ml_check']['s_max_no_fail'], 0)
ref['N:code:ml_check:p_fail_s48'] = (N['code']['ml_check']['p_fail_s48'], 1e-9)
L = str(round(N['lambda_eV']['minus'], 3))
ref['N:lambda_eV:minus'] = (N['lambda_eV']['minus'], 1e-6)
ref['N:tilt:dE_DFT_charge_part'] = (N['tilt']['dE_DFT_charge_part'], 1e-6)
for lam in (L, '1.0', '1.5'):
    for m in ('phone', 'laptop', 'desktop'):
        for mode in ('fixed', 'adaptive'):
            ref[f'AD:{lam}/{m}:{mode}_max_ambient_gaming_C'] = (AD[f'{lam}/{m}'][f'{mode}_max_ambient_gaming_C'], 0.02)
json.dump(ref, open('reference_r7.json', 'w'), indent=1)
print(len(ref), 'reference values')
