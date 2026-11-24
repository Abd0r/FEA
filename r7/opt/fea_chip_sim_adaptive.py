"""Adaptive FEA chip vs fixed (non-adaptive) chip, under real-world conditions.

Builds on fea_chip_sim.py (same device physics, error code, chip, energies).

Real-world effects added:
  * Hot spots: the die is split into 16 Zone groups whose temperatures spread around the chip mean.
    Spread = 2 K + 60% of the chip's rise above ambient (edge I/O transistors and busy Zones run hotter).
  * Memory occupancy: only part of the memory holds live data (OS + apps), varying with the workload.
  * Ambient that changes through the day (phone in a pocket, outdoors in the afternoon; laptop on a lap).
  * Access locality: program reads/writes land on the live data.

FIXED chip (what a simple controller must do to be safe): refresh every word, at the interval required
by the HOTTEST Zone, majority-of-5 sensing everywhere.

ADAPTIVE chip (all inside the FEA fabric, done by each Zone's FZC):
  A1 per-Zone temperature: each Zone group refreshes at its own temperature's interval
  A2 skip empty Zones: unallocated words are not refreshed (they hold no data)
  A3 skip fresh words: words read or written by programs since the last pass are already restored
  A4 adaptive sensing: majority-of-3 for program reads when the error budget allows, else majority-of-5
     (refresh passes always choose whichever sensing gives the lower total energy within the target)

Both chips must meet the same target: <= 1 uncorrectable word per year chip-wide.
Usage: python3 fea_chip_sim_adaptive.py [lambda_eV ...]
"""
import math, sys, json
from functools import lru_cache
import fea_chip_sim as B

G = 16                                                   # Zone groups (temperature bins)
# Sensing options (silent read error, energy per bit read+restore). R7 electron read = one logic-mode pass;
# a majority of 3 passes lowers the silent error further at 3x the read injections (A4 adaptive sensing).
_e1 = B.READ_ERR
SENSE = {1: (_e1, B.E_REFRESH_NOM),
         3: (3*_e1**2 - 2*_e1**3, B.E_REFRESH_NOM + 2*B.INJ_J)}

# live-data fraction per workload (share of memory holding data), typical OS + app footprints
OCC = {'sleep': 0.30, 'office/web': 0.50, 'video call': 0.55, 'compile': 0.65, 'gaming': 0.75, 'AI inference': 0.90}
# ambient through the day per machine (C): realistic use
AMB = {'phone':   {'sleep': 28, 'office/web': 33, 'video call': 33, 'compile': 33, 'gaming': 35, 'AI inference': 40},
       'laptop':  {'sleep': 25, 'office/web': 33, 'video call': 33, 'compile': 35, 'gaming': 35, 'AI inference': 35},
       'desktop': {'sleep': 22, 'office/web': 25, 'video call': 25, 'compile': 27, 'gaming': 27, 'AI inference': 27}}

@lru_cache(maxsize=None)
def pick(T, lam, dE, reads, sense, frac_words):
    """Longest refresh interval meeting this group's share of the yearly target, for a given sensing."""
    err, _ = SENSE[sense]
    tau = B.tau_bit(T, lam, dE); words = B.WORDS*frac_words
    budget = B.TARGET_FAILS_PER_YEAR/B.YEAR/G
    def fails(Tr):
        p = -math.expm1(-Tr/tau)
        return B.word_fail(p + B.DISTURB, err)*words/Tr + B.word_fail(p/2 + B.DISTURB, err)*words*reads
    # U-shaped in log(Tr): find the minimum by golden search, then bisect the long-interval branch
    lo, hi = -10.0, math.log10(10*tau)
    g = (math.sqrt(5) - 1)/2
    a_, b_ = lo, hi
    for _ in range(40):
        c, d = b_ - g*(b_ - a_), a_ + g*(b_ - a_)
        if fails(10**c) < fails(10**d): b_ = d
        else: a_ = c
    xmin = 0.5*(a_ + b_)
    ok = fails(10**xmin) < budget
    if not ok:
        Tr = 10**xmin
    else:
        l, h = xmin, hi
        if fails(10**h) < budget: l = h
        for _ in range(40):
            m_ = 0.5*(l + h)
            l, h = (m_, h) if fails(10**m_) < budget else (l, m_)
        Tr = 10**l
    return Tr, fails(Tr), bool(ok)

def zone_temps(T, Tamb):
    spread = 2.0 + 0.6*max(0.0, T - Tamb)
    return [T + spread*(g/(G - 1) - 0.5) for g in range(G)]

def chip_power(T, Tamb, lam, dE, wname, e_scale, adaptive):
    w = B.WORK[wname]; occ = OCC[wname]
    reads_live = w['traffic_GBs']*8e9/64/(B.WORDS*occ)                # accesses per LIVE word per s
    temps = zone_temps(T, Tamb)
    P_ref = 0.0; fails = 0.0; ok = True; info = []; senses = []
    for Tg in temps:
        Tq = round(Tg*2)/2
        if not adaptive:
            Tr, f, k = pick(round(max(temps)*2)/2, lam, dE, float(f'{reads_live*occ:.2g}'), 1, 1.0/G)
            words = B.WORDS/G; fresh = 0.0; e_bit = SENSE[1][1]
        else:
            best = None
            for s in (1, 3):
                Tr_s, f_s, k_s = pick(Tq, lam, dE, float(f'{reads_live:.2g}'), s, occ/G)
                fr = -math.expm1(-reads_live*Tr_s)                     # fraction refreshed by program traffic
                cost = (B.WORDS*occ/G)*92*(1 - fr)*SENSE[s][1]/Tr_s
                if k_s and (best is None or cost < best[0]): best = (cost, Tr_s, f_s, k_s, s, fr)
            if best is None:
                Tr, f, k = pick(Tq, lam, dE, float(f'{reads_live:.2g}'), 3, occ/G); s, fresh = 3, 0.0
            else:
                _, Tr, f, k, s, fresh = best
            words = B.WORDS*occ/G; e_bit = SENSE[s][1]; senses.append(s)
        P_ref += words*92*(1 - fresh)*e_bit*e_scale/Tr
        fails += f
        ok &= k; info.append(Tr)
    sense_reads = 1 if (not adaptive or all(x == 1 for x in senses)) else 3
    P_acc = w['traffic_GBs']*8e9*SENSE[sense_reads][1]
    P_ops = w['ops']*B.E_GATE
    return P_ref, P_acc, P_ops, fails, ok, min(info), max(info)

def run_day(machine, lam, adaptive, e_scale=1.0, dt=30.0):
    m = B.MACHINES[machine]; dE = B.best_level(lam)
    T = AMB[machine]['sleep'] + 273.15; rows = []; E = 0.0; F = 0.0; runaway = None; budget_ok = True
    for wname, hours in B.DAY:
        Tamb = AMB[machine][wname] + 273.15; w = B.WORK[wname]
        P_io = w['io_GBs']*8e9*B.IO_J_PER_BIT + m['P_edge']
        for _ in range(int(hours*3600/dt)):
            P_ref, P_acc, P_ops, f, ok, trmin, trmax = chip_power(T, Tamb, lam, dE, wname, e_scale, adaptive)
            P = P_ref + P_acc + P_ops + P_io
            if T > 273.15 + 150 or P > 50*m['power_cap']: runaway = wname; break
            T += (Tamb + m['Rth']*P - T)*(1 - math.exp(-dt/m['tau_th']))
            E += P*dt; F += f*dt; budget_ok &= ok
        if runaway: break
        rows.append(dict(work=wname, T_C=T - 273.15, hot_C=max(zone_temps(T, Tamb)) - 273.15, P=P, P_ref=P_ref,
                         P_acc=P_acc, P_ops=P_ops, P_io=P_io, refresh_fast=trmin, refresh_slow=trmax))
    return dict(rows=rows, Wh=E/3600, fails_day=F, runaway=runaway, target_met=budget_ok)

def max_ambient(machine, lam, adaptive, wname='gaming'):
    """Highest steady ambient (C) at which the chip still meets the target without runaway."""
    m = B.MACHINES[machine]; dE = B.best_level(lam); lo, hi = 0.0, 120.0
    for _ in range(25):
        mid = 0.5*(lo + hi); Tamb = mid + 273.15; T = Tamb; good = True
        w = B.WORK[wname]; P_io = w['io_GBs']*8e9*B.IO_J_PER_BIT + m['P_edge']
        for _ in range(400):
            P_ref, P_acc, P_ops, f, ok, *_ = chip_power(T, Tamb, lam, dE, wname, 1.0, adaptive)
            Tn = Tamb + m['Rth']*(P_ref + P_acc + P_ops + P_io)
            if Tn > 273.15 + 150: good = False; break
            if abs(Tn - T) < 1e-3: break
            T = 0.5*(T + Tn)
        good &= ok and abs(Tn - T) < 0.1
        lo, hi = (mid, hi) if good else (lo, mid)
    return lo

if __name__ == '__main__':
    lams = [float(x) for x in sys.argv[1:]] or [1.0, 1.5]
    out = {}
    for lam in lams:
        print(f'\n================ lambda = {lam} eV ================')
        for mach in ('phone', 'laptop', 'desktop'):
            res = {}
            for ad in (False, True):
                r = run_day(mach, lam, ad); res['adaptive' if ad else 'fixed'] = r
                lab = 'ADAPTIVE' if ad else 'fixed   '
                if r['runaway']:
                    print(f'{mach:8s} {lab} THERMAL RUNAWAY during {r["runaway"]}'); continue
                pk = max(r['rows'], key=lambda x: x['P'])
                print(f'{mach:8s} {lab} day {r["Wh"]:7.2f} Wh | peak {pk["P"]:6.3f} W, hottest Zone {pk["hot_C"]:5.1f} C '
                      f'({pk["work"]}) | refresh share of day energy '
                      f'{sum(x["P_ref"] for x in r["rows"])/max(1e-12, sum(x["P"] for x in r["rows"]))*100:4.1f}% | '
                      f'uncorrectable/day {r["fails_day"]:.2g} | target {"met" if r["target_met"] else "MISSED"}')
                for x in r['rows']:
                    print(f'         {x["work"]:13s} mean {x["T_C"]:5.1f} C hot {x["hot_C"]:5.1f} C  P {x["P"]:7.3f} W '
                          f'(refresh {x["P_ref"]:.3g}, access {x["P_acc"]:.3g}, compute {x["P_ops"]:.3g}, io {x["P_io"]:.3g})  '
                          f'refresh every {x["refresh_fast"]:.2g}-{x["refresh_slow"]:.2g} s')
            for ad in (False, True):
                res[('adaptive' if ad else 'fixed') + '_max_ambient_gaming_C'] = a = max_ambient(mach, lam, ad)
                print(f'{mach:8s} {"ADAPTIVE" if ad else "fixed   "} highest ambient while gaming: {a:5.1f} C')
            out[f'{lam}/{mach}'] = res
    json.dump(out, open('chip_sim_adaptive_results.json', 'w'), indent=1, default=float)
