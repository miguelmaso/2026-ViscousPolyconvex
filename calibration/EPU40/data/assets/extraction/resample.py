"""Resample the vector-extracted EPU40 curves (raw.json, nominal strain / MPa)
onto equally spaced stretch grids and write CSVs in the VHB49 format (stress in Pa)."""
import json, numpy as np, csv, sys, os
raw = json.load(open(sys.argv[1])); out = sys.argv[2]; os.makedirs(out, exist_ok=True)
D_QS, D_CYC = 0.02, 0.04          # stretch increments (0.04 = native sampling of the paper)

def grid(a, b, d):                 # a -> b inclusive, step d (sign from b-a)
    n = int(round(abs(b - a) / d)); return a + np.sign(b - a) * d * np.arange(n + 1)

# quasi-static: Fig. 4b, SR = 0.00032 /s (drop the closing segment back to the origin)
qs = np.array(raw['4b']['0.00032']); qs = qs[:np.argmax(qs[:, 0]) + 1]
lam = 1 + qs[:, 0]; lam_q = grid(1.0, lam.max(), D_QS)
with open(f'{out}/quasi-static.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['id', 'stretch', 'stress'])
    for l, s in zip(lam_q, np.interp(lam_q, lam, qs[:, 1])):
        w.writerow([1, round(l, 4), round(max(s, 0.0) * 1e6, 3)])
    # id 2: equilibrium end points of the multi-step relaxation test, Fig. 10a (plus the origin)
    for e, s in [(0.0, 0.0)] + raw['10a']['multistep']:
        w.writerow([2, round(1 + e, 4), round(s * 1e6, 3)])

# cyclic: Figs. 5a,5b,6a,6b
rows, tid, summary = [], 10, []
for panel, lmax in [('5a', 2), ('5b', 3), ('6a', 4), ('6b', 5)]:
    for rate in ['0.384', '0.128', '0.032']:
        a = np.array(raw[panel][rate]); lam, P = 1 + a[:, 0], a[:, 1]
        turn = np.where(np.abs(lam - lam.max()) < 1e-3)[0]
        ld_l, ld_P = lam[:turn[0] + 1], P[:turn[0] + 1]          # loading up to (incl.) the peak
        ul_l, ul_P = lam[turn[-1]:][::-1], P[turn[-1]:][::-1]     # unloading from the post-reversal point
        g_ld = grid(1.0, float(lmax), D_CYC)
        g_ul = grid(float(lmax), ul_l.min(), D_CYC)[1:]
        g_ul = g_ul[g_ul >= ul_l.min() - 1e-9]
        s_ld = np.interp(g_ld, ld_l, ld_P); s_ul = np.interp(g_ul, ul_l, ul_P)
        for l, s in list(zip(g_ld, s_ld)) + list(zip(g_ul, s_ul)):
            rows.append([tid, round(l, 4), round(s * 1e6, 3), float(rate)])
        summary.append((tid, lmax, rate, len(g_ld) + len(g_ul), s_ld.max(), P[turn[-1]] if len(turn) > 1 else None, g_ul.min()))
        tid += 1
with open(f'{out}/cyclic-loading.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['id', 'stretch', 'stress', 'rate']); w.writerows(rows)
for s in summary: print(s)
