"""E5: impact of each identified physics defect, one at a time, at the central point of
Buttazzo et al. (arXiv:1307.3536): M_h = 125.15, M_t = 173.34, alpha_s = 0.1184.
Compares couplings at M_Pl = 1.22e19 GeV with their eq. (61) and reports lambda_eff,min,
the lambda_eff zero crossing and the FL action.  Output: audit/results/defect_impact.csv"""
import math, csv
import numpy as np
from scipy.optimize import brentq
from audit import legacy_port as L, independent as I

MH, MT = 125.15, 173.34
REF = dict(g1=0.6154, g2=0.5055, g3=0.4873, yt=0.3825, lam=-0.0143)  # eq. (61)
CONFIGS = [
    ("legacy", {}),
    ("4-loop QCD sign fixed", dict(flip_4loop_sign=True)),
    ("3-loop beta_lambda y_t^2 factor fixed", dict(fix_3loop_lambda=True)),
    ("3-loop beta_yt y_t^2 factor fixed", dict(fix_3loop_yt=True)),
    ("matching: Buttazzo eqs. 55-60", dict(matching="buttazzo2013")),
    ("Z-boson coupling in lambda_eff fixed", dict(fix_z_coupling=True)),
    ("all of the above", dict(flip_4loop_sign=True, fix_3loop_lambda=True, fix_3loop_yt=True,
                              matching="buttazzo2013", fix_z_coupling=True)),
]
rows = []
for name, kw in CONFIGS:
    opt = L.Options(**kw)
    rg = I.RGSolution(MH, MT, opt)
    g1, g2, g3, yt, yb, ytau, lam = rg.c(rg.tP)
    pot = I.FieldPotential(rg)
    lmin, lnphi = I.lambda_min(pot)
    # first zero crossing of lambda_eff above max(M_t, v)
    ts = np.linspace(pot.t_lo, pot.t_hi, 4000); le = rg.lambda_eff(ts)
    k = np.nonzero(le < 0)[0]
    t_inst = brentq(lambda t: float(rg.lambda_eff(t)), ts[k[0] - 1], ts[k[0]], xtol=1e-12) if k.size and k[0] > 0 else float('nan')
    # MS-bar lambda minimum (for the scale-choice systematic)
    lam_ms = rg.c(ts)[6]; i = int(np.argmin(lam_ms))
    row = dict(config=name, g1_MPl=g1, g2_MPl=g2, g3_MPl=g3, yt_MPl=yt, lambda_MPl=lam,
               lambda_eff_min=lmin, log10_mu_at_lambda_eff_min=lnphi / math.log(10),
               log10_mu_lambda_eff_zero=t_inst / 2 / math.log(10),
               S_FL_lambda_eff=8 * math.pi**2 / 3 / abs(lmin),
               lambda_MSbar_min=float(lam_ms[i]), S_FL_lambda_MSbar=8 * math.pi**2 / 3 / abs(lam_ms[i]))
    rows.append(row)
print(f"{'config':40s} {'g1':>7s} {'g2':>7s} {'g3':>7s} {'yt':>7s} {'lam(MPl)':>9s} {'lameff_min':>10s} {'log10Linst':>10s} {'S_FL':>8s} {'lamMS_min':>9s}")
print(f"{'Buttazzo et al. 2013 eq.(61)':40s} {REF['g1']:7.4f} {REF['g2']:7.4f} {REF['g3']:7.4f} {REF['yt']:7.4f} {REF['lam']:9.4f}")
for r in rows:
    print(f"{r['config']:40s} {r['g1_MPl']:7.4f} {r['g2_MPl']:7.4f} {r['g3_MPl']:7.4f} {r['yt_MPl']:7.4f} {r['lambda_MPl']:9.5f} {r['lambda_eff_min']:10.6f} {r['log10_mu_lambda_eff_zero']:10.3f} {r['S_FL_lambda_eff']:8.1f} {r['lambda_MSbar_min']:9.5f}")
# integrator check: legacy fixed-step RK4 (dt=0.1) vs DOP853 for the legacy config
run = L.run_rge(MH, MT)
rg = I.RGSolution(MH, MT)
ydop = rg.c(rg.tP); yrk = run.table.y[-1]
print("\nRK4(dt=0.1) vs DOP853(rtol=1e-12) at M_Pl, max |diff| over couplings:", max(abs(a - b) for a, b in zip(ydop, yrk)))
for dt in (0.2, 0.1, 0.05):
    r2 = L.run_rge(MH, MT, L.Options(dt=dt))
    print(f"   dt={dt}: lambda(M_Pl) RK4 - DOP853 = {r2.table.y[-1][6] - ydop[6]:+.3e},  lambda_eff_min(grid) = {r2.min_lambda_eff:.10f}")
with open('audit/results/defect_impact.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
