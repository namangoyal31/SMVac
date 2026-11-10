"""E4b: zoomed view of the legacy objective S(R) near its minimum, plus fine-scale smoothness probe.
Outputs audit/results/action_landscape_zoom.png and prints roughness diagnostics."""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from audit import legacy_port as L

POINTS = [(125.1, 173.1), (115.0, 180.0), (134.75, 176.5), (120.0, 172.0)]
fig, axes = plt.subplots(2, len(POINTS), figsize=(4.3 * len(POINTS), 7))
for j, (Mh, Mt) in enumerate(POINTS):
    res = L.classify(Mh, Mt)
    run = res['run']
    x0 = res['logR_min']
    # coarse zoom: +-2.5 in ln R
    g = np.linspace(x0 - 2.5, x0 + 2.5, 301)
    S = np.array([L.action_at_R(run.table, run.mu_inst, math.exp(x)) for x in g])
    mu = np.log10(run.mu_inst) - g / np.log(10)
    ax = axes[0, j]
    ax.plot(mu, S, lw=1)
    ax.axhline(res['S_approx'], color='k', ls='--', lw=0.8, label=r'$8\pi^2/3|\lambda_{\rm eff,min}|$')
    ax.plot([math.log10(res['mu_R_min'])], [res['S_exact']], 'r*', ms=10, label='GSS result')
    ax.set_title(f"$M_h$={Mh}, $M_t$={Mt}", fontsize=10)
    ax.set_xlabel(r'$\log_{10}(\mu_R/{\rm GeV})$')
    ax.set_ylim(res['S_approx'] * 0.97, res['S_approx'] * 1.25)
    ax.grid(alpha=0.3)
    # fine probe: +-0.02 in ln R with 2001 points
    gf = np.linspace(x0 - 0.02, x0 + 0.02, 2001)
    Sf = np.array([L.action_at_R(run.table, run.mu_inst, math.exp(x)) for x in gf])
    d2 = np.diff(Sf, 2)
    ax2 = axes[1, j]
    ax2.plot(gf - x0, Sf - res['S_exact'], lw=0.7)
    ax2.set_xlabel(r'$\ln R - \ln R_{\rm GSS}$')
    ax2.set_ylabel(r'$S(R) - S_{\rm GSS}$')
    ax2.grid(alpha=0.3)
    # smooth reference: quadratic fit residual
    c = np.polyfit(gf - x0, Sf, 2)
    resid = Sf - np.polyval(c, gf - x0)
    print(f"(Mh={Mh}, Mt={Mt}): S_exact={res['S_exact']:.10g}  S_approx={res['S_approx']:.10g}  rel diff={res['S_exact']/res['S_approx']-1:+.4%}")
    print(f"   fine window +-0.02 in lnR: max|S - quadratic fit| = {np.abs(resid).max():.3e}  (rel {np.abs(resid).max()/res['S_exact']:.2e});"
          f"  #sign changes of 1st diff = {int((np.diff(np.sign(np.diff(Sf)))!=0).sum())}")
    print(f"   mu_inst={run.mu_inst:.4e} GeV, mu(lambda_eff min)={res['mu_min']:.4e} GeV, mu_R at GSS min={res['mu_R_min']:.4e} GeV")
axes[0, 0].legend(fontsize=7)
fig.tight_layout()
fig.savefig('audit/results/action_landscape_zoom.png', dpi=150)
