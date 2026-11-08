"""
Overlay / comparison plot between analytical and numerical stability maps.

Colouring rules:
  - Points where both methods AGREE on stability → plotted in that stability's colour
      Stable      = green
      Metastable  = yellow / goldenrod
      Unstable    = red
      Non-pert.   = black / grey
  - Points where the two methods DISAGREE → plotted in BLUE
      (these mark the boundary shift between analytical and numerical)

Both the full plot and a close-up around the experimental point are saved.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── colour map ────────────────────────────────────────────────────────────────
COLOUR = {
    1: 'green',   # stable
    2: 'yellow',  # metastable
    3: 'red',     # unstable
    4: 'black',   # non-perturbative
}
DIFF_COLOUR = 'blue'   # disagreement between methods

# ── helper: align grids ───────────────────────────────────────────────────────
from scipy.spatial import cKDTree

def align_data(df_ana, df_num):
    df_ana['_Mt']  = df_ana['Mt'].round(6)
    df_ana['_Mh']  = df_ana['Mh_calc'].round(6)
    df_num['_Mt']  = df_num['Mt'].round(6)
    df_num['_Mh']  = df_num['Mh_calc'].round(6)
    
    tree = cKDTree(df_num[['_Mt', '_Mh']].values)
    dists, indices = tree.query(df_ana[['_Mt', '_Mh']].values)
    
    # Accept matches within a small distance
    valid = dists < 1.0
    
    merged = pd.DataFrame({
        '_Mt': df_ana['_Mt'].values[valid],
        '_Mh': df_ana['_Mh'].values[valid],
        'S_ana': df_ana['Stability'].values[valid],
        'S_num': df_num['Stability'].values[indices[valid]]
    })
    return merged

# ── helper: build per-stability masks for agreeing points ─────────────────────
def plot_stability_diagram(ax, merged, s_agree=3, s_disagree=8, xlim=None, ylim=None):
    Mt  = merged['_Mt'].values
    Mh  = merged['_Mh'].values
    S_a = merged['S_ana'].values
    S_n = merged['S_num'].values

    agree    = S_a == S_n
    disagree = ~agree

    # Plot agreeing points – grouped by stability value for legend clarity
    for stab_val, label in [(1, 'Stable'), (2, 'Metastable'), (3, 'Unstable'), (4, 'Non-perturbative')]:
        mask = agree & (S_n == stab_val)
        if mask.sum() == 0:
            continue
        ax.scatter(Mh[mask], Mt[mask],
                   c=COLOUR[stab_val], alpha=0.5, s=s_agree,
                   label=f'{label} (both agree)', rasterized=True)

    # Plot disagreeing points on top in blue
    if disagree.sum() > 0:
        ax.scatter(Mh[disagree], Mt[disagree],
                   c=DIFF_COLOUR, alpha=0.9, s=s_disagree, zorder=10,
                   label='Disagreement (boundary shift)', rasterized=True)

    # Experimental point
    ax.scatter(125.10, 173.1, color='blue', marker='*', s=200,
               edgecolor='black', zorder=15,
               label='Experimental (125.1, 173.1)')

    ax.set_xlabel('Higgs Mass $M_h$ (GeV)', fontsize=12)
    ax.set_ylabel('Top Mass $M_t$ (GeV)', fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    ax.legend(loc='lower right', fontsize=9, markerscale=2)

# ── full plot ─────────────────────────────────────────────────────────────────
print("Loading broad analytical data for full plot …")
df_ana_full = pd.read_csv('data/analytical_data.csv')

print("Loading broad numerical data for full plot …")
df_num_full = pd.read_csv('data/numerical_data.csv')

print("Aligning broad grids …")
merged_full = align_data(df_ana_full, df_num_full)
print(f"Merged full points: {len(merged_full):,}")

print("Generating full overlay plot …")
fig, ax = plt.subplots(figsize=(11, 9))
plot_stability_diagram(ax, merged_full, s_agree=3, s_disagree=8)
ax.set_title('Vacuum Stability: Analytical vs Numerical Overlay\n'
             '(Blue = disagreement between methods)', fontsize=13)
fig.tight_layout()
fig.savefig('figures/overlay_plot.png', dpi=200, bbox_inches='tight')
plt.close(fig)
print("Saved figures/overlay_plot.png")

# ── close-up plot ─────────────────────────────────────────────────────────────
print("\nLoading dense analytical closeup data for closeup plot …")
df_ana_dense = pd.read_csv('data/analytical_dense_closeup.csv')

print("Loading dense numerical closeup data for closeup plot …")
df_num_dense = pd.read_csv('data/numerical_dense_closeup.csv')

print("Aligning dense closeup grids …")
merged_dense = align_data(df_ana_dense, df_num_dense)
print(f"Merged dense closeup points: {len(merged_dense):,}")

print("Generating close-up overlay plot …")
fig, ax = plt.subplots(figsize=(11, 9))
plot_stability_diagram(ax, merged_dense, s_agree=1, s_disagree=2, xlim=(110, 140), ylim=(160, 185))
ax.set_title('Vacuum Stability Overlay – Close-up\n'
             '(Blue = disagreement between methods)', fontsize=13)
fig.tight_layout()
fig.savefig('figures/overlay_closeup_plot.png', dpi=200, bbox_inches='tight')
plt.close(fig)
print("Saved figures/overlay_closeup_plot.png")

print("\nDone!")

