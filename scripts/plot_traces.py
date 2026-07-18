import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Read the summary
df_sum = pd.read_csv("data/sweep_summary.csv")
print("Summary of trajectory statuses:")
print(df_sum['status'].value_counts())

plt.figure(figsize=(12, 8))

colors = ['blue', 'green', 'red']
labels = ['Near Barrier', 'Mid-region', 'Near Minimum']

for i in range(3):
    fname = f"results/trace_{i}.csv"
    if os.path.exists(fname):
        df = pd.read_csv(fname)
        # Plot phi(r)
        plt.plot(df['r'], df['phi'], label=f"Trace {i} ({labels[i]})", color=colors[i])

plt.xscale('log')
plt.yscale('symlog', linthresh=1e2)
plt.xlabel('Radius r')
plt.ylabel('Field \u03d5(r)')
plt.title('Shooting Trajectories in the Instability Interval')
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)

plt.savefig("figures/trace_plot.png", dpi=300)
print("Saved plot to figures/trace_plot.png")
