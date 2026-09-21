#!/usr/bin/env python3
import os
from pathlib import Path

# Anchor every relative path below to the repository root, so this script
# can be run from any working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT_ROOT)

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# ====================================================
# STEP 1: Read CSV
# ====================================================
csv_file = "data/processed/affinity_corrarray.csv"
df_raw = pd.read_csv(csv_file, sep=",", engine="python", index_col=0)

# Drop rows/columns with null labels
df_raw = df_raw.loc[df_raw.index.notnull(), df_raw.columns.notnull()]

print("df_raw head:\n", df_raw.head())

methods = df_raw.columns  # method names
df_corr = pd.DataFrame(index=df_raw.index, columns=methods, dtype=float)
df_n    = pd.DataFrame(index=df_raw.index, columns=methods, dtype=float)
df_sig  = pd.DataFrame(index=df_raw.index, columns=methods, dtype=float)
df_col_sample = pd.DataFrame(index=df_raw.index, columns=methods, dtype=float)

# ====================================================
# STEP 2: Parse each cell of the correlation matrix.
#
# calibrate.py has written two formats over the life of this project:
#   5 fields: r  p  q  n  n_unique   (current; q is the BH-FDR value)
#   4 fields: r  p  n  n_unique      (before q was added)
# Accept both, so the shipped matrix and a freshly generated one both load.
# ====================================================
for r in df_raw.index:
    for c in df_raw.columns:
        raw_val = str(df_raw.loc[r, c])
        split_vals = raw_val.split()
        if len(split_vals) == 5:
            corr_str, sig_str, _q_str, row_samp_str, col_samp_str = split_vals
        elif len(split_vals) == 4:
            corr_str, sig_str, row_samp_str, col_samp_str = split_vals
        else:
            corr_str = sig_str = row_samp_str = col_samp_str = None

        if corr_str is not None:
            try:
                df_corr.loc[r, c] = float(corr_str)
                df_sig.loc[r, c]  = float(sig_str)
                df_n.loc[r, c]    = float(row_samp_str)
                df_col_sample.loc[r, c] = float(col_samp_str)
            except ValueError:
                df_corr.loc[r, c] = np.nan
                df_sig.loc[r, c]  = np.nan
                df_n.loc[r, c]    = np.nan
                df_col_sample.loc[r, c] = np.nan
        else:
            df_corr.loc[r, c] = np.nan
            df_sig.loc[r, c]  = np.nan
            df_n.loc[r, c]    = np.nan
            df_col_sample.loc[r, c] = np.nan

# Optionally remove row/column combos that have n <= 5
df_corr[df_n <= 5] = np.nan

# ====================================================
# STEP 3: Flatten into df_pairs
# ====================================================
methods_list = df_corr.index  # or df_corr.columns
pairs = []
for row_m in methods_list:
    for col_m in methods_list:
        if row_m == col_m:
            continue
        r_val = df_corr.loc[row_m, col_m]
        n_val = df_n.loc[row_m, col_m]
        p_val = df_sig.loc[row_m, col_m]
        pairs.append({
            "row_method": row_m,
            "col_method": col_m,
            "r": r_val,
            "n": n_val,
            "p": p_val
        })
df_pairs = pd.DataFrame(pairs)

# Drop any row with missing r, n, p
df_pairs.dropna(subset=["r","n","p"], inplace=True)

# ====================================================
# STEP 4: Clip extremely small p-values, build color array
# ====================================================
min_p = 1e-15
p = df_pairs["p"].values
p_clipped = np.maximum(p, min_p)
craw = -np.log10(p_clipped)  # raw -log10(p)

val_min = craw.min()
val_max = craw.max()

if np.isclose(val_min, val_max):
    # no variation
    cvals = np.full_like(craw, 0.5)
    print("Warning: 'craw' had no variation, using constant color=0.5")
else:
    # normalize to [0,1]
    cvals = (craw - val_min)/(val_max - val_min)

# Print debug
print("\n===== DF_PAIRS INFO =====")
print("df_pairs shape:", df_pairs.shape)
print(df_pairs.head(15))
print(f"\nMin/Max of r: {df_pairs['r'].min()}, {df_pairs['r'].max()}")
print(f"Min/Max of n: {df_pairs['n'].min()}, {df_pairs['n'].max()}")
print(f"Min/Max of p: {df_pairs['p'].min()}, {df_pairs['p'].max()}")
print(f"Range of cvals: [{cvals.min()}, {cvals.max()}]")

# ====================================================
# STEP 5: Create scatter plot
# ====================================================
plt.style.use("default")
fig, ax = plt.subplots(figsize=(6,4))  # bigger figure

sc = ax.scatter(
    x=df_pairs["r"],
    y=df_pairs["n"],
    c=cvals,         # normalized [0,1]
    cmap="viridis",
    alpha=1.0,       # no partial transparency
    s=120,           # bigger marker size
    edgecolors="none"
)

# Force axis limits to ensure data is visible
ax.set_xlim(-0.1, 1.1)   # r in range [0.0..1.0], give margin
ax.set_ylim(0, 60)       # n in range [9..46], give margin

# Build colorbar with even spacing in cvals
cbar = plt.colorbar(sc, ax=ax)
n_ticks = 5
tick_positions = np.linspace(0,1, n_ticks) # cvals from 0..1
cbar.set_ticks(tick_positions)

# Convert each cval back to p = 10^(- [unscaled cval]) using a linear map
# unscaled c = val_min + (val_max - val_min)* cvals
labels = []
for cval in tick_positions:
    unscaled = val_min + cval*(val_max - val_min)
    p_val = 10**(-unscaled)
    labels.append(f"{p_val:.2e}")
cbar.set_ticklabels(labels)
cbar.set_label("p-value   (color scale = -log10(p))")

ax.set_xlabel("Correlation Coefficient (r)")
ax.set_ylabel("Sample Size (n)")
ax.set_title("Correlation vs. Sample Size\n(color = p-value)")

plt.show()

# Save to PDF
out_pdf = "results/summary/correlation_scatter.pdf"
fig.savefig(out_pdf, bbox_inches="tight")
print(f"\nPlot saved to: {out_pdf}")

# ====================================================
# STEP 8: Bubble Plot on a Grid (Method x Method)
# ====================================================
# Build bubble size and x,y coords
df_pairs["marker_size"] = np.sqrt(df_pairs["n"]) * 50  # you can tweak scale
method_to_num = {m: i for i, m in enumerate(methods_list)}
df_pairs["x"] = df_pairs["col_method"].map(method_to_num)
df_pairs["y"] = df_pairs["row_method"].map(method_to_num)
df_pairs["r_clipped"] = df_pairs["r"].clip(-1, 1)

# Filter for significance
df_sig = df_pairs[df_pairs["p"] < 0.05]

fig2, ax2 = plt.subplots(figsize=(8,8), facecolor="white")
ax2.set_facecolor("white")

scatter = ax2.scatter(
    x=df_pairs["x"],
    y=df_pairs["y"],
    s=df_pairs["marker_size"],
    c=df_pairs["r_clipped"],
    cmap="coolwarm",
    alpha=1.0
)

cbar = plt.colorbar(scatter, ax=ax2)
cbar.set_label("Correlation Coefficient (r)")

# Outline significant pairs
ax2.scatter(
    x=df_sig["x"],
    y=df_sig["y"],
    s=df_sig["marker_size"],
    facecolors="none",
    edgecolors="black",
    linewidths=1.2
)

ax2.set_xticks(range(len(methods_list)))
ax2.set_xticklabels(methods_list, rotation=90)
ax2.set_yticks(range(len(methods_list)))
ax2.set_yticklabels(methods_list)
ax2.set_xlim(-0.5, len(methods_list) - 0.5)
ax2.set_ylim(len(methods_list) - 0.5, -0.5)  # invert
ax2.set_aspect("equal", "box")
ax2.set_title(
    "Correlation Bubble Plot\n"
    "(Size = sample size, Color = correlation, Black edge = p<0.05)",
    pad=20
)
plt.tight_layout()
fig2.savefig("results/summary/bubble.pdf")
# PNG alongside the PDF: GitHub cannot render PDFs in README.md.
fig2.savefig("results/summary/bubble.png", dpi=150, bbox_inches="tight")

# ====================================================
# STEP 9: Heatmap with Annotations
# ====================================================
fig3, ax3 = plt.subplots(figsize=(6,6), facecolor="white")
ax3.set_facecolor("white")

cax = ax3.matshow(df_corr, cmap="coolwarm", vmin=-1, vmax=1)
plt.colorbar(cax)

ax3.set_xticks(range(len(methods_list)))
ax3.set_yticks(range(len(methods_list)))
ax3.set_xticklabels(methods_list, rotation=90)
ax3.set_yticklabels(methods_list)

for i, row_m in enumerate(methods_list):
    for j, col_m in enumerate(methods_list):
        if row_m == col_m:
            continue
        r_val = df_corr.loc[row_m, col_m]
        n_val = df_n.loc[row_m, col_m]
        if pd.notna(r_val) and pd.notna(n_val):
            text = f"{r_val:.2f}\n(n={int(n_val)})"
            ax3.text(j, i, text, ha="center", va="center", fontsize=6, color="white")

ax3.set_title("Correlation Heatmap with (r, n)")
plt.tight_layout()
fig3.savefig("results/summary/heatmap.pdf")



import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Suppose df_corr is your NxN correlation matrix (numeric, no NaNs).
# For example, from df.corr() or your final correlation DataFrame.
print(df_corr)
import numpy as np
import pandas as pd

# 1) Drop any row or column that is entirely NaN in df_corr:
df_corr_dropped = df_corr.dropna(axis=0, how="all").dropna(axis=1, how="all")
df_corr_dropped = df_corr_dropped.fillna(0)

# 2) Optionally, also remove matching rows/columns in df_n, so shapes remain consistent
df_n_dropped = df_n.loc[df_corr_dropped.index, df_corr_dropped.columns]
print(df_corr_dropped)
if df_corr_dropped.isnull().values.any():
    print("Still have NaNs. Must drop or fill them before clustermap.")
else:
    print("No NaNs left. Ready for clustermap.")
sns.set(style="white")  # optional styling
g = sns.clustermap(
    df_corr_dropped,
    method="average",        # how to cluster (average, complete, etc.)
    metric="euclidean",      # distance metric (euclidean, 1 - correlation, etc.)
    cmap="vlag",             # colormap
    center=0,               # so 0 correlation is white, positive is red, negative is blue
    linewidths=.75,
    figsize=(8, 8)
)
g.savefig("results/summary/clusteredheatmap.pdf")
