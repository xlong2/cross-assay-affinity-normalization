import os
from pathlib import Path

# Anchor every relative path below to the repository root, so this script
# can be run from any working directory.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT_ROOT)

import pandas as pd

# Using readlines()
import numpy as np
from scipy.stats import pearsonr

from scipy.stats import pearsonr, spearmanr
# keep this import; remove any 'scipy.stats.stats' import

def pearson_permutation_test(x, y, n_permutations=10000, alternative='two-sided', random_seed=None):
    """
    Permutation test for Pearson correlation. Returns (r_obs, p_value).
    Guards against constant vectors and uses RNG.permutation for speed.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    rng = np.random.default_rng(random_seed)

    # Not enough variation → no correlation possible
    if x.size < 3 or np.all(x == x[0]) or np.all(y == y[0]):
        return 0.0, 1.0

    r_obs, _ = pearsonr(x, y)
    r_obs_use = abs(r_obs) if alternative == 'two-sided' else r_obs

    count_extreme = 0
    for _ in range(n_permutations):
        y_perm = rng.permutation(y)
        r_perm, _ = pearsonr(x, y_perm)
        if alternative == 'two-sided':
            if abs(r_perm) >= r_obs_use: count_extreme += 1
        elif alternative == 'greater':
            if r_perm >= r_obs_use:       count_extreme += 1
        elif alternative == 'less':
            if r_perm <= r_obs_use:       count_extreme += 1
        else:
            raise ValueError("alternative must be 'two-sided', 'greater', or 'less'.")
    p_value = (count_extreme + 1) / (n_permutations + 1)
    return r_obs, p_value

#skempi_table=pd.read_table("data/raw/skempi_v2.csv", delimiter=",")
skempi_table=pd.read_table("data/raw/skempi_v2_temperature_parsed.csv", delimiter=",")

skempi_table["pdb_mutations"] = skempi_table['#Pdb']+ skempi_table['Mutation(s)_PDB']


skempi_table_non_binding_filtered = skempi_table.loc[skempi_table['Affinity_mut (M)'] !="n.b."]
group_by_method = skempi_table_non_binding_filtered.groupby("Method")

method_names = list(group_by_method.groups.keys())  # e.g. ["SPR","ITC",...]


import numpy as np

import matplotlib.pyplot as plt

# ADD RIGHT BELOW:
import os
for _d in ("data/processed", "results/summary", "results/calibration",
           "results/correlation_affinity", "results/correlation_ddg",
           "results/correlation_ddg_no_nonbinding"):
    os.makedirs(_d, exist_ok=True)

import itertools


fig = plt.figure()
index=0
ax = fig.add_subplot(1, 1, 1)


ax.set(xlabel="methods", ylabel="affinity")

vals= [group['Delta_G_mut '] for name, group in group_by_method]
xs= [group['Method'] for name, group in group_by_method]

names = [name for name, group in group_by_method]
plt.boxplot(vals, labels=names)
for i, arr in enumerate(vals, start=1):
    plt.scatter([i]*len(arr), arr, alpha=0.4)
#plt.xticks(np.asarray([x for x in range(20)]), [name for name, group in group_by_method])
plt.savefig(f"results/summary/mut_delta_G_no_nonbinding_valueplot_fig.pdf",
            bbox_inches='tight', dpi=100)
plt.close()






df1 = dict(tuple(skempi_table.groupby(skempi_table['pdb_mutations'])))
interested_tuple = []
# all possible methods
import numpy as np


for key, value in df1.items():
    if value.shape[0]>1:
        if len(np.unique(value['Method']))>1:
            print("has more than one method")
            interested_tuple.append((key, value))




all_methods = np.unique(skempi_table['Method'].to_numpy().astype("str"))
dummyarray = np.empty((len(interested_tuple),len(all_methods)))
dummyarray[:] = np.nan

affinity_difference_df = pd.DataFrame(dummyarray, [x for x,_ in interested_tuple], all_methods)
for key, value in interested_tuple:
    for id, each in value.iterrows():
        assert isinstance(key, str)
        assert isinstance(each['Method'], str)
        if each['Method']=="SPR":
            print("this")
        
        
        
        
        
        
        wt = each['Affinity_wt_parsed']
        mut = each['Affinity_mut_parsed']

        if (wt is not None) and (not np.isnan(wt)) and (wt != 0.0):
            ratio = mut / wt
            affinity_difference_df.loc[key, each['Method']] = ratio
        else:
            affinity_difference_df.loc[key, each['Method']] = np.nan
affinity_difference_df = affinity_difference_df.where(affinity_difference_df > 0, np.nan)
        
affinity_difference_df = np.log10(affinity_difference_df)
         
# --- delta-delta-G branch: build the ddG matrices -------------------
dummyarray_ddg = np.empty((len(interested_tuple),len(all_methods)))
dummyarray_ddg[:] = np.nan

delta_delta_G_difference_df = pd.DataFrame(dummyarray_ddg, [x for x,_ in interested_tuple], all_methods)
for key, value in interested_tuple:
    for id, each in value.iterrows():
        assert isinstance(key, str)
        assert isinstance(each['Method'], str)
        if each['Method']=="SPR":
            print("check")
        delta_delta_G_difference_df.loc[key][each['Method']]= each['Delta_Delta_G']




dummyarray_ddg_nononbinding = np.empty((len(interested_tuple),len(all_methods)))
dummyarray_ddg_nononbinding[:] = np.nan

delta_delta_G_difference_nononbinding_df = pd.DataFrame(dummyarray_ddg_nononbinding, [x for x,_ in interested_tuple], all_methods)
for key, value in interested_tuple:
    for id, each in value.iterrows():
        assert isinstance(key, str)
        assert isinstance(each['Method'], str)
        if each['Method']=="SPR":
            print("check")
        if each["Affinity_mut (M)"]!="n.b":

            delta_delta_G_difference_nononbinding_df.loc[key][each['Method']]= each['Delta_Delta_G']
        else:
            print(each)
# --- end ddG matrix construction ------------------------------------

# import itertools package
import itertools
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# Pearson r between assay columns of the log10 fold-change matrix.
# NOT an energy quantity - the ddG equivalent is written further below to
# data/processed/energy_difference_correlation.csv.
rho = affinity_difference_df.corr()
rho.to_csv("data/processed/fold_change_correlation.csv")
from itertools import permutations
permut = itertools.combinations(affinity_difference_df.columns, 2)
corrarray = pd.DataFrame(index=affinity_difference_df.columns,
                         columns=affinity_difference_df.columns,
                         dtype=object)



import matplotlib.pyplot as plt

permut = itertools.combinations(affinity_difference_df.columns, 2)


import itertools
import numpy as np
import matplotlib.pyplot as plt
def _safe(s: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-+" else "_" for ch in str(s))

# --- BH FDR helper (pure Python; no statsmodels needed) ---
def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    m = p.size
    order = np.argsort(p)
    p_sorted = p[order]
    q_sorted = p_sorted * m / (np.arange(1, m+1))
    # monotone (from largest to smallest), then clip to [0,1]
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty_like(q_sorted)
    q[order] = np.clip(q_sorted, 0.0, 1.0)
    return q

# ----------------------------------------------------------
# Build records, then apply BH across the actually-tested pairs
pair_records = []  # will hold tuples: (m1, m2, r, p, n, n_unique)

permut = itertools.combinations(affinity_difference_df.columns, 2)
stats_map = {}  # add this before the loop

for method1, method2 in permut:
    # vectorized overlap mask
    v1 = affinity_difference_df[method1].to_numpy()
    v2 = affinity_difference_df[method2].to_numpy()
    idx = affinity_difference_df.index.to_numpy()

    mask = ~(np.isnan(v1) | np.isnan(v2))
    val = np.column_stack((v1[mask], v2[mask]))  # shape (n, 2)
    mutations = idx[mask]

    val_len = val.shape[0]
    unique_mutation_len = np.unique(mutations).size

    if val_len > 1:
        r_value, p_value = pearson_permutation_test(
            val[:, 0], val[:, 1],
            n_permutations=10000,
            alternative='greater'
        )
        
        
        pair_records.append((
            method1, method2,
            float(r_value), float(p_value),
            int(val_len), int(unique_mutation_len)
        ))
        stats_map[(method1, method2)] = (float(r_value), float(p_value), int(val_len), int(unique_mutation_len))


        # (optional) quick scatter as before
        #fig = plt.figure()
        #ax = fig.add_subplot(1, 1, 1)
        #ax.scatter(val[:, 0], val[:, 1])
        #ax.set_title(f'{method1} vs {method2}')
        #ax.set(xlabel=method1, ylabel=method2)
        #plt.savefig(f"results/correlation_affinity/{_safe(method1)}_{_safe(method2)}_correlation_fig.pdf",
        #    bbox_inches='tight', dpi=100)
        #plt.close()
    else:
        corrarray.at[method1, method2] = "NA"
        corrarray.at[method2, method1] = "NA"


# --- Apply Benjamini–Hochberg across tested pairs only ---
if pair_records:
    pvals = [(rec[3] if np.isfinite(rec[3]) else 1.0) for rec in pair_records]

    qvals = bh_fdr(pvals)  # same order as pair_records

    for (method1, method2, r, p, n, nuniq), q in zip(pair_records, qvals):
        entry = f"{r:.7g}  {p:.7g}  {q:.7g}  {n}  {nuniq}"
        corrarray.at[method1, method2] = entry
        corrarray.at[method2, method1] = entry

# Save matrix (now r p q n n_unique per tested cell)
corrarray.to_csv("data/processed/affinity_corrarray.csv")


# Minimal reader to extract r/p/q/n/nuniq fields into matrices
import pandas as pd, numpy as np

csv_path = "data/processed/affinity_corrarray.csv"
df_raw = pd.read_csv(csv_path, index_col=0, dtype=str)
def _get_field(cell, k):  # k=0:r, 1:p, 2:q, 3:n_overlap, 4:n_unique_mut
    parts = str(cell).split()
    if len(parts) >= 5:
        try:
            return float(parts[k])
        except Exception:
            return np.nan
    return np.nan

df_r = df_raw.map(lambda s: _get_field(s, 0))
df_p = df_raw.map(lambda s: _get_field(s, 1))
df_q = df_raw.map(lambda s: _get_field(s, 2))
df_n = df_raw.map(lambda s: _get_field(s, 3))
df_r = df_r.loc[df_r.index.notna(), :]
df_r = df_r.loc[:, df_r.columns.notna()]
df_p = df_p.loc[df_p.index.notna(), df_p.columns.notna()]
df_q = df_q.loc[df_q.index.notna(), df_q.columns.notna()]
df_n = df_n.loc[df_n.index.notna(), df_n.columns.notna()]

# 2) Optionally drop empty/placeholder labels like '', 'nan', 'None', 'Unnamed:*'
def _bad_label(x):
    s = str(x).strip().lower()
    return (s == "") or (s == "nan") or (s == "none") or s.startswith("unnamed")

good_rows = [i for i in df_q.index if not _bad_label(i)]
good_cols = [c for c in df_q.columns if not _bad_label(c)]

df_r = df_r.loc[good_rows, good_cols]
df_p = df_p.loc[good_rows, good_cols]
df_q = df_q.loc[good_rows, good_cols]
df_n = df_n.loc[good_rows, good_cols]
# (optional) save for debugging/inspection
df_q.to_csv("data/processed/affinity_q_matrix.csv")


# ==== Build long-form pairs table from df_r/df_q/df_n ====
methods = list(df_q.index)
pairs = []
for i, m1 in enumerate(methods):
    for j, m2 in enumerate(methods):
        if j <= i:
            continue  # upper triangle only
        q = df_q.at[m1, m2]
        r = df_r.at[m1, m2]
        n = df_n.at[m1, m2]
        if np.isnan(q):
            continue
        pairs.append({"i": i, "j": j, "m1": m1, "m2": m2, "q": q, "r": r, "n": n})

pairs = pd.DataFrame(pairs)
if pairs.empty:
    print("No pairwise entries with q-values; skipping bubble plot.")
else:
    print(f"Pairs with q-values: {len(pairs)}")



for method1, method2 in itertools.combinations(methods, 2):
    # recompute the overlap & values for plotting
    v1 = affinity_difference_df[method1].to_numpy()
    v2 = affinity_difference_df[method2].to_numpy()
    mask = ~(np.isnan(v1) | np.isnan(v2))
    val = np.column_stack((v1[mask], v2[mask]))
    if val.shape[0] <= 1:
        continue

    # fetch stats
    r_plot, p_plot, n_plot, nuniq_plot = stats_map.get((method1, method2), (np.nan, np.nan, np.nan, np.nan))
    q_plot = df_q.at[method1, method2] if (method1 in df_q.index and method2 in df_q.columns) else np.nan

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(val[:, 0], val[:, 1])
    ax.set(xlabel=method1, ylabel=method2)
    ax.set_title(f'{method1} vs {method2}\nr={r_plot:.2f}, p={p_plot:.2g}, q={q_plot:.2g}, n={int(n_plot) if np.isfinite(n_plot) else "NA"}')
    plt.savefig(f"results/correlation_affinity/{_safe(method1)}_{_safe(method2)}_correlation_fig.pdf",
                bbox_inches='tight', dpi=100)
    plt.close()


# Example: make a -log10(q) matrix for coloring
neglog10_q = -np.log10(df_q.clip(lower=1e-15))
# ----------------------------------------------------------------

# ==== Bubble plot using FDR ====




def bubble_sizes(n_array, s_min=120, s_max=2400, clip=(0.10, 0.95), gamma=0.6):
    """
    Map overlap counts n to scatter areas (pt^2).
    - clip: lower/upper quantiles to clip extreme values
    - gamma<1 expands mid-range (more separation), gamma>1 compresses
    """
    n = np.asarray(n_array, float)
    n = np.nan_to_num(n, nan=0.0)

    # quantile clipping
    lo = np.nanquantile(n, clip[0]) if n.size else 0.0
    hi = np.nanquantile(n, clip[1]) if n.size else 1.0
    n_clip = np.clip(n, lo, hi)

    # normalize to 0..1
    if hi > lo:
        u = (n_clip - lo) / (hi - lo)
    else:
        u = np.zeros_like(n_clip)

    # power stretch and map to area range
    u = u**gamma
    sizes = s_min + u * (s_max - s_min)
    return sizes

# Visual encodings
cvals = -np.log10(np.clip(pairs["q"].to_numpy(), 1e-15, None))

sizes = bubble_sizes(pairs["n"].to_numpy(),
                     s_min=120, s_max=2400, clip=(0.10, 0.95), gamma=0.6)


edges = np.where(pairs["q"].to_numpy() < 0.05, "k", "none")

# sizes, colors, edges
sizes = bubble_sizes(pairs["n"].to_numpy(), s_min=120, s_max=2400, clip=(0.10, 0.95), gamma=0.6)
cvals = -np.log10(np.clip(pairs["q"].to_numpy(), 1e-15, None))
sig = pairs["q"].to_numpy() < 0.05
edgecols = np.where(sig, "k", "none")
edgeline = np.where(sig, 1.8, 0.0)

fig, ax = plt.subplots(figsize=(0.6*len(methods)+2, 0.6*len(methods)+2), facecolor="white")

sc = ax.scatter(pairs["j"].to_numpy(), pairs["i"].to_numpy(),
                s=sizes, c=cvals, cmap="viridis",
                edgecolors=edgecols, linewidths=edgeline, alpha=0.9)

# axes cosmetics
ax.set_xticks(range(len(methods))); ax.set_yticks(range(len(methods)))
ax.set_xticklabels(methods, rotation=45, ha="right"); ax.set_yticklabels(methods)
ax.set_xlim(-0.5, len(methods)-0.5); ax.set_ylim(-0.5, len(methods)-0.5); ax.invert_yaxis()
ax.set_title("Assay pair agreement (color = -log10(q), edge = q<0.05, size ∝ overlap)")

# colorbar in its own axis (no overlap with legend)
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="4%", pad=0.6)
cbar = plt.colorbar(sc, cax=cax); cbar.set_label("-log10(q)")

# size legend INSIDE the plot
n_ticks = np.array([np.nanmin(pairs["n"]), np.nanmedian(pairs["n"]), np.nanmax(pairs["n"])], float)
n_ticks = np.unique(np.round(n_ticks).astype(int))
s_ticks = bubble_sizes(n_ticks, s_min=120, s_max=2400, clip=(0.10, 0.95), gamma=0.6)
handles = [ax.scatter([], [], s=s, facecolors='none', edgecolors='k', linewidths=1.2) for s in s_ticks]
labels  = [f"n={t}" for t in n_ticks]
ax.legend(handles, labels, title="Overlap (n)", scatterpoints=1, loc="upper left", frameon=False)

# plain integer labels on larger bubbles only
label_thresh = np.quantile(sizes, 0.40) if sizes.size else np.inf
for x, y, n, s in zip(pairs["j"], pairs["i"], pairs["n"], sizes):
    if s >= label_thresh and np.isfinite(n):
        ax.text(float(x), float(y), f"{int(n)}", ha="center", va="center",
                fontsize=8, color="black", zorder=3)

plt.tight_layout()
plt.savefig("results/summary/correlation_bubble_FDR.pdf", dpi=150)
plt.close(fig)

# --- delta-delta-G branch: cross-assay ddG agreement -----------------


# Pearson r between assay columns of the ddG matrix - the actual energy quantity.
rho = delta_delta_G_difference_df.corr()
rho.to_csv("data/processed/energy_difference_correlation.csv")
from itertools import permutations
permut = itertools.combinations(delta_delta_G_difference_df.columns, 2)


corrarray_spearman  = np.empty((len(delta_delta_G_difference_df.columns), len(delta_delta_G_difference_df.columns)))

corrarray_spearman[:]= None
corrarray_spearman = pd.DataFrame(corrarray_spearman, [x for x in delta_delta_G_difference_df.columns],[x for x in  delta_delta_G_difference_df.columns] )

import matplotlib.pyplot as plt

permut = itertools.combinations(delta_delta_G_difference_df.columns, 2)


index=0
for method1, method2 in permut:
    print([method1, method2])
    mutations=[]
    for v1, v2,mut_name in zip(delta_delta_G_difference_df[method1].to_numpy(), delta_delta_G_difference_df[method2].to_numpy(), delta_delta_G_difference_df.index):
        #print([a, b])
        if not (np.isnan(v1) or np.isnan(v2)):

            #print([v1,v2,mut_name])
            mutations.append(mut_name)
    val = [(a,b, c) for a, b, c in zip(delta_delta_G_difference_df[method1].to_numpy(), delta_delta_G_difference_df[method2].to_numpy(), delta_delta_G_difference_df.index) if not (np.isnan(a) or np.isnan(b)) ]

    val_len = len(val)
    unique_mutation_len = len(np.unique(mutations))
    if val_len>1:
        val= np.asarray(val)[:,0:2].astype(float)
        corrarray[method1][method2] =str('  '.join([str(round(x,3)) for x in [pearsonr(np.asarray(val)[:,0], np.asarray(val)[:,1])[0], val_len, unique_mutation_len]]))
        corrarray[method2][method1] =str('  '.join([str(round(x,3)) for x in [pearsonr(np.asarray(val)[:,0], np.asarray(val)[:,1])[0], val_len, unique_mutation_len]]))
        corrarray_spearman[method2][method1] =str('  '.join([str(round(x,3)) for x in [spearmanr(np.asarray(val)[:,0], np.asarray(val)[:,1])[0], val_len, unique_mutation_len]]))

        corrarray_spearman[method1][method2] =str('  '.join([str(round(x,3)) for x in [spearmanr(np.asarray(val)[:,0], np.asarray(val)[:,1])[0], val_len, unique_mutation_len]]))

        if True:
            fig = plt.figure()

            ax=fig.add_subplot(1, 1, 1)
            total= np.asarray(val)[:, 0].tolist() + np.asarray(val)[:, 1].tolist()
            max_v = max(total)
            min_v = min(total)

            ax.set(xlabel=method1, ylabel=method2)

            ax.scatter(np.asarray(val)[:, 0], np.asarray(val)[:, 1])
            #ax.axis('equal')
            ax.set_xlim(left=min_v, right=max_v)
            ax.set_ylim(bottom=min_v, top=max_v)


            ax.set_title('Axis [0, 0]')
            plt.savefig(f"results/correlation_ddg/{method1}_{method2}_correlation_fig.pdf", bbox_inches='tight', dpi=100)
            plt.close()
        index = index + 1
    else:
        corrarray[method1][method2] ="NA"#str('/'.join([str(x) for x in ["NA", val_len,unique_mutation_len]]))

corrarray.to_csv("data/processed/delta_delta_G_pearson_corrarray.csv")

corrarray_spearman.to_csv("data/processed/delta_delta_G_spearman_corrarray.csv")







corrarray_spearman  = np.empty((len(delta_delta_G_difference_nononbinding_df.columns), len(delta_delta_G_difference_nononbinding_df.columns)))

corrarray_spearman[:]= None
corrarray_spearman = pd.DataFrame(corrarray_spearman, [x for x in delta_delta_G_difference_nononbinding_df.columns],[x for x in  delta_delta_G_difference_nononbinding_df.columns] )

import matplotlib.pyplot as plt

permut = itertools.combinations(delta_delta_G_difference_nononbinding_df.columns, 2)


index=0
for method1, method2 in permut:
    print([method1, method2])
    mutations=[]
    for v1, v2,mut_name in zip(delta_delta_G_difference_nononbinding_df[method1].to_numpy(), delta_delta_G_difference_nononbinding_df[method2].to_numpy(), delta_delta_G_difference_nononbinding_df.index):
        #print([a, b])
        if not (np.isnan(v1) or np.isnan(v2)):

            #print([v1,v2,mut_name])
            mutations.append(mut_name)
    val = [(a,b, c) for a, b, c in zip(delta_delta_G_difference_nononbinding_df[method1].to_numpy(), delta_delta_G_difference_nononbinding_df[method2].to_numpy(), delta_delta_G_difference_nononbinding_df.index) if not (np.isnan(a) or np.isnan(b)) ]

    val_len = len(val)
    unique_mutation_len = len(np.unique(mutations))
    if val_len>1:
        val= np.asarray(val)[:,0:2].astype(float)
        corrarray[method1][method2] =str('  '.join([str(round(x,3)) for x in [pearsonr(np.asarray(val)[:,0], np.asarray(val)[:,1])[0], val_len, unique_mutation_len]]))
        corrarray[method2][method1] =str('  '.join([str(round(x,3)) for x in [pearsonr(np.asarray(val)[:,0], np.asarray(val)[:,1])[0], val_len, unique_mutation_len]]))
        corrarray_spearman[method2][method1] =str('  '.join([str(round(x,3)) for x in [spearmanr(np.asarray(val)[:,0], np.asarray(val)[:,1])[0], val_len, unique_mutation_len]]))

        corrarray_spearman[method1][method2] =str('  '.join([str(round(x,3)) for x in [spearmanr(np.asarray(val)[:,0], np.asarray(val)[:,1])[0], val_len, unique_mutation_len]]))

        if True:
            fig = plt.figure()

            ax=fig.add_subplot(1, 1, 1)
            total= np.asarray(val)[:, 0].tolist() + np.asarray(val)[:, 1].tolist()
            max_v = max(total)
            min_v = min(total)

            ax.set(xlabel=method1, ylabel=method2)

            ax.scatter(np.asarray(val)[:, 0], np.asarray(val)[:, 1])
            #ax.axis('equal')
            ax.set_xlim(left=min_v, right=max_v)
            ax.set_ylim(bottom=min_v, top=max_v)


            ax.set_title('Axis [0, 0]')
            plt.savefig(f"results/correlation_ddg_no_nonbinding/{method1}_{method2}_correlation_fig.pdf", bbox_inches='tight', dpi=100)
            plt.close()
        index = index + 1
    else:
        corrarray[method1][method2] ="NA"#str('/'.join([str(x) for x in ["NA", val_len,unique_mutation_len]]))

corrarray.to_csv("data/processed/delta_delta_G_no_nonbinding_pearson_corrarray.csv")

corrarray_spearman.to_csv("data/processed/delta_delta_G_no_nonbinding_spearman_corrarray.csv")


# --- end ddG agreement analysis --------------------------------------

###############################################################################
# After building 'affinity_difference_df' in your script:
#   Rows = unique (pdb_mutations),
#   Columns = different methods,
#   Values = numeric affinity differences
###############################################################################

import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

# Suppose the user picks "SPR" as the reference method (edit as desired).
REFERENCE_METHOD = "SPR"

# Optionally, do a log transform if your data vary over orders of magnitude:
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

REFERENCE_METHOD = "SPR"

# (1) Suppose we've done log10 transform:

df_calibrated = pd.DataFrame(index=affinity_difference_df.index,
                             columns=affinity_difference_df.columns,
                             dtype=float)

# The reference remains the same
df_calibrated[REFERENCE_METHOD] = affinity_difference_df[REFERENCE_METHOD]

for method in df_calibrated.columns:
    if method == REFERENCE_METHOD:
        continue

    overlap_mask = (~affinity_difference_df[method].isna()) & \
                   (~affinity_difference_df[REFERENCE_METHOD].isna())
    overlap_count = overlap_mask.sum()
    if overlap_count < 3:
        print(f"Skipping '{method}', only {overlap_count} overlap with {REFERENCE_METHOD}.")
        df_calibrated[method] = np.nan
        continue

    # **Swap the roles**:
    # X = method’s raw values
    # y = reference’s values
    x_vals = affinity_difference_df.loc[overlap_mask, method].values.reshape(-1,1)
    y_vals = affinity_difference_df.loc[overlap_mask, REFERENCE_METHOD].values

    reg = LinearRegression()
    reg.fit(x_vals, y_vals)
    alpha, beta = reg.intercept_, reg.coef_[0]
    print(f"Calibrating {method} -> {REFERENCE_METHOD}:  Reference = {alpha:.3f} + {beta:.3f} * {method}")

    # Now create new_col = alpha + beta * old_method_value
    new_col = []
    for row_i in affinity_difference_df.index:
        method_val = affinity_difference_df.loc[row_i, method]
        if pd.notna(method_val):
            new_val = alpha + beta * method_val
            new_col.append(new_val)
        else:
            new_col.append(np.nan)

    df_calibrated[method] = new_col

df_calibrated.to_csv("data/processed/affinity_calibrated.csv")

# Re-check correlation on the calibrated DataFrame
corr_calibrated = df_calibrated.corr()
corr_calibrated.to_csv("data/processed/calibrated_corr.csv")

print("Multi-method calibration done. Calibrated data saved to 'affinity_calibrated.csv'.")
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

REFERENCE_METHOD = "SPR"  # or whichever method is your reference

# df: original data
# df_calibrated: same shape, new "reference-like" scale

for method in affinity_difference_df.columns:
    if method == REFERENCE_METHOD or all([np.isnan(x) for x in affinity_difference_df[method].values]):
        continue  # skip the reference itself

    # 1) Overlap for "before" plot
    # require at least some non-NaN after-calibration values
    if df_calibrated[method].notna().sum() < 5:
        continue

    mask_before = (~affinity_difference_df[method].isna()) & \
                (~affinity_difference_df[REFERENCE_METHOD].isna())
    mask_after  = (~df_calibrated[method].isna()) & \
                (~df_calibrated[REFERENCE_METHOD].isna())

    # compare on the common set to make r-before/after comparable
    mask = mask_before & mask_after

    print("Before overlap:", mask_before.sum())
    print("After  overlap:",  mask_after.sum())

    if mask.sum() < 3:
        print(f"No data to plot for method '{method}', skipping.")
        continue

    x_before = affinity_difference_df.loc[mask, method].values
    x_after  = df_calibrated.loc[mask, method].values
    y_ref    = df_calibrated.loc[mask, REFERENCE_METHOD].values  # same values as original for ref

    # correlations
    r_before, _ = pearsonr(x_before, y_ref) if x_before.size > 2 else (np.nan, None)
    r_after,  _ = pearsonr(x_after,  y_ref) if x_after.size  > 2 else (np.nan, None)








    # 4) Create a figure with 2 side-by-side subplots
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), facecolor="white")
    fig.suptitle(f"Comparison for Method {method} vs. {REFERENCE_METHOD}")

    # -- Left subplot: original data
    axes[0].scatter(x_before, y_ref, alpha=0.7, edgecolors='k')
    axes[0].set_title(f"Before Calibration\nPearson r={r_before:.2f}")
    axes[0].set_xlabel(f"{method} (original)")
    axes[0].set_ylabel(f"{REFERENCE_METHOD} (original)")

    # -- Right subplot: calibrated data
    axes[1].scatter(x_after, y_ref, alpha=0.7, edgecolors='k')
    axes[1].set_title(f"After Calibration\nPearson r={r_after:.2f}")
    axes[1].set_xlabel(f"{method} (calibrated)")
    axes[1].set_ylabel(f"{REFERENCE_METHOD} (reference)")

    plt.tight_layout(rect=[0, 0, 1, 0.95])  # leave room for suptitle
    plt.savefig(f"results/calibration/{_safe(method)}_vs_{_safe(REFERENCE_METHOD)}_before_after.pdf", dpi=120)

    plt.close(fig)

print("Created 'before and after' plots for all non-reference methods.")




vals_before = []
selected_methods=[]
for method in method_names:
    if sum(~df_calibrated[method].isna())<5:
        continue
    if method in df_calibrated.columns:
        selected_methods.append(method)
        # Extract the 'Delta_G_mut ' values for that method, drop NaNs
        arr = affinity_difference_df[method].dropna().values
        vals_before.append(arr)

# 2) Gather 'after' data from your calibrated DataFrame
#    Assuming 'df_calibrated' has columns = method_names, rows = mutations
vals_after = []
for method in method_names:
    if sum(~df_calibrated[method].isna())<5:
        continue
    if method in df_calibrated.columns:
        arr = df_calibrated[method].dropna().values
        vals_after.append(arr)


# 3) Create side-by-side boxplots
fig, axes = plt.subplots(ncols=2, figsize=(12, 6), sharey=True)

### LEFT PLOT: BEFORE calibration
axes[0].boxplot(vals_before, labels=selected_methods)
axes[0].set_title("Before Calibration")
axes[0].set_xlabel("Methods")
axes[0].set_ylabel("affinity_fold_change")
# Scatter the same points over each box
for i, arr in enumerate(vals_before, start=1):
    axes[0].scatter([i]*len(arr), arr, alpha=0.4)
plt.setp(axes[0].get_xticklabels(), rotation=45, ha="right")

### RIGHT PLOT: AFTER calibration
axes[1].boxplot(vals_after, labels=selected_methods)
axes[1].set_title("After Calibration")
axes[1].set_xlabel("Methods")
axes[1].set_ylabel("Calibrated affinity_fold_change")
plt.setp(axes[1].get_xticklabels(), rotation=45, ha="right")

# Scatter the same points over each box
for i, arr in enumerate(vals_after, start=1):
    axes[1].scatter([i]*len(arr), arr, alpha=0.4)

plt.tight_layout()
plt.savefig("results/summary/comparison_before_after_calibration_boxplot.pdf", dpi=120)
# PNG alongside the PDF: GitHub cannot render PDFs in README.md.
plt.savefig("results/summary/comparison_before_after_calibration_boxplot.png",
            dpi=150, bbox_inches="tight")
