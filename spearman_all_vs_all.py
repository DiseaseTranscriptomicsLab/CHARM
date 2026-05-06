#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from multiprocessing import Pool, cpu_count
from itertools import combinations
import warnings
warnings.filterwarnings("ignore")

# =====================================================
# CONFIG
# =====================================================
INPUT_CSV = os.path.expanduser(
    "~/Projects/CHARM/intermediary_files/Charm.binding.filtered.both_new.csv"
)

OUT_INC = os.path.expanduser(
    "~/Projects/CHARM/intermediary_files/Spearman_oddsratinc_new.csv"
)

OUT_DEC = os.path.expanduser(
    "~/Projects/CHARM/intermediary_files/Spearman_oddsratdec_new.csv"
)

NPROC = max(1, cpu_count() - 4)

# =====================================================
# LOAD DATA
# =====================================================
print("Loading oddsrat dataframe...")
df = pd.read_csv(INPUT_CSV, sep="\t")

pos_cols = [c for c in df.columns if c.startswith("pos_")]

# Ensure numeric
df[pos_cols] = df[pos_cols].apply(pd.to_numeric, errors="coerce")

print(f"Loaded {len(df)} spectra with {len(pos_cols)} positions")

# =====================================================
# SPLIT INC / DEC
# =====================================================
df_inc = df[df["Direction"] == "oddsratinc"].reset_index(drop=True)
df_dec = df[df["Direction"] == "oddsratdec"].reset_index(drop=True)

print(f"INC spectra: {len(df_inc)}")
print(f"DEC spectra: {len(df_dec)}")

# =====================================================
# PREPARE NUMPY MATRICES
# =====================================================
def prepare_matrix(df):
    mat = df[pos_cols].to_numpy(dtype=np.float32)
    meta = df[["RBP", "Target", "dPSI"]].to_numpy()
    return mat, meta

inc_mat, inc_meta = prepare_matrix(df_inc)
dec_mat, dec_meta = prepare_matrix(df_dec)

# =====================================================
# WORKER FUNCTION (SPEARMAN)
# =====================================================
def spearman_worker(args):
    i, j, mat, meta = args
    x = mat[i]
    y = mat[j]

    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 10:
        return None

    rho, _ = spearmanr(x[mask], y[mask])
    if np.isnan(rho):
        return None

    return {
        "RBP1": meta[i][0],
        "Target1": meta[i][1],
        "dPSI1": meta[i][2],
        "RBP2": meta[j][0],
        "Target2": meta[j][1],
        "dPSI2": meta[j][2],
        "SpearmanR": rho
    }

# =====================================================
# RUN ALL-PAIR COMPARISONS
# =====================================================
def run_all_pairs(mat, meta, outfile):
    n = mat.shape[0]
    pairs = combinations(range(n), 2)

    print(f"Running {n*(n-1)//2:,} comparisons ? {outfile}")

    with open(outfile, "w") as f:
        f.write("RBP1,Target1,dPSI1,RBP2,Target2,dPSI2,SpearmanR\n")

        with Pool(NPROC) as pool:
            for res in pool.imap_unordered(
                spearman_worker,
                ((i, j, mat, meta) for i, j in pairs),
                chunksize=500
            ):
                if res is None:
                    continue
                f.write(
                    f"{res['RBP1']},{res['Target1']},{res['dPSI1']},"
                    f"{res['RBP2']},{res['Target2']},{res['dPSI2']},"
                    f"{res['SpearmanR']}\n"
                )

# =====================================================
# EXECUTION
# =====================================================
if __name__ == "__main__":
    print(f"Using {NPROC} cores")

    run_all_pairs(inc_mat, inc_meta, OUT_INC)
    run_all_pairs(dec_mat, dec_meta, OUT_DEC)

    print("? All Spearman correlations completed.")
