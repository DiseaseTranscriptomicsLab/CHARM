#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
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
    "~/Projects/CHARM/intermediary_files/FullSimilarity_oddsratinc_nonorm.csv"
)

OUT_DEC = os.path.expanduser(
    "~/Projects/CHARM/intermediary_files/FullSimilarity_oddsratdec_nonorm.csv"
)

NPROC = max(1, cpu_count() - 4)

# =====================================================
# LOAD DATA
# =====================================================
print("Loading oddsrat dataframe...")
df = pd.read_csv(INPUT_CSV, sep="\t")

pos_cols = [c for c in df.columns if c.startswith("pos_")]

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
# PREPARE MATRICES
# =====================================================
def prepare_matrix(df):
    mat = df[pos_cols].to_numpy(dtype=np.float32)
    meta = df[["RBP", "Target", "dPSI"]].to_numpy()
    return mat, meta

inc_mat, inc_meta = prepare_matrix(df_inc)
dec_mat, dec_meta = prepare_matrix(df_dec)

# =====================================================
# METRICS
# =====================================================
def cosine_similarity(x, y):
    denom = np.linalg.norm(x) * np.linalg.norm(y)
    if denom == 0:
        return np.nan
    return np.dot(x, y) / denom

def jsd_similarity(p, q):
    # Normalize to probability distributions
    p_sum = p.sum()
    q_sum = q.sum()
    if p_sum == 0 or q_sum == 0:
        return np.nan

    p = p / p_sum
    q = q / q_sum

    m = 0.5 * (p + q)

    def kl(a, b):
        mask = (a > 0) & (b > 0)
        return np.sum(a[mask] * np.log2(a[mask] / b[mask]))

    jsd = 0.5 * kl(p, m) + 0.5 * kl(q, m)
    return 1 - jsd  # convert to similarity

def emd_similarity(p, q):
    # Normalize distributions
    p_sum = p.sum()
    q_sum = q.sum()
    if p_sum == 0 or q_sum == 0:
        return np.nan

    p = p / p_sum
    q = q / q_sum

    cdf_p = np.cumsum(p)
    cdf_q = np.cumsum(q)

    emd = np.sum(np.abs(cdf_p - cdf_q))
    max_dist = len(p) - 1
    emd_norm = emd / max_dist
    return 1 - emd_norm

# =====================================================
# WORKER FUNCTION
# =====================================================
def similarity_worker(args):
    i, j, mat, meta = args
    x = mat[i]
    y = mat[j]

    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 10:
        return None

    x = x[mask]
    y = y[mask]

    # Pearson
    r = np.corrcoef(x, y)[0, 1]
    if np.isnan(r):
        return None

    # Cosine
    cos = cosine_similarity(x, y)

    # Shift negative values for distribution metrics
    shift = min(np.min(x), np.min(y))
    x_shifted = x - shift + 1e-12
    y_shifted = y - shift + 1e-12

    # JSD similarity
    jsd_sim = jsd_similarity(x_shifted, y_shifted)

    # EMD similarity
    emd_sim = emd_similarity(x_shifted, y_shifted)

    # Combined similarity: geometric mean of Cosine and EMD similarity
    if np.isnan(cos) or np.isnan(emd_sim):
        combined = np.nan
    else:
        combined = np.sqrt(max(cos, 0) * max(emd_sim, 0))

    return {
        "RBP1": meta[i][0],
        "Target1": meta[i][1],
        "dPSI1": meta[i][2],
        "RBP2": meta[j][0],
        "Target2": meta[j][1],
        "dPSI2": meta[j][2],
        "PearsonR": r,
        "Cosine": cos,
        "JSD_sim": jsd_sim,
        "EMD_sim": emd_sim,
        "Combined_CosEMD": combined
    }

# =====================================================
# RUN ALL-PAIR COMPARISONS
# =====================================================
def run_all_pairs(mat, meta, outfile):
    n = mat.shape[0]
    pairs = combinations(range(n), 2)

    print(f"Running {n*(n-1)//2:,} comparisons ? {outfile}")

    with open(outfile, "w") as f:
        f.write(
            "RBP1,Target1,dPSI1,"
            "RBP2,Target2,dPSI2,"
            "PearsonR,Cosine,JSD_sim,EMD_sim,Combined_CosEMD\n"
        )

        with Pool(NPROC) as pool:
            for res in pool.imap_unordered(
                similarity_worker,
                ((i, j, mat, meta) for i, j in pairs),
                chunksize=500
            ):
                if res is None:
                    continue

                f.write(
                    f"{res['RBP1']},{res['Target1']},{res['dPSI1']},"
                    f"{res['RBP2']},{res['Target2']},{res['dPSI2']},"
                    f"{res['PearsonR']},{res['Cosine']},"
                    f"{res['JSD_sim']},{res['EMD_sim']},{res['Combined_CosEMD']}\n"
                )

# =====================================================
# EXECUTION
# =====================================================
if __name__ == "__main__":
    print(f"Using {NPROC} cores")
    run_all_pairs(inc_mat, inc_meta, OUT_INC)
    run_all_pairs(dec_mat, dec_meta, OUT_DEC)
    print("All similarity metrics completed.")
