#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings("ignore")

# =====================================================
# LOAD eCLIPSE BigExon file
# =====================================================
ECLIPSE_PATH = os.path.expanduser("~/Projects/CHARM/data/eCLIPSE_BigExon_HEPG2.txt")

eCLIPSE_BigExon_full = pd.read_csv(ECLIPSE_PATH, sep=" ")
rbps = eCLIPSE_BigExon_full["RBP"].unique()

print(f"Loaded eCLIPSE file with {len(eCLIPSE_BigExon_full)} rows and {len(eCLIPSE_BigExon_full.columns)} columns.")
print(f"Found {len(rbps)} unique RBPs.\n")


# =====================================================
# eCLIPSE_heatmap equivalent function
# =====================================================
def eCLIPSE_heatmap(rnamapfile, ASfile, rnaBP, PSIthreshold=0.05):
    """
    Computes per-position chi-square enrichment and odds ratios between 
    inclusion/exclusion (ASfile) and eCLIP binding (rnamapfile) for a given RBP.
    Now NaN-safe and division-by-zero-safe.
    """

    # ---------------------------------------------------------------
    # 1. Filter relevant AS events
    # ---------------------------------------------------------------
    ASfile = ASfile[ASfile["Event.ID"].str.contains("HsaEX", na=False)]
    if ASfile.empty:
        return None

    ASfile_filt_increased = ASfile[ASfile["dPSI"] > PSIthreshold]
    ASfile_filt_decreased = ASfile[ASfile["dPSI"] < -PSIthreshold]
    ASfile_filt_maintained = ASfile[
        (~(ASfile["dPSI"] > PSIthreshold)) & (~(ASfile["dPSI"] < -PSIthreshold))
    ]

    def subset_events(df, events):
        return df[df["EVENTS"].isin(events)]

    # ---------------------------------------------------------------
    # 2. Split into increased / decreased / maintained subsets
    # ---------------------------------------------------------------
    inc = subset_events(rnamapfile, ASfile_filt_increased["Event.ID"])
    dec = subset_events(rnamapfile, ASfile_filt_decreased["Event.ID"])
    main = subset_events(rnamapfile, ASfile_filt_maintained["Event.ID"])

    # ---------------------------------------------------------------
    # 3. Separate target RBP vs non-targets
    # ---------------------------------------------------------------
    def split_rbp(df, rnaBP):
        rbp_df = df[df["RBP"] == rnaBP]
        nrbp_df = df[df["RBP"] != rnaBP].drop_duplicates(subset=["EVENTS"]).replace(1, 0)
        return pd.concat([rbp_df, nrbp_df], axis=0)

    inc = split_rbp(inc, rnaBP)
    dec = split_rbp(dec, rnaBP)
    main = split_rbp(main, rnaBP)

    # Skip if not enough data
    if len(inc) < 10 or len(dec) < 10 or len(main) < 10:
        return None

    # ---------------------------------------------------------------
    # 4. Prepare numeric matrices (ignore last 2 metadata columns)
    # ---------------------------------------------------------------
    def prepare_matrix(df):
        mat = df.iloc[:, :-2].apply(pd.to_numeric, errors="coerce").to_numpy()
        comp = np.nansum(mat, axis=0)
        incomp = np.sum(~np.isnan(mat), axis=0)
        return mat, comp, incomp

    inc_mat, inc_comp, inc_incomp = prepare_matrix(inc)
    dec_mat, dec_comp, dec_incomp = prepare_matrix(dec)
    main_mat, main_comp, main_incomp = prepare_matrix(main)

    # ---------------------------------------------------------------
    # 4b. Compute normalized raw binding profiles (like R version)
    # ---------------------------------------------------------------
    def safe_div(a, b):
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(b != 0, a / b, np.nan)

    increased_norm = safe_div(inc_comp + 1, inc_incomp)
    decreased_norm = safe_div(dec_comp + 1, dec_incomp)
    maintained_norm = safe_div(main_comp + 1, main_incomp)

    df_raw = pd.DataFrame({
        "pos": np.arange(1, 1001),
        "IncreasedEvents": increased_norm,
        "DecreasedEvents": decreased_norm,
        "MaintainedEvents": maintained_norm
    })

    positions = np.arange(1, 1001)

    pval_inc, pval_dec, odds_inc, odds_dec = [], [], [], []

    # ---------------------------------------------------------------
    # 5. Safe helper functions
    # ---------------------------------------------------------------
    def safe_ratio(a, b):
        return a / b if b != 0 else np.nan

    def safe_chi2(table):
        # if any cell = 0, skip
        if np.any(table == 0):
            return np.nan, np.nan
        try:
            chi2, p, _, _ = stats.chi2_contingency(table)
            return chi2, p
        except Exception:
            return np.nan, np.nan

    # ---------------------------------------------------------------
    # 6. Loop through positions
    # ---------------------------------------------------------------
    for i in range(1000):
        inc_fish = inc_comp[i] + 1
        dec_fish = dec_comp[i] + 1
        main_fish = main_comp[i] + 1

        inc_fish_rev = inc_incomp[i] - inc_comp[i]
        dec_fish_rev = dec_incomp[i] - dec_comp[i]
        main_fish_rev = main_incomp[i] - main_comp[i]

        inc_table = np.array([[inc_fish, inc_fish_rev], [main_fish, main_fish_rev]])
        dec_table = np.array([[dec_fish, dec_fish_rev], [main_fish, main_fish_rev]])

        inc_chi2, inc_p = safe_chi2(inc_table)
        dec_chi2, dec_p = safe_chi2(dec_table)

        pval_inc.append(inc_p)
        pval_dec.append(dec_p)
        odds_inc.append(inc_chi2)
        odds_dec.append(dec_chi2)

    # ---------------------------------------------------------------
    # 7. BH correction with NaN handling
    # ---------------------------------------------------------------
    def safe_multipletests(pvals):
        pvals = np.array(pvals, dtype=float)
        mask = ~np.isnan(pvals)
        corrected = np.full_like(pvals, np.nan)
        if mask.sum() > 0:
            _, p_corr, _, _ = multipletests(pvals[mask], method="fdr_bh")
            corrected[mask] = p_corr
        return -np.log10(corrected)

    pval_inc_corr = safe_multipletests(pval_inc)
    pval_dec_corr = safe_multipletests(pval_dec)

    # ---------------------------------------------------------------
    # 8. Sign flipping (with safe ratios)
    # ---------------------------------------------------------------
    for i in range(1000):
        inc_ratio = safe_ratio(inc_comp[i], inc_incomp[i])
        main_ratio = safe_ratio(main_comp[i], main_incomp[i])
        dec_ratio = safe_ratio(dec_comp[i], dec_incomp[i])

        if np.isfinite(inc_ratio) and np.isfinite(main_ratio) and inc_ratio <= main_ratio:
            pval_inc_corr[i] = -abs(pval_inc_corr[i]) if not np.isnan(pval_inc_corr[i]) else np.nan
            odds_inc[i] = -abs(odds_inc[i]) if not np.isnan(odds_inc[i]) else np.nan

        if np.isfinite(dec_ratio) and np.isfinite(main_ratio) and dec_ratio <= main_ratio:
            pval_dec_corr[i] = -abs(pval_dec_corr[i]) if not np.isnan(pval_dec_corr[i]) else np.nan
            odds_dec[i] = -abs(odds_dec[i]) if not np.isnan(odds_dec[i]) else np.nan

    # ---------------------------------------------------------------
    # 9. Combine results into DataFrames
    # ---------------------------------------------------------------
    df_pval = pd.DataFrame({
        "pos": positions,
        "pvalinc": pval_inc_corr,
        "pvaldec": pval_dec_corr
    })

    df_odds = pd.DataFrame({
        "pos": positions,
        "oddsratinc": odds_inc,
        "oddsratdec": odds_dec
    })

    return {"pval_data": df_pval, "oddsratio_data": df_odds, "raw_data": df_raw}


# =====================================================
# Worker function for multiprocessing
# =====================================================
def process_rbp(rbp):
    base_path = os.path.expanduser("~/Projects/StressGranules/AS.WC_Transcriptome/shRNAExp")
    file_path = os.path.join(base_path, rbp, f"{rbp}_VulcanTable_HEPG2.txt")

    if not os.path.exists(file_path):
        print(f"⚠️ Skipping missing file for {rbp}")
        return None

    print(f"Processing {rbp} ...")
    ASfile = pd.read_csv(file_path, sep="\t", header=0)
    results_all = []

    targets = eCLIPSE_BigExon_full["RBP"].unique()
    dpsi_values = np.arange(0.01, 0.11, 0.01)

    for target in targets:
        temp_results = {}
        for dpsi in dpsi_values:
            try:
                res = eCLIPSE_heatmap(eCLIPSE_BigExon_full, ASfile, target, PSIthreshold=dpsi)
                if res is not None:
                    temp_results[round(dpsi, 2)] = res
            except Exception as e:
                print(f"Error for {rbp}-{target} dPSI={dpsi}: {e}")
                continue

        if not temp_results:
            continue

        # Summarize metrics
        stats_df = pd.DataFrame({
            "dpsi": list(temp_results.keys()),
            "sum_pvalinc": [res["pval_data"]["pvalinc"].sum() for res in temp_results.values()],
            "sum_pvaldec": [res["pval_data"]["pvaldec"].sum() for res in temp_results.values()],
            "sum_oddsinc": [res["oddsratio_data"]["oddsratinc"].sum() for res in temp_results.values()],
            "sum_oddsdec": [res["oddsratio_data"]["oddsratdec"].sum() for res in temp_results.values()],
        })

        # Use absolute values for maximization
        best_pvalinc = stats_df.loc[stats_df["sum_pvalinc"].abs().idxmax(), "dpsi"]
        best_pvaldec = stats_df.loc[stats_df["sum_pvaldec"].abs().idxmax(), "dpsi"]
        best_oddsinc = stats_df.loc[stats_df["sum_oddsinc"].abs().idxmax(), "dpsi"]
        best_oddsdec = stats_df.loc[stats_df["sum_oddsdec"].abs().idxmax(), "dpsi"]

        keep_dpsis_info = {
            0.05: "default",
            0.1: "default",
            float(best_pvalinc): "maximised pvalinc",
            float(best_pvaldec): "maximised pvaldec",
            float(best_oddsinc): "maximised oddsratinc",
            float(best_oddsdec): "maximised oddsratdec"
        }

        keep_dpsis = np.unique(list(keep_dpsis_info.keys()))
        print(f"  Target {target}: keeping dPSIs {[f'{v:.2f} ({keep_dpsis_info[v]})' for v in keep_dpsis]}")

        for dpsi_sel in keep_dpsis:
            res = temp_results.get(round(dpsi_sel, 2))
            if res is None:
                continue

            label = f"{dpsi_sel:.2f} ({keep_dpsis_info.get(float(dpsi_sel), 'NA')})"

            df_long = pd.concat([
                res["pval_data"].assign(Metric="pvalinc", dPSI=label, Target=target)
                    .rename(columns={"pvalinc": "Value", "pos": "Pos"})[["Pos", "Metric", "Value", "dPSI", "Target"]],
                res["pval_data"].assign(Metric="pvaldec", dPSI=label, Target=target)
                    .rename(columns={"pvaldec": "Value", "pos": "Pos"})[["Pos", "Metric", "Value", "dPSI", "Target"]],
                res["oddsratio_data"].assign(Metric="oddsratinc", dPSI=label, Target=target)
                    .rename(columns={"oddsratinc": "Value", "pos": "Pos"})[["Pos", "Metric", "Value", "dPSI", "Target"]],
                res["oddsratio_data"].assign(Metric="oddsratdec", dPSI=label, Target=target)
                    .rename(columns={"oddsratdec": "Value", "pos": "Pos"})[["Pos", "Metric", "Value", "dPSI", "Target"]],
                res["raw_data"].melt(id_vars=["pos"], var_name="Metric", value_name="Value")
                    .assign(dPSI=label, Target=target)
                    .rename(columns={"pos": "Pos"})[["Pos", "Metric", "Value", "dPSI", "Target"]]
            ], ignore_index=True)

            results_all.append(df_long)

    if results_all:
        combined = pd.concat(results_all, ignore_index=True)
        out_path = os.path.join(base_path, rbp, f"{rbp}_BindingValues_HEPG2.txt")
        combined.to_csv(out_path, sep="\t", index=False)
        print(f"✅ Saved {len(combined)} rows for {rbp} → {out_path}")
    else:
        print(f"⚠️ No results to save for {rbp}")

    return rbp


# =====================================================
# Run multiprocessing
# =====================================================
if __name__ == "__main__":
    nproc = max(1, cpu_count() - 10)
    print(f"Starting multiprocessing with {nproc} workers ...\n")

    with Pool(nproc) as pool:
        pool.map(process_rbp, rbps)

    print("\n✅ All RBPs processed.")
