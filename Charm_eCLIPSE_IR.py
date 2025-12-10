#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool
import warnings
warnings.filterwarnings("ignore")

# === Load Intron file ===
print("Loading eCLIPSE_Intron_full ...")
ECLIPSE_INTRON_PATH = os.path.expanduser("~/Projects/CHARM/data/eCLIPSE_Intron_HEPG2.txt")
eCLIPSE_Intron_full = pd.read_csv(ECLIPSE_INTRON_PATH, sep=" ")
print(f"✅ Loaded {eCLIPSE_Intron_full.shape[0]} rows, {eCLIPSE_Intron_full.shape[1]} columns.")


# =====================================================
# eCLIPSE_heatmap_IR (safe + normalized raw values + event counts)
# =====================================================
def eCLIPSE_heatmap_IR(metric, rnamapfile, pvalthreshold, PSIthreshold, ASfile, rnaBP, plot=False):
    """
    Safe version of eCLIPSE_heatmap_IR.
    Handles zeros, NaNs, and invalid contingency tables gracefully.
    Returns pval_data, oddsratio_data, raw_data, counts_data.
    """

    # 1. Filter intron retention events
    ASfile = ASfile[ASfile["Event.ID"].str.contains("HsaINT", na=False)]
    if ASfile.empty:
        return None

    ASfile_filt_increased = ASfile[ASfile["dPSI"] > PSIthreshold]
    ASfile_filt_decreased = ASfile[ASfile["dPSI"] < -PSIthreshold]
    ASfile_filt_maintained = ASfile[
        (~(ASfile["dPSI"] > PSIthreshold)) & (~(ASfile["dPSI"] < -PSIthreshold))
    ]

    # 2. Subset RNA map by events
    inc = rnamapfile[rnamapfile["EVENTS"].isin(ASfile_filt_increased["Event.ID"])]
    dec = rnamapfile[rnamapfile["EVENTS"].isin(ASfile_filt_decreased["Event.ID"])]
    mai = rnamapfile[rnamapfile["EVENTS"].isin(ASfile_filt_maintained["Event.ID"])]

    # 3. Separate target RBP vs non-RBP
    def prepare_subsets(df, rnaBP):
        rbp = df[df["RBP"] == rnaBP]
        nrbp = df[df["RBP"] != rnaBP].drop_duplicates(subset="EVENTS").copy()
        nrbp.replace(1, 0, inplace=True)
        return pd.concat([rbp, nrbp])

    inc = prepare_subsets(inc, rnaBP)
    dec = prepare_subsets(dec, rnaBP)
    mai = prepare_subsets(mai, rnaBP)

    if any(len(x) < 10 for x in [inc, dec, mai]):
        return None

    # 4. Convert numeric data
    def to_numeric(df):
        return df.iloc[:, :-2].apply(pd.to_numeric, errors="coerce")

    inc_num = to_numeric(inc)
    dec_num = to_numeric(dec)
    mai_num = to_numeric(mai)

    # 5. Compute counts and completeness
    inc_comp = np.nansum(inc_num.values, axis=0)
    dec_comp = np.nansum(dec_num.values, axis=0)
    mai_comp = np.nansum(mai_num.values, axis=0)

    inc_incomp = np.sum(~np.isnan(inc_num.values), axis=0)
    dec_incomp = np.sum(~np.isnan(dec_num.values), axis=0)
    mai_incomp = np.sum(~np.isnan(mai_num.values), axis=0)

    # 6. Normalized raw values
    def safe_div(a, b):
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(b != 0, a / b, np.nan)

    increased_norm = safe_div(inc_comp + 1, inc_incomp)
    decreased_norm = safe_div(dec_comp + 1, dec_incomp)
    maintained_norm = safe_div(mai_comp + 1, mai_incomp)

    df_raw = pd.DataFrame({
        "pos": np.arange(1, len(mai_comp) + 1),
        "IncreasedEvents": increased_norm,
        "DecreasedEvents": decreased_norm,
        "MaintainedEvents": maintained_norm
    })

    # 6b. Compute event counts
    inc_target_events = inc[inc["RBP"] == rnaBP]["EVENTS"].nunique()
    dec_target_events = dec[dec["RBP"] == rnaBP]["EVENTS"].nunique()
    mai_target_events = mai[mai["RBP"] == rnaBP]["EVENTS"].nunique()

    total_increased = ASfile_filt_increased.shape[0]
    total_decreased = ASfile_filt_decreased.shape[0]
    total_maintained = ASfile_filt_maintained.shape[0]

    df_counts = pd.DataFrame([
        {"Pos": np.nan, "Metric": "IncreasedTargetEvents", "Value": inc_target_events},
        {"Pos": np.nan, "Metric": "TotalIncreasedEvents", "Value": total_increased},
        {"Pos": np.nan, "Metric": "MaintainedTargetEvents", "Value": mai_target_events},
        {"Pos": np.nan, "Metric": "TotalMaintainedEvents", "Value": total_maintained},
        {"Pos": np.nan, "Metric": "DecreasedTargetEvents", "Value": dec_target_events},
        {"Pos": np.nan, "Metric": "TotalDecreasedEvents", "Value": total_decreased}
    ])

    # 7. Safe helpers
    def safe_chi2(table):
        if np.any(np.isnan(table)) or np.any(table <= 0):
            return np.nan, np.nan
        try:
            chi2, p, _, _ = chi2_contingency(table)
            return chi2, p
        except Exception:
            return np.nan, np.nan

    def safe_ratio(a, b):
        return a / b if b != 0 else np.nan

    # 8. Run chi-square tests across positions
    pval_inc, pval_dec, odds_inc, odds_dec = [], [], [], []
    n_positions = min(500, len(df_raw))

    for i in range(n_positions):
        inc_fish = inc_comp[i] + 1
        dec_fish = dec_comp[i] + 1
        mai_fish = mai_comp[i] + 1

        inc_table = np.array([[inc_fish, inc_incomp[i] - inc_fish],
                              [mai_fish, mai_incomp[i] - mai_fish]])
        dec_table = np.array([[dec_fish, dec_incomp[i] - dec_fish],
                              [mai_fish, mai_incomp[i] - mai_fish]])

        inc_chi, inc_p = safe_chi2(inc_table)
        dec_chi, dec_p = safe_chi2(dec_table)

        pval_inc.append(inc_p)
        pval_dec.append(dec_p)
        odds_inc.append(inc_chi)
        odds_dec.append(dec_chi)

    # 9. BH correction
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

    # 10. Sign flipping
    for i in range(n_positions):
        inc_ratio = safe_ratio(inc_comp[i], inc_incomp[i])
        main_ratio = safe_ratio(mai_comp[i], mai_incomp[i])
        dec_ratio = safe_ratio(dec_comp[i], dec_incomp[i])

        if np.isfinite(inc_ratio) and np.isfinite(main_ratio) and inc_ratio <= main_ratio:
            pval_inc_corr[i] = -abs(pval_inc_corr[i]) if not np.isnan(pval_inc_corr[i]) else np.nan
            odds_inc[i] = -abs(odds_inc[i]) if not np.isnan(odds_inc[i]) else np.nan
        if np.isfinite(dec_ratio) and np.isfinite(main_ratio) and dec_ratio <= main_ratio:
            pval_dec_corr[i] = -abs(pval_dec_corr[i]) if not np.isnan(pval_dec_corr[i]) else np.nan
            odds_dec[i] = -abs(odds_dec[i]) if not np.isnan(odds_dec[i]) else np.nan

    # 11. Build output DataFrames
    pval_data = pd.DataFrame({
        "pos": np.arange(1, n_positions + 1),
        "pvalinc": pval_inc_corr,
        "pvaldec": pval_dec_corr
    })

    oddsratio_data = pd.DataFrame({
        "pos": np.arange(1, n_positions + 1),
        "oddsratinc": odds_inc,
        "oddsratdec": odds_dec
    })

    raw_data = df_raw.melt(id_vars=["pos"], var_name="Metric", value_name="Value")

    return {"pval_data": pval_data, "oddsratio_data": oddsratio_data, "raw_data": raw_data, "counts_data": df_counts}


# =====================================================
# Worker function for multiprocessing
# =====================================================
def process_rbp(rbp):
    base_path = os.path.expanduser("~/Projects/StressGranules/AS.WC_Transcriptome/shRNAExp")
    file_path = os.path.join(base_path, rbp, f"{rbp}_VulcanTable_HEPG2.txt")

    if not os.path.exists(file_path):
        print(f"Skipping missing file for {rbp}")
        return None

    print(f"Processing {rbp} ...")
    betAS = pd.read_csv(file_path, sep="\t")
    all_metrics = []

    for target in eCLIPSE_Intron_full["RBP"].unique():
        dpsi_values = np.arange(0.01, 0.11, 0.01)
        temp_results = {}

        for dpsi in dpsi_values:
            try:
                res = eCLIPSE_heatmap_IR(
                    metric=None,
                    rnamapfile=eCLIPSE_Intron_full,
                    pvalthreshold=0.05,
                    PSIthreshold=dpsi,
                    ASfile=betAS,
                    rnaBP=target
                )
                if res is not None:
                    temp_results[round(dpsi, 3)] = res
            except Exception as e:
                print(f"    Error for {target} dPSI={dpsi}: {e}")

        if not temp_results:
            continue

        # summarize metrics
        stats = pd.DataFrame({
            "dpsi": list(temp_results.keys()),
            "sum_pvalinc": [res["pval_data"]["pvalinc"].sum() for res in temp_results.values()],
            "sum_pvaldec": [res["pval_data"]["pvaldec"].sum() for res in temp_results.values()],
            "sum_oddsinc": [res["oddsratio_data"]["oddsratinc"].sum() for res in temp_results.values()],
            "sum_oddsdec": [res["oddsratio_data"]["oddsratdec"].sum() for res in temp_results.values()]
        })

        best_pvalinc = stats.loc[stats["sum_pvalinc"].abs().idxmax(), "dpsi"]
        best_pvaldec = stats.loc[stats["sum_pvaldec"].abs().idxmax(), "dpsi"]
        best_oddsinc = stats.loc[stats["sum_oddsinc"].abs().idxmax(), "dpsi"]
        best_oddsdec = stats.loc[stats["sum_oddsdec"].abs().idxmax(), "dpsi"]

        keep_dpsis_info = {
            0.05: "default",
            0.1: "default",
            float(best_pvalinc): "maximised pvalinc",
            float(best_pvaldec): "maximised pvaldec",
            float(best_oddsinc): "maximised oddsratinc",
            float(best_oddsdec): "maximised oddsratdec"
        }
        keep_dpsis = np.unique(list(keep_dpsis_info.keys()))
        print(f"  Target {target}: keeping {[f'{v:.2f} ({keep_dpsis_info[v]})' for v in keep_dpsis]}")

        for dpsi_sel in keep_dpsis:
            res = temp_results.get(round(dpsi_sel, 3))
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
                res["raw_data"].assign(dPSI=label, Target=target)
                    .rename(columns={"pos": "Pos"})[["Pos", "Metric", "Value", "dPSI", "Target"]],
                res["counts_data"].assign(dPSI=label, Target=target)[["Pos", "Metric", "Value", "dPSI", "Target"]]
            ], ignore_index=True)

            all_metrics.append(df_long)

    if all_metrics:
        combined_df = pd.concat(all_metrics, ignore_index=True)
        outfile = os.path.join(base_path, rbp, f"{rbp}_BindingValues_IR_HEPG2.txt")
        combined_df.to_csv(outfile, sep="\t", index=False)
        print(f"✅ Saved {len(combined_df)} rows for {rbp} → {outfile}")
    else:
        print(f"⚠️ No results to save for {rbp}")


# =====================================================
# Run multiprocessing
# =====================================================
if __name__ == "__main__":
    rbps = eCLIPSE_Intron_full["RBP"].unique()
    print(f"Running multiprocessing for {len(rbps)} RBPs ...")

    with Pool(os.cpu_count() - 10) as pool:
        pool.map(process_rbp, rbps)

    print("\n✅ All RBPs processed.")
