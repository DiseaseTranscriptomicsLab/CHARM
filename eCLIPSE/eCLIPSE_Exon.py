import argparse
import pandas as pd
import numpy as np
import time, os
from multiprocessing import Pool

starttime = time.time()

# -----------------------------------------------------------------------------
# CLI arguments
#
# Originally this script hardcoded its two input filenames:
#   ./Final_BIGEX_Table.txt   (the splicing-genome / event-region table)
#   ./BenFile_Processed.txt   (the eCLIP peaks table)
# and its output filename:
#   ./BenFile_ICLIP.txt
#
# Both preprocessing paths in this folder produce differently-named event
# tables (Eclip_preprocessing_vttools.R writes Final_MIC_Table.txt,
# Eclip_preprocessing_rmats.R writes Final_AS_Table_Rmats.txt), and the raw
# combined eCLIP peaks file is named all_data_combined_38.txt - none of these
# match the old hardcoded names. Rather than requiring you to rename files
# by hand each time, the three filenames are now CLI flags (defaults match
# the original hardcoded names, so old invocations without flags still work).
# See the ECLIPSE README for exactly which file to pass for each pipeline
# path (VastDB/VAST-TOOLS vs. rMATS).
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Intersect eCLIP peaks against exon-skipping splicing-genome regions "
                "(the 'binding metagenome assembly' step of the eCLIPSE pipeline)."
)
parser.add_argument("--events", default="./Final_BIGEX_Table.txt",
                     help="Splicing-genome / event-region table produced by "
                          "Eclip_preprocessing_vttools.R (Final_MIC_Table.txt) or "
                          "Eclip_preprocessing_rmats.R (Final_AS_Table_Rmats.txt). "
                          "Default: ./Final_BIGEX_Table.txt")
parser.add_argument("--peaks", default="./BenFile_Processed.txt",
                     help="eCLIP peaks table, space-delimited with columns "
                          "Chrom/StartCord/EndCord/RBP (e.g. all_data_combined_38.txt). "
                          "Default: ./BenFile_Processed.txt")
parser.add_argument("--output", default="./BenFile_ICLIP.txt",
                     help="Output path for the per-RBP, per-event overlap table. "
                          "Default: ./BenFile_ICLIP.txt")
args = parser.parse_args()

#Establishing the file that will be used

map = pd.read_csv(args.events)
map = map[map['GENE'].astype(str) != ""]
NameGiver = map.dropna()

all_data_combined = pd.read_csv(args.peaks, sep=" ")

# This is the definition of the scaling function, so assure that RNA binding maps are possible to do.
def find_common_available_range_scaled_ups(range1_start, range1_end, range2_start, range2_end):
    common_start = max(range1_start, range2_start)
    common_end = min(range1_end, range2_end)
    if common_start <= common_end:
        scale_factor = range1_start
        common_start_scaled = common_start - scale_factor + 1
        common_end_scaled = common_end - scale_factor
        return (common_start, common_end), (common_start_scaled, common_end_scaled)
    else:
        return None, None

def find_common_available_range_scaled_dow(range1_start, range1_end, range2_start, range2_end):
    common_start = max(range1_start, range2_start)
    common_end = min(range1_end, range2_end)
    if common_start <= common_end:
        scale_factor = range1_end
        common_start_scaled = scale_factor - common_start
        common_end_scaled = scale_factor - common_end + 1
        return (common_start, common_end), (common_start_scaled, common_end_scaled)
    else:
        return None, None


# --- NEW: precompute per-chromosome NumPy arrays ONCE instead of filtering NameGiver
# for every single one of the 9.25M eCLIP rows. On Linux, Pool uses fork by default,
# so this is built once in the parent and inherited (copy-on-write) by workers for free.

META_COLS = ['GENE', 'EVENT', 'STRAND']
REGION_COLS = [
    'UPS.EX.ST', 'EX1',
    'UPS.INT.ST.ST', 'UPS.INT.ST.EN',
    'UPS.INT.EN.ST', 'UPS.INT.EN.EN',
    'AE.ST', 'AER.ST',
    'AER.EN', 'AE.EN',
    'DOW.INT.ST.ST', 'DOW.INT.ST.EN',
    'DOW.INT.EN.ST', 'DOW.INT.EN.EN',
    'EX2', 'DOW.EX.EN'
]

NameGiver_by_chrom = {}
for chrom, df in NameGiver.groupby('CHROM'):
    NameGiver_by_chrom[chrom] = {
        'meta': df[META_COLS].to_numpy(dtype=object),
        'coords': df[REGION_COLS].to_numpy(dtype=float),
    }


# This is the function that calculates the binding (overlap of coordinates) of the RBPs
# to the regulatory regions defined. Specific for exonic coordinates.
def process_row(args):
    index_B, row_B = args  # NEW: no longer carries the full 500MB dataframe
    cromo = row_B['Chrom']

    arrs = NameGiver_by_chrom.get(cromo)
    rbp_name = row_B['RBP']

    if arrs is None:
        return rbp_name, index_B, []

    meta = arrs['meta']
    coords = arrs['coords']
    n = coords.shape[0]

    startB = row_B['StartCord']
    endB = row_B['EndCord']

    results = []
    # NEW: raw NumPy array access instead of DataFrame.iterrows() (much faster per-row)
    for i in range(n):
        (ups_ex_st, ex1,
         ups_in_st_st, ups_in_st_en,
         ups_in_en_st, ups_in_en_en,
         ae_st, aer_st,
         aer_en, ae_en,
         dow_in_st_st, dow_in_st_en,
         dow_in_en_st, dow_in_en_en,
         ex2, dow_ex_en) = coords[i]
        gene, event, strand = meta[i]

        overlap_UPS_EX_occurred = False
        overlap_UPS_IN_ST_occurred = False
        overlap_UPS_IN_EN_occurred = False
        overlap_AE_ST_occurred = False
        overlap_AE_EN_occurred = False
        overlap_DOW_IN_ST_occurred = False
        overlap_DOW_IN_EN_occurred = False
        overlap_DOW_EX_occurred = False

        overlap_UPS_EX = (pd.NA, pd.NA); scaled_values_UPS_EX = (pd.NA, pd.NA); length_UPS_EX = pd.NA
        overlap_UPS_IN_ST = (pd.NA, pd.NA); scaled_values_UPS_IN_ST = (pd.NA, pd.NA); length_UPS_IN_ST = pd.NA
        overlap_UPS_IN_EN = (pd.NA, pd.NA); scaled_values_UPS_IN_EN = (pd.NA, pd.NA); length_UPS_IN_EN = pd.NA
        overlap_AE_ST = (pd.NA, pd.NA); scaled_values_AE_ST = (pd.NA, pd.NA); length_AE_ST = pd.NA
        overlap_AE_EN = (pd.NA, pd.NA); scaled_values_AE_EN = (pd.NA, pd.NA); length_AE_EN = pd.NA
        overlap_DOW_IN_ST = (pd.NA, pd.NA); scaled_values_DOW_IN_ST = (pd.NA, pd.NA); length_DOW_IN_ST = pd.NA
        overlap_DOW_IN_EN = (pd.NA, pd.NA); scaled_values_DOW_IN_EN = (pd.NA, pd.NA); length_DOW_IN_EN = pd.NA
        overlap_DOW_EX = (pd.NA, pd.NA); scaled_values_DOW_EX = (pd.NA, pd.NA); length_DOW_EX = pd.NA

        if ups_ex_st <= endB and ex1 >= startB:
            overlap_UPS_EX, scaled_values_UPS_EX = find_common_available_range_scaled_ups(
                ups_ex_st, ex1, startB, endB
            )
            length_UPS_EX = ex1 - ups_ex_st
            overlap_UPS_EX_occurred = True

        if ups_in_st_st <= endB and ups_in_st_en >= startB:
            overlap_UPS_IN_ST, scaled_values_UPS_IN_ST = find_common_available_range_scaled_ups(
                ups_in_st_st, ups_in_st_en, startB, endB
            )
            length_UPS_IN_ST = ups_in_st_en - ups_in_st_st
            overlap_UPS_IN_ST_occurred = True

        if ups_in_en_st <= endB and ups_in_en_en >= startB:
            overlap_UPS_IN_EN, scaled_values_UPS_IN_EN = find_common_available_range_scaled_dow(
                ups_in_en_st, ups_in_en_en, startB, endB
            )
            length_UPS_IN_EN = ups_in_en_en - ups_in_en_st
            overlap_UPS_IN_EN_occurred = True

        if ae_st <= endB and aer_st >= startB:
            overlap_AE_ST, scaled_values_AE_ST = find_common_available_range_scaled_ups(
                ae_st, aer_st, startB, endB
            )
            length_AE_ST = aer_st - ae_st
            overlap_AE_ST_occurred = True

        if aer_en <= endB and ae_en >= startB:
            overlap_AE_EN, scaled_values_AE_EN = find_common_available_range_scaled_dow(
                aer_en, ae_en, startB, endB
            )
            length_AE_EN = ae_en - aer_en
            overlap_AE_EN_occurred = True

        if dow_in_st_st <= endB and dow_in_st_en >= startB:
            overlap_DOW_IN_ST, scaled_values_DOW_IN_ST = find_common_available_range_scaled_ups(
                dow_in_st_st, dow_in_st_en, startB, endB
            )
            length_DOW_IN_ST = dow_in_st_en - dow_in_st_st
            overlap_DOW_IN_ST_occurred = True

        if dow_in_en_st <= endB and dow_in_en_en >= startB:
            overlap_DOW_IN_EN, scaled_values_DOW_IN_EN = find_common_available_range_scaled_dow(
                dow_in_en_st, dow_in_en_en, startB, endB
            )
            length_DOW_IN_EN = dow_in_en_en - dow_in_en_st
            overlap_DOW_IN_EN_occurred = True

        if ex2 <= endB and dow_ex_en >= startB:
            overlap_DOW_EX, scaled_values_DOW_EX = find_common_available_range_scaled_ups(
                ex2, dow_ex_en, startB, endB
            )
            length_DOW_EX = dow_ex_en - ex2
            overlap_DOW_EX_occurred = True

        if (overlap_UPS_EX_occurred or overlap_UPS_IN_ST_occurred or overlap_UPS_IN_EN_occurred
                or overlap_AE_ST_occurred or overlap_AE_EN_occurred or overlap_DOW_IN_ST_occurred
                or overlap_DOW_IN_EN_occurred or overlap_DOW_EX_occurred):
            results.append({
                'GENE': gene, 'EVENT': event, 'STRAND': strand, 'RBP': rbp_name,

                'overlap_UPS_EX': overlap_UPS_EX, 'scaled_values_UPS_EX': scaled_values_UPS_EX, 'length_UPS_EX': length_UPS_EX,
                'overlap_UPS_IN_ST': overlap_UPS_IN_ST, 'scaled_values_UPS_IN_ST': scaled_values_UPS_IN_ST, 'length_UPS_IN_ST': length_UPS_IN_ST,
                'overlap_UPS_IN_EN': overlap_UPS_IN_EN, 'scaled_values_UPS_IN_EN': scaled_values_UPS_IN_EN, 'length_UPS_IN_EN': length_UPS_IN_EN,
                'overlap_AE_ST': overlap_AE_ST, 'scaled_values_AE_ST': scaled_values_AE_ST, 'length_AE_ST': length_AE_ST,
                'overlap_AE_EN': overlap_AE_EN, 'scaled_values_AE_EN': scaled_values_AE_EN, 'length_AE_EN': length_AE_EN,
                'overlap_DOW_IN_ST': overlap_DOW_IN_ST, 'scaled_values_DOW_IN_ST': scaled_values_DOW_IN_ST, 'length_DOW_IN_ST': length_DOW_IN_ST,
                'overlap_DOW_IN_EN': overlap_DOW_IN_EN, 'scaled_values_DOW_IN_EN': scaled_values_DOW_IN_EN, 'length_DOW_IN_EN': length_DOW_IN_EN,
                'overlap_DOW_EX': overlap_DOW_EX, 'scaled_values_DOW_EX': scaled_values_DOW_EX, 'length_DOW_EX': length_DOW_EX,
            })

    return rbp_name, index_B, results


if __name__ == "__main__":
    total_rows = len(all_data_combined)
    print(f"[{time.strftime('%H:%M:%S')}] Starting processing of {total_rows} rows...", flush=True)
    print(f"  events file: {args.events}", flush=True)
    print(f"  peaks file:  {args.peaks}", flush=True)
    print(f"  output file: {args.output}", flush=True)

    # NEW: only pass (index, row) - not the whole 500MB dataframe - per task
    args_list = [(index, row) for index, row in all_data_combined.iterrows()]

    flat_results = []
    completed = 0
    PRINT_EVERY = 2000  # NEW: printing+flushing on every single row (of 9.25M) adds real overhead

    n_workers = max(1, os.cpu_count() - 2)  # NEW: os.cpu_count()-10 can be 0/negative on smaller machines
    print(f"Using {n_workers} worker processes", flush=True)

    with Pool(n_workers) as pool:
        # NEW: chunksize batches tasks instead of dispatching 9.25M individually
        for rbp_name, index_B, row_results in pool.imap_unordered(process_row, args_list, chunksize=500):
            completed += 1
            flat_results.extend(row_results)

            if completed % PRINT_EVERY == 0 or completed == total_rows:
                elapsed = time.time() - starttime
                avg_per_row = elapsed / completed
                remaining = avg_per_row * (total_rows - completed)
                pct = completed / total_rows * 100
                print(
                    f"[{time.strftime('%H:%M:%S')}] Done {completed}/{total_rows} ({pct:.2f}%) "
                    f"| last RBP: {rbp_name} (row {index_B}) "
                    f"| elapsed: {elapsed/60:.1f} min | est. remaining: {remaining/60:.1f} min",
                    flush=True
                )

    EclipSE_Table = pd.DataFrame(flat_results)
    EclipSE_Table.to_csv(args.output, index=False)

    print('That took {} seconds'.format(time.time() - starttime))
