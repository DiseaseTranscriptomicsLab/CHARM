import argparse
import pandas as pd
import numpy as np
import time, os
from multiprocessing import Pool

starttime = time.time()

# -----------------------------------------------------------------------------
# CLI arguments
#
# Originally this script hardcoded its two input filenames -
#   ../Final_IR_Table_hg19.txt (splicing-genome / intron-retention event table)
#   ./RBM39.txt                (eCLIP peaks table, named after one specific
#                                example RBP rather than being a generic name)
# - and its output filename, ./Intron_events_hg19_RBM39.txt. These are now
# CLI flags (defaults match the original hardcoded names, so old invocations
# without flags still work). See the ECLIPSE README for exactly which file
# to pass. NOTE: there is currently no rMATS intron-retention splicing-genome
# builder in this folder, so this script's --events input should be the
# output of Eclip_preprocessing_vttools.R (VastDB/VAST-TOOLS path only).
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Intersect eCLIP peaks against intron-retention splicing-genome regions "
                "(the 'binding metagenome assembly' step of the eCLIPSE pipeline)."
)
parser.add_argument("--events", default="../Final_IR_Table_hg19.txt",
                     help="Intron-retention splicing-genome / event-region table produced "
                          "by Eclip_preprocessing_vttools.R. Default: ../Final_IR_Table_hg19.txt")
parser.add_argument("--peaks", default="./RBM39.txt",
                     help="eCLIP peaks table, space-delimited with columns "
                          "Chrom/StartCord/EndCord/RBP (e.g. all_data_combined_38.txt). "
                          "Default: ./RBM39.txt")
parser.add_argument("--output", default="./Intron_events_hg19_RBM39.txt",
                     help="Output path for the per-RBP, per-event overlap table. "
                          "Default: ./Intron_events_hg19_RBM39.txt")
args = parser.parse_args()

#Establishing the file that will be used
map = pd.read_csv(args.events)

map = map[map['GENE'].astype(str) != ""]
NameGiver = map.dropna()

all_data_combined = pd.read_csv(args.peaks, sep=" ")


# Define a scaling function
def find_common_available_range_scaled_ups(range1_start, range1_end, range2_start, range2_end):
    # Determine the common available range
    common_start = max(range1_start, range2_start)
    common_end = min(range1_end, range2_end)

    # Check if there is an overlap
    if common_start <= common_end:
        # Calculate scaled values
        scale_factor = range1_start # Adjust the scale factor to start from day 1
        common_start_scaled = common_start - scale_factor + 1
        common_end_scaled = common_end - scale_factor

        # Return actual and scaled values
        return (common_start, common_end), (common_start_scaled, common_end_scaled)
    else:
        return None, None  # No common available range

def find_common_available_range_scaled_dow(range1_start, range1_end, range2_start, range2_end):
    # Determine the common available range
    common_start = max(range1_start, range2_start)
    common_end = min(range1_end, range2_end)

    # Check if there is an overlap
    if common_start <= common_end:
        # Calculate scaled values
        scale_factor = range1_end # Adjust the scale factor to start from day 1
        common_start_scaled = scale_factor - common_start
        common_end_scaled = scale_factor - common_end + 1

        # Return actual and scaled values
        return (common_start, common_end), (common_start_scaled, common_end_scaled)
    else:
        return None, None  # No common available range


# --- Ported from eCLIPSE_Exon.py: precompute per-chromosome NumPy arrays ONCE
# instead of filtering NameGiver for every single one of the ~9M eCLIP rows.
# On Linux, Pool uses fork by default, so this is built once in the parent
# and inherited (copy-on-write) by workers for free.

META_COLS = ['GENE', 'EVENT', 'STRAND']
REGION_COLS = [
    'UPS.ST', 'UPS.EN',
    'COORD.ST', 'IN.ST',
    'IN.EN', 'COORD.EN',
    'DOW.ST', 'DOW.EN'
]

NameGiver_by_chrom = {}
for chrom, df in NameGiver.groupby('CHROM'):
    NameGiver_by_chrom[chrom] = {
        'meta': df[META_COLS].to_numpy(dtype=object),
        'coords': df[REGION_COLS].to_numpy(dtype=float),
    }


# This is the function that calculates the binding (overlap of coordinates) of the RBPs
# to the regulatory regions defined. Specific for intron-retention coordinates.
def process_row(args):
    index_B, row_B = args  # no longer carries the full eCLIP peaks dataframe
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
    # Raw NumPy array access instead of DataFrame.iterrows() (much faster per-row)
    for i in range(n):
        (ups_st, ups_en,
         coord_st, in_st,
         in_en, coord_en,
         dow_st, dow_en) = coords[i]
        gene, event, strand = meta[i]

        overlap_UPS_EX_occurred = False
        overlap_IN_ST_occurred = False
        overlap_IN_EN_occurred = False
        overlap_DOW_EX_occurred = False

        overlap_UPS_EX = (pd.NA, pd.NA); scaled_values_UPS_EX = (pd.NA, pd.NA); length_UPS_EX = pd.NA
        overlap_IN_ST = (pd.NA, pd.NA); scaled_values_IN_ST = (pd.NA, pd.NA); length_IN_ST = pd.NA
        overlap_IN_EN = (pd.NA, pd.NA); scaled_values_IN_EN = (pd.NA, pd.NA); length_IN_EN = pd.NA
        overlap_DOW_EX = (pd.NA, pd.NA); scaled_values_DOW_EX = (pd.NA, pd.NA); length_DOW_EX = pd.NA

        # Check if eCLIP region is in the upstream exon
        if ups_st <= endB and ups_en >= startB:
            overlap_UPS_EX, scaled_values_UPS_EX = find_common_available_range_scaled_ups(
                ups_st, ups_en, startB, endB
            )
            length_UPS_EX = ups_en - ups_st
            overlap_UPS_EX_occurred = True

        # Check if eCLIP region is in the intron start
        if coord_st <= endB and in_st >= startB:
            overlap_IN_ST, scaled_values_IN_ST = find_common_available_range_scaled_ups(
                coord_st, in_st, startB, endB
            )
            length_IN_ST = in_st - coord_st
            overlap_IN_ST_occurred = True

        # Check if eCLIP region is in the intron end
        if in_en <= endB and coord_en >= startB:
            overlap_IN_EN, scaled_values_IN_EN = find_common_available_range_scaled_dow(
                in_en, coord_en, startB, endB
            )
            length_IN_EN = coord_en - in_en
            overlap_IN_EN_occurred = True

        # Check if eCLIP region is in downstream exon
        if dow_st <= endB and dow_en >= startB:
            overlap_DOW_EX, scaled_values_DOW_EX = find_common_available_range_scaled_ups(
                dow_st, dow_en, startB, endB
            )
            length_DOW_EX = dow_en - dow_st
            overlap_DOW_EX_occurred = True

        if overlap_UPS_EX_occurred or overlap_IN_ST_occurred or overlap_IN_EN_occurred or overlap_DOW_EX_occurred:
            results.append({
                'GENE': gene, 'EVENT': event, 'STRAND': strand, 'Name': rbp_name,
                'overlap_UPS_EX': overlap_UPS_EX, 'scaled_values_UPS_EX': scaled_values_UPS_EX, 'length_UPS_EX': length_UPS_EX,
                'overlap_IN_ST': overlap_IN_ST, 'scaled_values_IN_ST': scaled_values_IN_ST, 'length_IN_ST': length_IN_ST,
                'overlap_IN_EN': overlap_IN_EN, 'scaled_values_IN_EN': scaled_values_IN_EN, 'length_IN_EN': length_IN_EN,
                'overlap_DOW_EX': overlap_DOW_EX, 'scaled_values_DOW_EX': scaled_values_DOW_EX, 'length_DOW_EX': length_DOW_EX,
            })

    return rbp_name, index_B, results


if __name__ == "__main__":
    total_rows = len(all_data_combined)
    print(f"[{time.strftime('%H:%M:%S')}] Starting processing of {total_rows} rows...", flush=True)
    print(f"  events file: {args.events}", flush=True)
    print(f"  peaks file:  {args.peaks}", flush=True)
    print(f"  output file: {args.output}", flush=True)

    # Only pass (index, row) - not the whole eCLIP peaks dataframe - per task
    args_list = [(index, row) for index, row in all_data_combined.iterrows()]

    flat_results = []
    completed = 0
    PRINT_EVERY = 2000  # printing+flushing on every single row adds real overhead at this scale

    n_workers = max(1, os.cpu_count() - 2)
    print(f"Using {n_workers} worker processes", flush=True)

    with Pool(n_workers) as pool:
        # chunksize batches tasks instead of dispatching them one at a time
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
