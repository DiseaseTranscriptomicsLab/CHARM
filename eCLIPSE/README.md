# Running eCLIPSE Locally: Building Your Own Splicing Genome, Binding Metagenome, and RNA Binding Maps

This is a tutorial for computational users who want to run CHARM's underlying eCLIPSE pipeline
directly - outside the Shiny app - to build a splicing genome from their own alternative-splicing
data (VAST-TOOLS/VastDB **or** rMATS), intersect it against their own eCLIP peak files, and
generate their own RNA binding maps (chi-squared statistics + the increased/decreased/maintained
density plots used throughout the CHARM paper).

This audience and workflow is different from the CHARM Shiny app itself: there is no upload
button here, and several of the scripts in this folder are working research-notebook scripts
rather than a packaged tool. This README explains, step by step, what each script does, what
file format it expects, and - importantly - **what was broken in the original notebook and what
was changed to fix it**, so you know exactly what you're running.

---

## 1. Overview: the three pipeline stages

Regardless of whether your splicing quantification comes from VAST-TOOLS/VastDB or rMATS, the
pipeline has the same three stages:

1. **Splicing-genome construction** - take a table of alternative-splicing events (exon-skipping
   or intron-retention) and, for each event, define a fixed set of genomic sub-regions around it
   (flanking exon ends, flanking intron ends, the alternative exon itself, etc. - see the CHARM
   paper's Methods, "Splicing-genome creation"). This produces one row per event with genomic
   coordinates for each sub-region.
2. **Binding metagenome assembly** - intersect a table of eCLIP peaks (RBP, chromosome, start,
   end) against every region of every event, recording where and how much of each region a given
   RBP's peaks cover. This produces a long-format, per-RBP-per-event overlap table.
3. **RNA binding maps** - reshape the overlap table into a wide "position matrix" (metagenomic
   coordinates 1-1000 for exon-skipping, 1-500 for intron-retention), split events by direction of
   splicing change (increased/decreased/maintained inclusion), and run a chi-squared test at every
   position to test whether the RBP's binding density differs between event classes. This is what
   `eCLIPSE_full()` / `eCLIPSE_full_IR()` do, and is exactly what the CHARM paper's Figures show.

Two independent event sources feed stage 1: **VastDB/VAST-TOOLS** (event IDs like
`HsaEX0024923`, `HsaINT...`) and **rMATS** (`SE.MATS.JC.txt`, numeric event IDs). Both converge
on the same stage-2/stage-3 scripts - that convergence is the whole point of rMATS compatibility,
and is what makes this folder useful if you don't use VAST-TOOLS.

```
 VAST-TOOLS align/combine  ──▶ Eclip_preprocessing_vttools.R ──▶ Final_MIC_Table.txt (exons)
                                                               └▶ Final_IR_Table_hg19.txt (introns)

 rMATS (SE.MATS.JC.txt)    ──▶ Eclip_preprocessing_rmats.R   ──▶ Final_AS_Table_Rmats.txt (exons)

                     (either exon table) ──▶ eCLIPSE_Exon.py   ──▶ overlap table (long format)
                     (VastDB intron table only) ──▶ eCLIPSE_Intron.py ──▶ overlap table (long format)

                     overlap table ──▶ [reshape to position matrix, §5] ──▶ eCLIPSE_full() / eCLIPSE_full_IR()
                                                                             (eCLIPSE_plotting_functions.R)
                                                                             ──▶ RNA binding map plot
```

---

## 2. Environment / dependencies

**R** (tested against R packages used throughout the notebook):
`data.table`, `tidyr`, `dplyr`, `stringr`, `ggplot2`, `ggpubr`, `scales`

**Python 3** (used by `eCLIPSE_Exon.py` / `eCLIPSE_Intron.py`):
`pandas`, `numpy`

**External tools**, only needed for the VAST-TOOLS/VastDB path:
[`vast-tools`](https://github.com/vastgroup/vast-tools) (align/combine), and `hg38`/`hg19`
reference indices as required by vast-tools.

**External tools**, only needed for the rMATS path:
[`rMATS-turbo`](https://github.com/Xinglab/rmats-turbo), run beforehand to produce
`SE.MATS.JC.txt`.

> **Important:** I was not able to run or test any R code in the environment used to write this
> tutorial (no R installation was available), so none of the R scripts below - including the
> bug fix in `eCLIPSE_plotting_functions.R` - have been executed end-to-end. The extraction was
> done mechanically (exact line ranges copied from `eCLIPSE_Analysis.Rmd`, not retyped), and the
> fix is a small, easily-reviewed change (see §6), but please test on a small example before
> trusting results for analysis.

---

## 3. Path A: VAST-TOOLS / VastDB splicing genome

### 3.1 Align and combine (`vt_align.sh`, `vt_combine.sh`)

`vt_align.sh` is a thin wrapper around `vast-tools align`. Run it from the directory containing
your paired-end FASTQs (named `<sample>_1.paired.fastq.gz` / `<sample>_2.paired.fastq.gz`):

```bash
./vt_align.sh /path/to/fastqs/
```

Then combine all aligned samples into one inclusion table:

```bash
./vt_combine.sh
```

This produces a `vast-tools combine` output (an `INCLUSION_LEVELS_FULL-hg38-*.tab.gz` file) -
this is the `vttable` input to the next script. **You will need to edit `vt_combine.sh` if your
sample count/design differs from the hardcoded example** (it currently names its log file after
one specific experiment, `VT_align_shRNA_RBFOX2.txt`, but the actual `vast-tools combine` command
itself is generic).

### 3.2 Build the splicing genome (`Eclip_preprocessing_vttools.R`)

This script reads a VAST-TOOLS inclusion table and eCLIP peak files, and writes the
exon-skipping and intron-retention splicing-genome tables.

**Before running it, open the script and edit these hardcoded paths** (this script was written
for one person's specific lab filesystem, so nothing here will resolve on your machine as-is):

| Line | What it reads/writes | You need to point it at |
|---|---|---|
| 10, 12 | `all_data_combined_38.txt` (hg38 combined eCLIP peaks) | your combined eCLIP peaks file, or the downloadable example data (§7) |
| 14, 31 | `all_combined_data.txt` (hg19 combined eCLIP peaks) | same, hg19 version - **only needed if you also want the legacy hg19 path**; most users can ignore lines 14-31 entirely |
| 33 | `INCLUSION_LEVELS_FULL-hg38-25-v251.tab.gz` | the `vast-tools combine` output from §3.1 |
| 94 | `vttable_IR_regions_hg19.RDS` | output path for the intermediate intron-retention regions object |
| 209 | `vttable_MIC_exons.RDS` | output path for the intermediate exon-skipping regions object |
| 315 | `Final_IR_Table_hg19.txt` | final intron-retention splicing-genome table (feeds `eCLIPSE_Intron.py`) |
| 316 | `Final_MIC_Table.txt` | final exon-skipping splicing-genome table (feeds `eCLIPSE_Exon.py`) |

Run with `Rscript Eclip_preprocessing_vttools.R` (after editing paths) or interactively.

> **Known issue - filename mismatch:** this script writes the exon-skipping table as
> `Final_MIC_Table.txt` (line 316), but `eCLIPSE_Exon.py` hardcodes its input default as
> `Final_BIGEX_Table.txt`, and `Eclip_preprocessing_rmats.R` (the rMATS equivalent, §4) writes
> `Final_AS_Table_Rmats.txt`. **None of these three names match each other.** This has been
> fixed on the `eCLIPSE_Exon.py` side (see §6) by making the input filename a `--events` CLI flag
> - so you no longer need to rename files by hand, just pass the exact filename each script
> actually wrote. The R scripts' internal output filenames were left as-is (not renamed) to avoid
> touching logic in scripts I could not execute/test; only the *documentation* + the
> already-safe, mechanical *Python-side* CLI parameterization were changed.

> **Known issue - stale "MIC pipeline" comment:** line 6 of this script reads
> `# ON THE 20/01/2025, you did ALTERATIONS FOR MIC PIPELINE. DONT FORGET TO RECHANGE THEM`. This
> is a note-to-self left in the script; from reading the code, the final exon-skipping output
> (`vttable_EX_bigexon_final`, written to `Final_MIC_Table.txt`) is built from `vttable_EX`
> filtered to `grepl("^HsaEX", EVENT)` (line 44) - i.e. **all** exon-skipping events, not only
> "MIC" (microexon) events - so as far as I could trace by reading the code, the current script
> state does **not** appear to be MIC-restricted. I flag the comment here rather than silently
> resolving it, since I could not execute the script to confirm this by testing. If you know what
> the "MIC pipeline alterations" refer to, it's worth double-checking before you rely on the
> output for a new analysis.

---

## 4. Path B: rMATS splicing genome

### 4.1 Run rMATS-turbo yourself

Run rMATS-turbo per its own documentation to produce `SE.MATS.JC.txt` (exon-skipping events).
There is currently **no intron-retention equivalent** on the rMATS path in this folder - only
exon-skipping is supported for rMATS today.

### 4.2 Build the splicing genome (`Eclip_preprocessing_rmats.R`)

This script converts `SE.MATS.JC.txt` into the same 20-column splicing-genome schema that
`Eclip_preprocessing_vttools.R` produces, so it can be fed straight into `eCLIPSE_Exon.py`
alongside (or instead of) the VastDB path.

Edit these hardcoded paths before running:

| Line | What it reads/writes | You need to point it at |
|---|---|---|
| 6 | `SE.MATS.JC.txt` | your rMATS-turbo output |
| 7 | `Final_BIGEX_Table.txt` | **an already-existing VastDB splicing-genome table** - see note below |
| 93, 97 | `rmatstable_as.RDS` | intermediate output/reload path |
| 104 | `Final_AS_Table_Rmats.txt` | final rMATS-derived splicing-genome table (feeds `eCLIPSE_Exon.py`) |

> **Important dependency you might miss:** line 7 (`vttable <- fread(".../Final_BIGEX_Table.txt")`)
> loads an existing VastDB-derived table purely to copy its column names (`newnames <-
> colnames(vttable)`, line 101) so the rMATS table ends up with matching column names. **You need
> to have already run the VAST-TOOLS path (§3) at least once**, or otherwise obtain a table with
> the correct 20 column names, before this script's renaming step (line 101-102) will work. This
> is not mentioned anywhere in the original script.

> **Data note:** this script converts rMATS's 0-based coordinates to 1-based (lines 16-21) before
> doing anything else - this is correct and necessary, since VAST-TOOLS/VastDB coordinates
> (which the rest of the pipeline assumes) are 1-based.

> **Minor known issue (not changed):** in the downstream-intron region-splitting block (around
> line 66-81), the minus-strand branch (`else` at line 73) computes the same split as the
> plus-strand branch, whereas the equivalent upstream-intron block (lines 50-64) correctly swaps
> the upstream/downstream half-lengths for minus-strand events. This is a small, easily-missed
> asymmetry (it affects, by at most a base or two due to floor/ceiling rounding, how a
> **downstream intron shorter than 400bp** is split on the minus strand) that I noticed while
> reading the script but have **not** changed, since I could not run the script to confirm the
> effect or regression-test a fix against your existing results. Flagging it here so you can
> decide whether it matters for your use case.

---

## 5. Binding metagenome assembly (`eCLIPSE_Exon.py`, `eCLIPSE_Intron.py`)

These scripts intersect your eCLIP peaks against the splicing-genome regions from §3/§4, per RBP
per event, and write a long-format overlap table.

Both scripts now take CLI flags instead of hardcoded filenames (this was the second filename
mismatch - see §6):

```bash
# Exon-skipping (works for both the VastDB and rMATS splicing-genome tables)
python3 eCLIPSE_Exon.py \
  --events Final_MIC_Table.txt \
  --peaks  all_data_combined_38.txt \
  --output BenFile_ICLIP.txt

# Intron retention (VastDB path only - no rMATS IR table exists yet)
python3 eCLIPSE_Intron.py \
  --events Final_IR_Table_hg19.txt \
  --peaks  all_data_combined_38.txt \
  --output Intron_events_hg19.txt
```

Omitting the flags falls back to the original hardcoded defaults, so any existing invocations
you had will still behave the same way.

`--peaks` expects a **space-delimited** table (as produced by R's `write.table(..., row.names =
F)`) with at minimum the columns `Chrom`, `StartCord`, `EndCord`, `RBP` (the combined eCLIP peaks
file in this folder, `all_data_combined_38.txt`, already has exactly this format - see §7).

**Performance note:** `eCLIPSE_Exon.py` is written for large-scale runs (multiprocessing,
per-chromosome NumPy pre-indexing - it can comfortably chew through ~9M eCLIP peak rows).
`eCLIPSE_Intron.py` is an older, unoptimized version (`DataFrame.iterrows()`-based, no
per-chromosome pre-indexing) - it will be substantially slower on large inputs. If you plan to
run the intron-retention path at scale, budget for that, or consider porting the same
optimizations from `eCLIPSE_Exon.py`.

---

## 6. RNA binding maps (`eCLIPSE_plotting_functions.R`)

`eCLIPSE_plotting_functions.R` (new in this tutorial) is a clean, sourceable extraction of the
`get_pos_df()`, `eCLIPSE_full()`, and `eCLIPSE_full_IR()` functions from
`eCLIPSE_Analysis.Rmd`'s "Function Based System" section. It was extracted mechanically (exact
`sed` line ranges, not retyped) to avoid transcription errors in ~400 lines of statistical/plotting
code.

```r
source("eCLIPSE_plotting_functions.R")
```

### 6.1 Reshaping `eCLIPSE_Exon.py` output into a position matrix

`eCLIPSE_full()`/`eCLIPSE_full_IR()` do **not** take the long-format overlap table from
`eCLIPSE_Exon.py`/`eCLIPSE_Intron.py` directly - they expect a wide `rnamapfile`: one row per
(event, RBP), one column per metagenomic position (1-1000 for exon-skipping, 1-500 for intron
retention), plus `EVENTS` and `RBP` columns. This reshaping step exists in
`eCLIPSE_Analysis.Rmd`'s "## Visualization" chunk but is not currently a standalone script - you
will need to adapt that chunk (or write your own) to go from the long overlap table to this wide
matrix, using `get_pos_df()` to define which output columns from the overlap table map to which
of the 1..1000 (or 1..500) metagenome positions. This is the one stage of the pipeline that is
still "notebook code" rather than a runnable script - flagging it clearly rather than guessing at
a rewrite I can't test.

### 6.2 Calling the plotting function

```r
plot <- eCLIPSE_full(
  rnamapfile = my_position_matrix,   # wide matrix, see 6.1
  ASfile     = my_splicing_table,    # needs columns Event.ID, dPSI
  rnaBP      = "RBM39",              # RBP name to plot, must match a value in rnamapfile$RBP
  metric     = "FDR",                # or "EffectSize"
  title      = "My experiment"
)
print(plot)
```

`ASfile` needs a `Event.ID` column (VastDB `HsaEX...`/`HsaINT...` IDs **or** rMATS numeric IDs -
both now work, see below) and a `dPSI` column (signed delta-PSI, used to classify events as
increased/decreased/maintained inclusion against `PSIthreshold`, default 0.05).

Use `eCLIPSE_full()` for exon-skipping data (VastDB or rMATS) and `eCLIPSE_full_IR()` for
intron-retention data (VastDB only).

> **Note:** `eCLIPSE_Analysis.Rmd`'s "## Running" example chunk calls a function named
> `eCLIPSE()`. That function does not exist anywhere in this codebase - it is simply an old name
> for `eCLIPSE_full()`. Always call `eCLIPSE_full()`.

### 6.3 Bug fix applied: event-ID filter was VastDB-only

The single functional bug found and fixed in this tutorial: both `eCLIPSE_full()` and
`eCLIPSE_full_IR()` started with a hardcoded filter,

```r
# eCLIPSE_full():
ASfile <- ASfile[grepl("HsaEX", Event.ID)]
# eCLIPSE_full_IR():
ASfile <- ASfile[grepl("HsaINT", Event.ID)]
```

`"HsaEX..."`/`"HsaINT..."` are VastDB/VAST-TOOLS-specific event-ID formats. rMATS event IDs are
plain integers (`1`, `2`, `3`, ...) - so for rMATS-derived `ASfile` tables, this `grepl()` matched
**zero rows**, silently emptying `ASfile` before any of the actual binding-map logic ran.
Critically, this does not raise an error - `eCLIPSE_full()` would proceed and simply report every
RBP as not having enough events (`"<RBP> does not have enough events."`), which looks like a data
problem rather than a code bug.

`eCLIPSE_plotting_functions.R` replaces both lines with:

```r
ASfile <- ASfile[ASfile$Event.ID %in% rnamapfile$EVENTS, ]
```

This keeps the original intent (drop `ASfile` rows unrelated to the event type `rnamapfile` was
built for) but works regardless of whether `Event.ID` uses VastDB or rMATS naming, since it
filters by "does this event actually have binding-metagenome data" rather than by ID string
pattern. This is the fix that makes the rMATS path usable end-to-end with `eCLIPSE_full()`.

---

## 7. The large data files in this folder

Several files here are large precomputed datasets rather than code:

| File | Size | Contents |
|---|---|---|
| `all_data_combined_38.txt` | ~450 MB | Combined ENCODE eCLIP peaks (hg38): `Chrom`, `StartCord`, `EndCord`, `Strand`, `Signal`, `PValue`, `RBP`, `Cell`, `Rep` |
| `BigExon_events.txt` | ~88 MB | Precomputed VastDB exon-skipping binding overlap table (`eCLIPSE_Exon.py` output) |
| `Intron_events.txt` | ~309 MB | Precomputed VastDB intron-retention binding overlap table (`eCLIPSE_Intron.py` output) |
| `Exon_events_Rmats.txt` | ~278 MB | Precomputed rMATS exon-skipping binding overlap table (numeric event IDs) |

These are included as **downloadable example data** so you can try the pipeline (or just stages
5-6) without first needing your own eCLIP peaks or splicing quantification - matching how the
CHARM paper already cites this data via its Zenodo DOI. If you're setting this repository up for
others to clone, these files are large enough that they likely belong hosted alongside the paper's
existing Zenodo deposit rather than committed to git directly.

---

## 8. Quick end-to-end example (rMATS path, using the example data)

```bash
# 1. Splicing genome (needs a VastDB Final_MIC_Table.txt already built once, for column names - §4.2)
Rscript Eclip_preprocessing_rmats.R

# 2. Binding metagenome assembly, using the downloadable example eCLIP peaks
python3 eCLIPSE_Exon.py \
  --events Final_AS_Table_Rmats.txt \
  --peaks  all_data_combined_38.txt \
  --output my_rmats_overlap_table.txt
```

```r
# 3. Reshape my_rmats_overlap_table.txt into a position matrix (§6.1 - adapt from
#    eCLIPSE_Analysis.Rmd's "## Visualization" chunk), then:
source("eCLIPSE_plotting_functions.R")
plot <- eCLIPSE_full(
  rnamapfile = my_position_matrix,
  ASfile     = my_rmats_splicing_table,   # from SE.MATS.JC.txt, with Event.ID + dPSI columns
  rnaBP      = "RBM39",
  metric     = "FDR"
)
print(plot)
```

---

## 9. Summary of changes made in this tutorial

| File | Change | Risk |
|---|---|---|
| `eCLIPSE_plotting_functions.R` (new) | Mechanical extraction of `get_pos_df`, `eCLIPSE_full`, `eCLIPSE_full_IR` from `eCLIPSE_Analysis.Rmd`; fixed the VastDB-only `grepl("HsaEX"/"HsaINT", Event.ID)` filter to be event-ID-format-agnostic (§6.3) | Extraction verified by exact line-range copy + brace/paren balance check; the fix itself is a 1-line, easily-reviewed change per function, but **not executed/tested** (no R available in the environment this was written in) |
| `eCLIPSE_Exon.py` | Added `--events`/`--peaks`/`--output` CLI flags (defaults unchanged, so existing calls still work) | Low - purely mechanical, logic untouched; confirmed to compile |
| `eCLIPSE_Intron.py` | Same CLI flags as above; fixed `Pool(os.cpu_count()-10)`, which crashes on any machine with ≤10 CPU cores, to `Pool(max(1, os.cpu_count()-2))` (matching `eCLIPSE_Exon.py`'s existing pattern) | Low - mechanical + a defensive one-line fix; confirmed to compile |
| `Eclip_preprocessing_vttools.R`, `Eclip_preprocessing_rmats.R` | **Not modified.** Hardcoded personal paths and the "MIC pipeline" comment are documented above (§3.2, §4.2) rather than edited, since these scripts run complex row-by-row biological coordinate logic I could not execute or regression-test | N/A - left as-is by design |

If anything above doesn't match what you observe when you actually run the pipeline, that's the
best signal that one of my readings of the (untested) code was wrong - please flag it back.
