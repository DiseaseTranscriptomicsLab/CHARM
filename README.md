# CHARM <a href="https://diseasetranscriptomicslab.github.io/CHARM/"><img src="Charm_logo.png" align="right" height="139"/></a>

[![DOI](https://zenodo.org/badge/1040091867.svg)](https://doi.org/10.5281/zenodo.21724194)

**Comprehensive Hub for Alternative Regulatory Mapping (CHARM)**  

CHARM is an R Shiny web application for exploring the regulatory roles of RNA-binding proteins (RBPs) in **gene expression, splicing regulation, and direct RNA binding**.  

The app integrates data from the **ENCODE project** (Luo et al., 2020), including eCLIP and RNA-seq datasets from RBP knockdown/knockout experiments in **HEPG2 (liver cancer)** and **K562 (leukaemia)** cell lines.  

CHARM serves as both a **repository** of this curated data and a **discovery tool**:  
- Explore how 168 RBPs from ENCODE affect expression, splicing, and binding.  
- Infer potential networks between these regulatory layers, revealing new pathways and mechanisms.  
- Upload your own expression, splicing, or binding datasets to identify RBPs most likely altered in your system.  

---

## Data Layers in CHARM  

CHARM integrates three main data types:  

- **Expression**: Differential expression upon RBP knockdown/knockout, with pathway-level insights from Gene Set Enrichment Analysis.  

- **Splicing**: Alternative splicing changes quantified with **betAS** (Ferreira et al., 2024), using **VastDB** nomenclature (Tapial et al., 2017).  

- **Binding**: Altered RNA binding patterns of the silenced RBP and other RBPs, characterized with the **eCLIPSE** tool (see below).  

---

## eCLIPSE (eCLIP for Splicing Evaluation)  

**eCLIPSE** is a companion tool that processes eCLIP data from ENCODE to generate **RNA splicing maps**. It aligns binding profiles to the genomic coordinates of splicing events defined in **VastDB** (Tapial et al., 2017).  

This allows assessment of how each RBP regulates splicing events both **directly** (through its own binding) and **indirectly** (through effects on other RBPs).  

eCLIPSE was used to generate the binding layer in CHARM. The tool is freely available and can also be applied to user-submitted eCLIP datasets to map splicing regulation of additional RBPs.  

---

## Documentation

Full tutorials are published at **[diseasetranscriptomicslab.github.io/CHARM](https://diseasetranscriptomicslab.github.io/CHARM/)**:

- **CHARM App Tutorial** — using the Shiny app itself, no code required.
- **eCLIPSE Local Pipeline** — running eCLIPSE outside the app: building your own splicing genome (VAST-TOOLS/VastDB or rMATS), intersecting it with your own eCLIP peaks, and generating RNA binding maps.

Source for both lives in [`docs/`](docs/).

---

## Repository Structure

| Path | What it is |
|---|---|
| `app.R` | The CHARM Shiny app itself (UI + server). |
| `helper_functions.R` | Shared R functions used by `app.R` — data loading, statistics, and plotting (RNA binding maps, heatmaps, network views, etc.). |
| `Dockerfile`, `.dockerignore` | Container build for deploying the app. |
| `www/` | Static assets served by the Shiny app (logo). |
| `Images/` | Screenshots and figures used in documentation/READMEs. |
| `data/` *(not on GitHub)* | Precomputed app data (ENCODE-derived expression/splicing/binding tables, `.qs2` objects). Too large for git — see [Large data files](#large-data-files) below. |
| `example_data/` | Small example input files for trying the app's Discovery-mode upload features. |
| `eCLIPSE/` | The eCLIPSE pipeline: preprocessing scripts (`Eclip_preprocessing_*.R`, `eCLIPSE_Exon.py`, `eCLIPSE_Intron.py`), the RNA-binding-map notebooks (`Eclip_position_matrix.Rmd`, `eCLIPSE_plotting_functions.Rmd`), and small reference tables. Some raw event tables here are too large for git — see [Large data files](#large-data-files). See [`eCLIPSE/README.md`](eCLIPSE/README.md) for details on running the pipeline. |
| `Markdowns/` | Supplementary R Markdown notebooks documenting each app module (Expression, Splicing, Binding) and the original eCLIPSE analysis notebook these pipeline scripts were extracted from. |
| `docs/` | The [GitHub Pages](https://diseasetranscriptomicslab.github.io/CHARM/) documentation site (app tutorial + eCLIPSE local-pipeline tutorial) and downloadable user-facing docs (`.docx`). |
| `LocalJob_File*.R`, `app_binding_similar_patches.R`, `*_all_vs_all.py` | One-off/personal analysis and patch scripts used during development; not part of the app or pipeline proper. |
| `LICENSE`, `DESCRIPTION`, `NAMESPACE`, `CHARM.Rproj` | R package/project metadata. |

### Large data files

A few files are intentionally excluded from git (see `.gitignore`) because they exceed or approach GitHub's 100 MB per-file limit: the `data/` folder (precomputed app data, several files close to 1 GB) and four raw eCLIPSE event tables (`eCLIPSE/all_data_combined_38.txt`, `eCLIPSE/BigExon_events.txt`, `eCLIPSE/Exon_events_Rmats.txt`, `eCLIPSE/Intron_events.txt`). These are distributed via Zenodo instead — see the [Zenodo record](https://doi.org/10.5281/zenodo.21724194) for download links, or regenerate them from ENCODE/VastDB following [`eCLIPSE/README.md`](eCLIPSE/README.md).
