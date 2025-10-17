# CHARM

**Comprehensive Hub for Alternative Regulatory Mapping (CHARM)**  

CHARM is an R Shiny web application for exploring the regulatory roles of RNA-binding proteins (RBPs) in **gene expression, splicing regulation, and direct RNA binding**.  

The app integrates data from the **ENCODE project** (Luo et al., 2020), including eCLIP and RNA-seq datasets from RBP knockdown/knockout experiments in **HEPG2 (liver cancer)** and **K562 (leukemia)** cell lines.  

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
