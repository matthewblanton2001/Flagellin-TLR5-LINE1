Flagellin–TLR5–LINE1 Analysis Pipeline
Overview

This repository contains all scripts used to analyze the role of bacterial flagellin in TLR5 signaling and LINE-1 activation in NSCLC.

Analyses include:

16S rRNA microbiome profiling (QIIME2)

RNA-seq analysis of LINE-1 expression (HOMER)

TCGA survival analysis (TLR5)

Repository Structure

NSCLC_16S/ – QIIME2 processing and downstream microbiome analysis

LINE1 Repeat Expression in HBEC treated with flagellin/ – RNA-seq alignment, TE repeat quantification via HOMER, and DE analysis

TCGA-LUAD-TLR5/ – Survival modeling using TLR5 from TCGA LUAD data

metadata/ – Sample metadata for microbiome analysis

Requirements
Software
QIIME2 (2024.10)
R (≥4.4.1)
HISAT2
HOMER
Reproducibility

Code Availability

All scripts are available in this repository.

Notes: If any part of the code needs further clarification please make an Issue and I will reach back as soon as possible.

Also, paths from manifest files may need to be adjusted depending on local file structure.
