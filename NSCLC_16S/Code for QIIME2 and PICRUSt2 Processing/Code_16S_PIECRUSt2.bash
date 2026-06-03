#!/bin/bash

# Activate environment and install/load picrust2
source ~/miniconda3/bin/activate
conda create -n picrust2 -c conda-forge -c bioconda picrust2=2.5.1 -y
conda activate picrust2
qiime tools export --input-path filtered-rep-seqs.qza --output-path exported-seqs
qiime tools export --input-path filtered-table.qza --output-path exported-table
picrust2_pipeline.py \
  -s exported-seqs/dna-sequences.fasta \
  -i exported-table/feature-table.biom \
  -o picrust2_out_pipeline \
  -p 4 \
  --stratified
 
# Export objects for RStudio Analysis
biom convert -i exported-table/feature-table.biom -o exported-table/feature-table.tsv --to-tsv
qiime tools export --input-path taxonomy.qza --output-path exported-taxonomy
 
