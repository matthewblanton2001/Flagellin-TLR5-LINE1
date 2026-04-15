#!/bin/bash

# Activate environment
source ~/miniconda3/bin/activate
conda activate qiime2-amplicon-2024.10

# Create Directory
mkdir NSCLC5
cd NSCLC5

# Import data
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest_5.txt \
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

# Summarize
qiime demux summarize \
  --i-data single-end-demux.qza \
  --o-visualization demux.qzv
qiime tools export --input-path demux.qzv --output-path Intratumoral

# DADA2 processing
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs single-end-demux.qza \
  --p-trim-left 17 \
  --p-trunc-len 250 \
  --o-representative-sequences asv-sequences-5.qza \
  --o-table feature-table-5.qza \
  --o-denoising-stats dada2-stats.qza
qiime metadata tabulate --m-input-file dada2-stats.qza --o-visualization dada2-stats-summ.qzv
qiime tools export --input-path dada2-stats-summ.qzv --output-path DADA2