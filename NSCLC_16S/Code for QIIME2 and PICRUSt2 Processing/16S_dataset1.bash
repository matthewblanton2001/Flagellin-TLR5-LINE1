#!/bin/bash

# Activate environment
source ~/miniconda3/bin/activate
conda activate qiime2-amplicon-2024.10

# Create Directory
mkdir NSCLC
cd NSCLC

# Import data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_1.txt \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# Summarize
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv
qiime tools export --input-path demux.qzv --output-path Intratumoral

# DADA2 processing
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 17 \
  --p-trim-left-r 21 \
  --p-trunc-len-f 280 \
  --p-trunc-len-r 220 \
  --o-representative-sequences asv-sequences.qza \
  --o-table feature-table.qza \
  --o-denoising-stats dada2-stats.qza
qiime metadata tabulate --m-input-file dada2-stats.qza --o-visualization dada2-stats-summ.qzv
qiime tools export --input-path dada2-stats-summ.qzv --output-path DADA2