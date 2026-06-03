# Traning Custom Silva Classifier 
wget https://data.qiime2.org/2024.2/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2024.2/common/silva-138-99-tax.qza
qiime feature-classifier extract-reads –i-sequences silva-138-99-seqs.qza –p-f-primer CCTACGGGNGGCWGCAG –p-r-primer GACTACHVGGGTATCTAATCC –p-trunc-len 0 –o-reads silva-138-99-v3v4-seqs.qza
qiime feature-classifier fit-classifier -naïve-bayes –i-reference-reads silva-138-99-v3v4-seqs.qza –i-reference-taxonomy silva-138-99-tax.qza –o-classifier silva-138-99-v3v4-classifier.qza

# Merge datasets:
qiime feature-table merge –i-tables feature-table.qza feature-table-3. feature-table-4.qza feature-table-5.qza –o-merged-table merged-table.qza
qiime feature-table merge-seqs –i-data asv-sequences.qza asv-sequences-3.qza asv-sequences-4.qza asv-sequences-5.qza –o-merged-data merged-rep-seqs.qza
qiime feature-classifier classify-consensus-vsearch   --i-query merged-rep-seqs.qza   --i-reference-reads silva-138-99-v3v4-seqs.qza   --i-reference-taxonomy silva-138-99-tax.qza   --o-search-results vsearch-search-results.qza   --o-classification taxonomy.qza   --p-perc-identity 0.8 --p-threads 6
qiime metadata tabulate –m-input-file taxonomy.qza –o-visualization taxonomy.qzv
qiime tools export –input-path taxonomy.qzv –output-path taxonomy
qiime taxa filter-table –i-table merged-table.qza –i-taxonomy taxonomy.qza –p-mode-contains –p-include p__ --p-exclude ‘p__;Chloroplast,Mitochondria’ –o-filtered-table filtered-table.qza
qiime feature-table filter-seqs –i-data merged-rep-seqs.qza –i-table filtered-table.qza –o-filtered-data filtered-rep-seqs.qza
qiime taxa barplot –i-table filtered-table.qza –i-taxonomy taxonomy.qza –m-metadata-file metadata_merged_NSCLC.txt –o-visualization taxa-bar-plots.qzv
qiime tools export –input-path taxa-bar-plots.qzv –output-path taxa
