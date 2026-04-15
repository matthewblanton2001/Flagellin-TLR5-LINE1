# Mapping to human genome
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422544.fastq.gz -S GSE164704_Untreated1.sam 
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422556.fastq.gz -S GSE164704_Untreated2.sam 
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422562.fastq.gz -S GSE164704_Untreated3.sam 
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422567.fastq.gz -S GSE164704_Untreated4.sam 
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422553.fastq.gz -S GSE164704_Flagellin_Treated1.sam 
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422559.fastq.gz -S GSE164704_Flagellin_Treated2.sam 
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422564.fastq.gz -S GSE164704_Flagellin_Treated3.sam 
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422570.fastq.gz -S GSE164704_Flagellin_Treated4.sam 

# Making Tag Directories: 
makeTagDirectory GSE164704_Untreated1 GSE164704_Untreated1.sam 
makeTagDirectory GSE164704_Untreated2 GSE164704_Untreated2.sam 
makeTagDirectory GSE164704_Untreated3 GSE164704_Untreated3.sam 
makeTagDirectroy GSE164704_Untreated4 GSE164704_Untreated4.sam 
makeTagDirectory GSE164704_Flagellin_Treated1 GSE164704_Flagellin_Treated1.sam 
makeTagDirectory GSE164704_Flagellin_Treated2 GSE164704_Flagellin_Treated2.sam 
makeTagDirectory GSE164704_Flagellin_Treated3 GSE164704_Flagellin_Treated3.sam 
makeTagDirectroy GSE164704_Flagellin_Treated4 GSE164704_Flagellin_Treated4.sam 

# Homer analysis: 
analyzeRepeats.pl repeats hg38 -d GS164704_Flagellin_Treated1 GSE164704_Flagellin_Treated2 GSE164704_Flagellin_Treated3 GSE164704_Flagellin_Treated4 GSE164704_Untreated1 GSE164704_Untreated2 GSE164704_Untreated3 GSE164704_Untreated4 -L3 L1 -noadj>outputfile.txt
