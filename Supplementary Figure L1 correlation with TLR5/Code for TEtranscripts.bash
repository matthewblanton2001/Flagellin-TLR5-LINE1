# Map to the human genome
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422544.fastq.gz -S GSE164704_Untreated1.sam
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422556.fastq.gz -S GSE164704_Untreated2.sam
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422562.fastq.gz -S GSE164704_Untreated3.sam
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422567.fastq.gz -S GSE164704_Untreated4.sam
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422553.fastq.gz -S GSE164704_Flagellin_Treated1.sam
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422559.fastq.gz -S GSE164704_Flagellin_Treated2.sam
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422564.fastq.gz -S GSE164704_Flagellin_Treated3.sam
hisat2 -p 8 –dta Index/hg38/genome -U SRR13422570.fastq.gz -S GSE164704_Flagellin_Treated4.sam
	
# Convert SAM to BAM
samtools view -@ 8 -bS GSE164704_Untreated1.sam | samtools sort -@ 8 -o GSE164704_Untreated1.sorted.bam
samtools view -@ 8 -bS GSE164704_Untreated2.sam | samtools sort -@ 8 -o GSE164704_Untreated2.sorted.bam
samtools view -@ 8 -bS GSE164704_Untreated3.sam | samtools sort -@ 8 -o GSE164704_Untreated3.sorted.bam
samtools view -@ 8 -bS GSE164704_Untreated4.sam | samtools sort -@ 8 -o GSE164704_Untreated4.sorted.bam
samtools view -@ 8 -bS GSE164704_Flagellin_Treated1.sam | samtools sort -@ 8 -o GSE164704_Treated1.sorted.bam
samtools view -@ 8 -bS GSE164704_Flagellin_Treated2.sam | samtools sort -@ 8 -o GSE164704_Treated2.sorted.bam
samtools view -@ 8 -bS GSE164704_Flagellin_Treated3.sam | samtools sort -@ 8 -o GSE164704_Treated3.sorted.bam
samtools view -@ 8 -bS GSE164704_Flagellin_Treated4.sam | samtools sort -@ 8 -o GSE164704_Treated4.sorted.bam
samtools index GSE164704_Untreated2.sorted.bam 
samtools index GSE164704_Untreated3.sorted.bam 
samtools index GSE164704_Untreated4.sorted.bam 
samtools index GSE164704_Treated1.sorted.bam 
samtools index GSE164704_Treated2.sorted.bam 
samtools index GSE164704_Treated3.sorted.bam 
samtools index GSE164704_Treated4.sorted.bam

# TEtranscripts:
TEtranscripts  --format BAM  --mode multi  --GTF gencode.v38.annotation.gtf --TE hg38_rmsk_TE.gtf  --project Flagellin_L1 -t GSE164704_Treated1.sorted.bam  GSE164704_Treated2.sorted.bam  GSE164704_Treated3.sorted.bam  GSE164704_Treated4.sorted.bam  -c GSE164704_Untreated1.sorted.bam  GSE164704_Untreated2.sorted.bam  GSE164704_Untreated3.sorted.bam GSE164704_Untreated4.sorted.bam 
