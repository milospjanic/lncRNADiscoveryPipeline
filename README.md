# lncRNADiscoveryPipeline

The lncRNADiscoveryPipeline is a bash pipeline for the discovery of novel and quantification of existing lncRNAs and mRNAs. 

First part of the pipe will take all bam files from your folder, merge them into a single bam file, sort and index it and then run Stringtie to assemble transcripts with the -G option to provide the reference gtf file. 

The second part of the pipe will take bams per condition and merge sort and index them, then run stringtie for each condition in the experiment with the -G option to provide the reference gtf file.

Thrid part will merge gtf files from per-condition analysis into a single gtf and then merge the output gtf with the gtf from merged-bam analysis. The resulting gtf output will contain de novo lncRNAs, known lncRNAs and known mRNA genes. Such lncRNA pipe will detect both **condition specific lncRNAs** captured using per-condition bam merge strategy, and **lncRNA that are expressed at low levels** accross all conditions detected using merge all bam strategy that will increase the overall number of sequence reads and the probality that you will detect low level lncRNA, hence the dual approach of this strategy.

![alt text](https://github.com/milospjanic/lncRNADiscoveryPipeline/blob/master/scheme.png)

# Requirements

Stringtie
Samtools

# Example run

Bam files must be in a same folder with *clean.bam* extension

<pre>
ls *clean.bam
CtrlCyto-1_clean.bam  CtrlRna-1_clean.bam  PDGFD-1_clean.bam  Serum-1_clean.bam  SMAD3KD-1_clean.bam  TCF21KD-1_clean.bam  TGFb-1_clean.bam  TNFa-1_clean.bam
CtrlCyto-2_clean.bam  CtrlRna-2_clean.bam  PDGFD-2_clean.bam  Serum-2_clean.bam  SMAD3KD-2_clean.bam  TCF21KD-2_clean.bam  TGFb-2_clean.bam  TNFa-2_clean.bam
CtrlCyto-3_clean.bam  CtrlRna-3_clean.bam  PDGFD-3_clean.bam  Serum-3_clean.bam  SMAD3KD-3_clean.bam  TCF21KD-3_clean.bam  TGFb-3_clean.bam  TNFa-3_clean.bam
CtrlCyto-4_clean.bam  CtrlRna-4_clean.bam  PDGFD-4_clean.bam  Serum-4_clean.bam  SMAD3KD-4_clean.bam  TCF21KD-4_clean.bam  TGFb-4_clean.bam  TNFa-4_clean.bam
CtrlCyto-5_clean.bam  CtrlRna-5_clean.bam  PDGFD-5_clean.bam  Serum-5_clean.bam  SMAD3KD-5_clean.bam  TCF21KD-5_clean.bam  TGFb-5_clean.bam  TNFa-5_clean.bam
CtrlCyto-6_clean.bam  CtrlRna-6_clean.bam  PDGFD-6_clean.bam  Serum-6_clean.bam  SMAD3KD-6_clean.bam  TCF21KD-6_clean.bam  TGFb-6_clean.bam  TNFa-6_clean.bam
</pre>



