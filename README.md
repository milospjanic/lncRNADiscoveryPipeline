# lncRNADiscoveryPipeline

The lncRNADiscoveryPipeline is a bash pipeline for the discovery of novel and quantification of existing lncRNAs and mRNAs. 

First part of the pipe will take all bam files from your folder, merge them into a single bam file, sort and index it and then run Stringtie to assemble transcripts with the -G option to provide the reference gtf file. 

The second part of the pipe will take bams per condition and merge sort and index them, then run stringtie for each condition in the experiment with the -G option to provide the reference gtf file.

Thrid part will merge gtf files from per-condition analysis into a single gtf and then merge the output gtf with the gtf from merged-bam analysis. The resulting gtf output will contain de novo lncRNAs, known lncRNAs and known mRNA genes. Such lncRNA pipe will detect both **condition specific lncRNAs** captured using per-condition bam merge strategy, and **lncRNA that are expressed at low levels** accross all conditions detected using merge all bam strategy that will increase the overall number of sequence reads and the probality that you will detect low level lncRNA, hence the dual approach of this strategy.

![alt text](lncRNADiscoveryPipeline/Screen Shot 2019-02-01 at 1.01.43 PM.png)
