# DrosoMapping by Martin

Handcrafted for Fabian Staubach!

See [shell/Pipeline.sh](shell/Pipeline.sh) for the full Pipelines


## comments add by Yun

See [shell/Pipeline_Yun.sh](shell/Pipeline_Yun.sh) for the new script


What was done so far:

0. toubleshoot samples with different reads between in R1 and R2, these samples should be redone
    list is here [data/R1R2.different.reads.list.txt](data/R1R2.different.reads.list.txt)

1. reads mappingn on fly genome via bowtie2

    Additionally, I created this R script [scripts/coverageplot_severalfiles_cluster.R](scripts/coverageplot_severalfiles_cluster.R) to visulaize the reads coverage of all bam files. The plots are here [data/bam_reads_coverage_plots/](data/bam_reads_coverage_plots/)

2. convert to sorted and indexed bam files

3. Pro-process the bam files following https://gatk.broadinstitute.org/hc/en-us/articles/360035535912

    The bam files handed to Mina are after step 3.2 not step 3.3. Step 3.3 did not finish (yet).

To do list:

SNP calling main step by GATK, https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

