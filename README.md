# DrosoMapping by Martin

Handcrafted for Fabian Staubach!

See [shell/Pipeline.sh](shell/Pipeline.sh) for the full Pipelines


## comments add by Yun

See [shell/Pipeline_Yun.sh](shell/Pipeline_Yun.sh) for the new script


What was done so far:

1. reads mappingn on fly genome via bowtie2

Additionally, I created this R script [scripts/coverageplot_severalfiles_cluster.R](scripts/coverageplot_severalfiles_cluster.R) to visulaize the reads coverage of all bam files. The plots are here [data/bam reads coverage plots](data/bam reads coverage plots) (will add this later)

2. convert to sorted and indexed bam files

3. Pro-process the bam files following https://gatk.broadinstitute.org/hc/en-us/articles/360035535912

To do list:

SNP calling by GATK, https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

toubleshoot samples with different between in R1 and R2?