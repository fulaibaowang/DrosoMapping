#!/usr/bin/Rscript
# module load math/R/4.1.2

#install.packages('fs',lib ='/home/fr/fr_fr/fr_yw1014/R/x86_64-pc-linux-gnu-library/4.1.2' )
#install.packages('tidyverse',lib ='/home/fr/fr_fr/fr_yw1014/R/x86_64-pc-linux-gnu-library/4.1.2' )
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager",lib ='/home/fr/fr_fr/fr_yw1014/R/x86_64-pc-linux-gnu-library/4.1.2')
#BiocManager::install(version = "3.14",lib ='/home/fr/fr_fr/fr_yw1014/R/x86_64-pc-linux-gnu-library/4.1.2')

#library("BiocManager",lib ='/home/fr/fr_fr/fr_yw1014/R/x86_64-pc-linux-gnu-library/4.1.2')
#BiocManager::install("karyoploteR",lib ='/home/fr/fr_fr/fr_yw1014/R/x86_64-pc-linux-gnu-library/4.1.2')
#BiocManager::install("Rsamtools",lib ='/home/fr/fr_fr/fr_yw1014/R/x86_64-pc-linux-gnu-library/4.1.2')

.libPaths('/home/fr/fr_fr/fr_yw1014/R/x86_64-pc-linux-gnu-library/4.1.2')

library(tidyverse)
library(fs)
library(karyoploteR)
library(Rsamtools)


setwd("/pfs/work7/workspace/scratch/fr_yw1014-minaproject/mapping_on_fly")

# read the genome file
seqs <- readDNAStringSet("/pfs/work7/workspace/scratch/fr_yw1014-minaproject/DrosoMapping/data/holo_dmel_6.12.fa")
# show the names of the chromosomes :
names(seqs)
# show the length of every chromosome
width(seqs)
# so we will use theses data to build a custom genome for our sequence using toGranges
# pick only the first 7 chromosomes "2L" "2R" "3L" "3R" "4"  "X"  "Y"
custom.genome <- toGRanges(data.frame(chr=c(names(seqs[1:7])), start=c(1), end=c(width(seqs[1:7]))))

# go the folder where our bam files are
file_paths <- fs::dir_ls("./", glob='*.bam')
file_paths

# then we make the plot for each bam file
for (i in 1:length(file_paths)){
  t=file_paths[i]
  pdf(file = paste('plot',substring(paste(t),1,13),'pdf',sep='.') ,width = 7,height = 10)
  kp <- plotKaryotype(genome=custom.genome,main=substring(paste(t),1,13))
#kp <- kpPlotBAMCoverage(kp, data=i,max.valid.region.size=1e10)         
  kp <- kpPlotBAMDensity(kp, data=t)
#kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage)
  kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density)
  dev.off()
}

