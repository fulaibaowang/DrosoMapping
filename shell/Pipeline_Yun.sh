## 1 mapping with bowtie

# 1.1 use bowtie2-build to generate the genome database
~/fr_yw50/downloads/bowtie2-2.4.2-linux-x86_64/bowtie2-build DrosoMapping/data/holo_dmel_6.12.fa.gz ./holodmel

# 1.2 do bowtie2 mapping for each sample in a parallel fashion
# in cluster, generate a .sh1 file for each sample by these commands:

cd raw_flyreads/
unset listA
unset listB
listA=(*_1.*q.gz)
listB=(*_2.*q.gz)
n=${#listA[@]}
for i in $(seq 0 $(($n-1))); do echo -e '#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=16\n#SBATCH --time=0:59:59\n#SBATCH --output=bowtie_out_'${listA[$i]}'.txt
date\ncd /pfs/work7/workspace/scratch/fr_yw1014-minaproject/raw_flyreads\n ~/fr_yw50/downloads/bowtie2-2.4.2-linux-x86_64/bowtie2 --threads 16 -x ../holodmel -1 '${listA[$i]}' -2 '${listB[$i]}' -S ../mapping_on_fly/Sample_'${listA[$i]}'.sam\ndate'  >  ${listA[$i]}.sh1; done

# in each file, it look something like this:

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=0:59:59
#SBATCH --output=bowtie_out_DE3035_1.fastq.gz.txt
date
cd /pfs/work7/workspace/scratch/fr_yw1014-minaproject/raw_flyreads
 ~/fr_yw50/downloads/bowtie2-2.4.2-linux-x86_64/bowtie2 --threads 16 -x ../holodmel -1 DE3035_1.fastq.gz -2 DE3035_2.fastq.gz -S ../mapping_on_fly/Sample_DE3035_1.fastq.gz.sam
date

# submit the sh1 files in cluster
for f in *sh1; do sbatch -p single $f;done 


## 2 convert sam -> bam, sort and index bam files

# here is the job I submitted in cluster:

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=20:59:59
date
cd /pfs/work7/workspace/scratch/fr_yw1014-minaproject/mapping_on_fly
for x in *.sam;
do ~/fr_yw50/downloads/samtools1.11/bin/samtools view -F 4 -bS $x > $x-RAW.bam;
done
#sort
for x in *RAW.bam;
do ~/fr_yw50/downloads/samtools1.11/bin/samtools sort $x > $x-sort.bam;
done
#index
for x in *RAW.bam;
do ~/fr_yw50/downloads/samtools1.11/bin/samtools index $x-sort.bam;
done
date


# 3 preprocess bam files with GATK

# 3.1 add labels of read groups using picard.jar AddOrReplaceReadGroups
# this is required for furthur analysis
# https://gatk.broadinstitute.org/hc/en-us/articles/5358911906459-AddOrReplaceReadGroups-Picard-

cd /pfs/work7/workspace/scratch/fr_yw1014-minaproject/mapping_on_fly
unset listA
listA=(*.bam)
n=${#listA[@]}

for i in $(seq 0 $(($n-1))); do str=${listA[$i]}; echo -e '#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=4\n#SBATCH --time=0:39:59\n#SBATCH --output=AddOrReplaceReadGroups_out_'${listA[$i]}'.txt
date\ncd /pfs/work7/workspace/scratch/fr_yw1014-minaproject/mapping_on_fly\n java  -Xmx3g -jar /pfs/work7/workspace/scratch/fr_yw1014-minaproject/gatk/picard.jar AddOrReplaceReadGroups I='${listA[$i]}' O='${listA[$i]}'_addlabel.bam RGID='${str:7:6}' RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM='${str:7:6}'\ndate'  >  ${listA[$i]}.sh2; done

# each .sh2 file looks like this:
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0:59:59
#SBATCH --output=AddOrReplaceReadGroups_out_Sample_DE6176_1.fastq.gz.sam-RAW.bam-sort.bam.txt
date
cd /pfs/work7/workspace/scratch/fr_yw1014-minaproject/mapping_on_fly
 java  -Xmx3g -jar /pfs/work7/workspace/scratch/fr_yw1014-minaproject/gatk/picard.jar AddOrReplaceReadGroups I=Sample_DE6176_1.fastq.gz.sam-RAW.bam-sort.bam O=Sample_DE6176_1.fastq.gz.sam-RAW.bam-sort.bam_addlabel.bam RGID=DE6176 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=DE6176
date

# the label added is like this: RGID=DE6176 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=DE6176
# RGID & RGSM are sample specific; the other labels are all the same

# submit these .sh2 files simultaneously:
for f in DE30*sh2; do sbatch -p single $f;done 
for f in DE36*sh2; do sbatch -p single $f;done 
for f in DE37*sh2; do sbatch -p single $f;done 
for f in DE55*sh2; do sbatch -p single $f;done 
for f in DE56*sh2; do sbatch -p single $f;done 
for f in DE57*sh2; do sbatch -p single $f;done 
for f in DE58*sh2; do sbatch -p single $f;done 
for f in DE59*sh2; do sbatch -p single $f;done 
for f in DE60*sh2; do sbatch -p single $f;done 
for f in DE61*sh2; do sbatch -p single $f;done 

# 3.2 using GATK MarkDuplicatesSpark, output files handed to Mina
unset listA
listA=(*_addlabel.bam)
n=${#listA[@]}

for i in $(seq 0 $(($n-1))); do echo -e '#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=8\n#SBATCH --time=0:59:59\n#SBATCH --output=MarkDuplicatesSpark_out_'${listA[$i]}'.txt
date\ncd /pfs/work7/workspace/scratch/fr_yw1014-minaproject/mapping_on_fly\n /pfs/work7/workspace/scratch/fr_yw1014-minaproject/gatk/gatk-4.2.6.1/gatk --java-options "-Xmx7G" MarkDuplicatesSpark -I '${listA[$i]}' -O '${listA[$i]}'_duplicates.bam\ndate'  >  ${listA[$i]}.sh3; done

# each .sh3 file looks like this:
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0:59:59
#SBATCH --output=MarkDuplicatesSpark_out_Sample_DE6177_1.fastq.gz.sam-RAW.bam-sort.bam_addlabel.bam.txt
date
cd /pfs/work7/workspace/scratch/fr_yw1014-minaproject/mapping_on_fly
 /pfs/work7/workspace/scratch/fr_yw1014-minaproject/gatk/gatk-4.2.6.1/gatk --java-options "-Xmx7G" MarkDuplicatesSpark -I Sample_DE6177_1.fastq.gz.sam-RAW.bam-sort.bam_addlabel.bam -O Sample_DE6177_1.fastq.gz.sam-RAW.bam-sort.bam_addlabel.bam_duplicates.bam
date

# submit these .sh3 files simultaneously


# 3.3 BaseRecalibrator and ApplyBQSR

# do 3 additional steps for prepariation, output files provided in data folder
# a) first, index referece fasta file using samtools
~/fr_yw50/downloads/samtools1.11/bin/samtools faidx /pfs/work7/workspace/scratch/fr_yw1014-minaproject/DrosoMapping/data/holo_dmel_6.12.fa 

# b) get dict file from fasta file
path/to/gatk --java-options "-Xmx4G" CreateSequenceDictionary -R /pfs/work7/workspace/scratch/fr_yw1014-minaproject/DrosoMapping/data/holo_dmel_6.12.fa

# c)index VCF file
path/to/gatk --java-options "-Xmx4G" IndexFeatureFile -I ../dest.all.PoolSNP.001.50.10Nov2020.ann.vcf

# now, ready to do BaseRecalibrator and ApplyBQSR
unset listA
listA=(*_duplicates.bam)
n=${#listA[@]}

for i in $(seq 0 $(($n-1))); do echo -e '#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=8\n#SBATCH --time=03:59:59\n#SBATCH --output=BQSR_out_'${listA[$i]}'.txt
date\ncd /pfs/work7/workspace/scratch/fr_yw1014-minaproject/mapping_on_fly\n /pfs/work7/workspace/scratch/fr_yw1014-minaproject/gatk/gatk-4.2.6.1/gatk --java-options "-Xmx7G" BaseRecalibrator -I '${listA[$i]}' -R /pfs/work7/workspace/scratch/fr_yw1014-minaproject/DrosoMapping/data/holo_dmel_6.12.fa --known-sites ../dest.all.PoolSNP.001.50.10Nov2020.ann.vcf -O '${listA[$i]}'_recal_data.table
/pfs/work7/workspace/scratch/fr_yw1014-minaproject/gatk/gatk-4.2.6.1/gatk ApplyBQSR -R /pfs/work7/workspace/scratch/fr_yw1014-minaproject/DrosoMapping/data/holo_dmel_6.12.fa -I '${listA[$i]}' --bqsr-recal-file '${listA[$i]}'_recal_data.table -O '${listA[$i]}'_BQSR.bam\ndate'  >  ${listA[$i]}.sh4; done

# submit these .sh4 files simultaneously:
