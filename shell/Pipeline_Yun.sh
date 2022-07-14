## mapping with bowtie

# use bowtie2-build to generate the genome database
~/fr_yw50/downloads/bowtie2-2.4.2-linux-x86_64/bowtie2-build DrosoMapping/data/holo_dmel_6.12.fa.gz ./holodmel

# in cluster, generate a .sh1 file for each sample
cd raw_flyreads/
unset listA
unset listB
listA=(*_1.*q.gz)
listB=(*_2.*q.gz)
n=${#listA[@]}

for i in $(seq 0 $(($n-1))); do echo -e '#!/bin/bash\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task=16\n#SBATCH --time=0:59:59\n#SBATCH --output=bowtie_out_'${listA[$i]}'.txt
date\ncd /pfs/work7/workspace/scratch/fr_yw1014-minaproject/raw_flyreads\n ~/fr_yw50/downloads/bowtie2-2.4.2-linux-x86_64/bowtie2 --threads 16 -x ../holodmel -1 '${listA[$i]}' -2 '${listB[$i]}' -S ../mapping_on_fly/Sample_'${listA[$i]}'.sam\ndate'  >  ${listA[$i]}.sh1; done

# in each file, it look something like this
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

## convert sam -> bam, sort and index bam files

# here is the job I submitted in cluster

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
