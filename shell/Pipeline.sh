### mapping pipeline for the FAB ;-)

## 1) trim with trimgalore

# a) install trim_galore:

conda create --name trim-galore-0.6.2 -c bioconda trim-galore
conda activate trim-galore-0.6.2

# b) Trim paired end reads

# -q 20         trim basequality <20
# --fastqc      run also FASTQC
# --gzip        compress output
# --length 75   discarD trimmed reads with length <75bp
# -o            output directory
# -j 20         20 cores to be used
# --paired      paired-end mode

mkdir trimmed/

# loop through all samples
for i in sample1 sample2 sample3 sample4

do

  trim_galore \
    -q 20 \
    --fastqc \
    --gzip \
    --length 75 \
    -o trimmed/ \
    untrimmed/${i}_R1.fq.gz \
    untrimmed/${i}_R2.fq.gz

done

### 2) mapping with minimap

# a) download precompiled minimap2 binary

curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -

# b) install samtools

wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
tar -jxvf - samtools-1.12.tar.bz2
cd samtools-1.12
./configure
make

# b) index reference

minimap2-2.17_x64-linux/minimap2 index data/holo_dmel6.12.fa.gz

# c) map all samples against reference and sort, discard unmapped reads and remove reads with mapping quality <20 with samtools. Then index the BAM file

for i in sample1 sample2 sample3 sample4

do

  ## minimap parameters
  # -a          output CIGAR String in BAM file format
  # -x sr       input are short reads
  # -t 20       20 cores

  ## samtools parameters
  # -h          print BAM header
  # -b          output in BAM file format
  # -q 20       discard reads with mapping quality <20

  minimap2 \
    -a \
    -x sr \
    -t 20 \
    data/holo_dmel6.12.fa.gz \
    untrimmed/${i}_R1.fq.gz \
    untrimmed/${i}_R2.fq.gz \
    | samtools view -b \
    -h \
    -q 20 \
    | samtools sort > ${i}.bam

  samtools index ${i}.bam

done

## 3) download SNP dataset from DEST

wget http://berglandlab.uvadcos.io/vcf/dest.all.PoolSNP.001.50.10Nov2020.ann.vcf.gz

## 4) convert to BED file format
gunzip -c dest.all.PoolSNP.001.50.10Nov2020.ann.vcf.gz \
  | grep -v "^#" \
  | awk '{print $1"\t"$2-1"\t"$2}' \
  > dest.all.PoolSNP.bed

## 5) synchronice BAM files and isolate SNP positions based on BED file

# a) make list of BAM files

for i in sample1 sample2 sample3 sample4

do

  echo ${i}.bam >> bamlist.txt

done

# b) synchronize BAM file and convert to VCF file

samtools mpileup \
  -B \
  -f data/holo_dmel6.12.fa.gz \
  -b bamlist.txt
