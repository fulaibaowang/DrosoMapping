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
  -o trimmed/
  untrimmed/${i}_R1.fq.gz
  untrimmed/${i}_R2.fq.gz

### 2) mapping with minimap

# a) download precompiled binary

curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -

# b) index reference


minimap2-2.17_x64-linux/
