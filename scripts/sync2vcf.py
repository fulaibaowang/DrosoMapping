from datetime import date
from optparse import OptionParser, OptionGroup
import sys
from collections import defaultdict
from operator import itemgetter
import gzip

# Author: Martin Kapun

# version 1.1 # 05/04/2018

parser = OptionParser()
group = OptionGroup(
    parser,
    """
""",
)

parser.add_option(
    "--reference",
    dest="ref",
    help="The name of the reference genome version used for alignemnt to be stored in the header",
    default="dmel v.5",
)
parser.add_option(
    "--source",
    dest="source",
    help="The name of the program used for SNP calling to be stored in the header",
    default="PoPoolation2",
)
parser.add_option(
    "--names",
    dest="names",
    help="comma-delimited list of the names of all libraries in the sync file",
)
parser.add_option("--input", dest="input", help="input sync file")
parser.add_option("--output", dest="output", help="output vcf file")
parser.add_option(
    "--biallelic",
    action="store_true",
    dest="biallelic",
    help="only consider the two most common alleles",
)

parser.add_option_group(group)
(options, args) = parser.parse_args()


def load_data(x):
    """ import data either from a gzipped or or uncrompessed file or from STDIN"""
    import gzip

    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


now = date.today()

Date = str(now).replace("-", "")
ref = options.ref
source = options.source
names = options.names.split(",")
if ".gz" in options.output:
    out = gzip.open(options.output, "w")
else:
    out = open(options.output, "w")
out.write("##fileformat=VCFv4.2\n")
out.write("##fileDate=" + Date + "\n")
out.write("##Source=" + source + "\n")
out.write("##reference=" + ref + "\n")
out.write('##FORMAT=<ID=AD,Number=1,Type=String,Description="Alternative Counts">\n')
out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
out.write(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
    "\t".join(names) + "\n"
)

alleles = ["A", "T", "C", "G"]

for l in load_data(options.input):
    a = l.split()
    CHR, POS, REF = a[:3]
    pops = a[3:]
    # find two most commmon alleles
    tot = []
    totalleles = defaultdict(int)
    for pop in pops:
        if ":" not in pop:
            continue
        allele = pop.split(":")[:4]
        for p in range(4):
            if allele[p] == "0" or alleles[p] == REF:
                continue
            tot.append(alleles[p])
            totalleles[alleles[p]] += int(allele[p])
    if len(totalleles) == 0:
        continue
    if options.biallelic:
        alt = list(zip(*sorted(totalleles.items(), key=itemgetter(1), reverse=True)))[
            0
        ][0]
    else:
        alt = sorted(set(tot))
    counts = []
    for pop in pops:
        if ":" not in pop:
            continue
        allele = pop.split(":")[:4]
        DP = sum([int(x) for x in allele])
        allelehash = dict(zip(alleles, allele))
        AP = [allelehash[x] for x in alt]
        counts.append(",".join([str(x) for x in AP]) + ":" + str(DP))

    out.write(
        CHR
        + "\t"
        + POS
        + "\t.\t"
        + REF
        + "\t"
        + ",".join(alt)
        + "\t.\t.\t.\tAD:DP\t"
        + "\t".join(counts)
        + "\n"
    )
