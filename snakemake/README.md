#snakemake notes by Yun

I played with snakemake with the fly data and created a snakemake workflow basically the same as in the script.

The workflow is in Snakefile.

Some key facts:

0. I did it in conda environment and the yml file is available. Java version has to be 8 not newer.

1. The Snakefile was tested on first 2 samples and looked good.

2. It is possible to create a workflow figure with the command: snakemake --dag all | dot -Tpng > dag.png

3. when running MarkDuplicatesSpark, there was a mixmal open file number issue, can be fixed when setting a higher number: ulimit -n 3000. There was also a java memory issue when setting the flag "-Xmx2G", and when setting "-Xmx4G" it worked fine.
