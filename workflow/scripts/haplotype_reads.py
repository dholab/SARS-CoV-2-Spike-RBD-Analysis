import subprocess
import math

# calculate number of reads for 1% of consensus_min_depth
haplotype_threshold = math.floor(
    float(snakemake.config["thresholds"]["target_depth"])
    * float(snakemake.params["haplo_freq"])
)

subprocess.call(
    "vsearch --fastx_uniques "
    + snakemake.input[0]
    + " -fastqout "
    + snakemake.output[0]
    + " -sizeout --minuniquesize "
    + str(haplotype_threshold)
    + " -strand both",
    shell=True,
)
