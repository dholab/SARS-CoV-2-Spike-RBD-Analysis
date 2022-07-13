import subprocess

# read amplicon BED line-by-line
with open(snakemake.input[1], "r") as a_file:
    for ct, line in enumerate(a_file):
        # write single line BED
        with open(snakemake.output[2], "w") as bed_file:
            bed_file.write(line)

        # make BAM of reads that span amplicon
        # first extract potentially matching reads with samtools view
        # this reduces the search space for the slower bedtools intersect program
        # then find sequences that fully span the interval with bedtools
        subprocess.call(
            "samtools view -b "
            + snakemake.input[0]
            + " --target-file "
            + snakemake.output[2]
            + " \
        | bedtools intersect -ubam -F 1 -a stdin -b "
            + snakemake.output[2]
            + " > "
            + snakemake.output[1],
            shell=True,
        )

        # downsample and create FASTQ.gz
        subprocess.call(
            "reformat.sh app=t in="
            + snakemake.output[1]
            + " ref="
            + snakemake.input[2]
            + " out="
            + snakemake.output[0]
            + " samplereadstarget="
            + snakemake.params.td,
            shell=True,
        )
