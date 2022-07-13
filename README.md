## Introduction

This workflow processes paired-end Illumina SARS-CoV-2 RBD amplicon sequencing data collected from air samples in [Ramuta et al. 2022, SARS-CoV-2 and other respiratory pathogens are detected in continuous air samples from congregate settings](https://doi.org/10.1101/2022.03.29.22272716).


## Getting started

This workflow uses the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow manager. If you don't have `snakemake` installed on your computer, you need to follow five steps:

1. Install the miniconda python distribution (https://docs.conda.io/en/latest/miniconda.html)
2. Install the `mamba` package installation tool:
`conda install -y -c conda-forge mamba`
3. Install `snakemake` into a new python virtual environment:
` mamba create -c conda-forge -c bioconda -n snakemake snakemake`
4. Activate the new environment:
`conda activate snakemake`

Try typing `snakemake -v`. If your terminal returns something like `6.12.3` snakemake is installed correctly.


### Bundled data

- SARS-CoV-2 RBD amplicon BED files from [Dr. Marc Johnson](https://doi.org/10.3390/v13081647) are in the `resources/amplicon` folder.
- SARS-CoV-2 RBD primer scheme from Dr. Marc Johnson is in the `resources/primers` folder. 
- An example SARS-CoV-2 genome reference in the `resources/genomes` folder.

## Workflow summary
- This workflow supports multiple interleaved FASTQ file inputs
- Paired-end reads are merged into synthetic reads spanning the entire RBD PCR amplicon
- Synthetic reads are mapped to a viral reference sequence
- Mapped read positions are compared to a BED file of amplicon sequences. The BED file shows amplicon coordinates after PCR primer trimming, representing the "full-length" amplicon sequences. Synthetic reads that span the entire amplicon are retained.
- These reads are downsampled on a per-amplicon basis
- The downsampled reads are mapped to a viral reference sequence and a consensus sequence is generated. The mapping BAM file can be loaded into Geneious or other tools for variant calling.
- Downsampled reads are deduplicated and those deduplicated clusters that exceed a minimum frequency are returned as putative haplotype sequences. These can be used in downstream analyses to determine linkage of variants within the same amplicon.

## Configuration parameters

`reads` - path to interleaved Illumina FASTQ files (expects .gz compressed FASTQ)
`amplicon_bed` - path to BED file that contains one row per amplicon, with coordinates of amplicon relative to reference after primers have been removed.
`primer_bed` - path to BED file that contains one row per primer pair.
`ref_fasta` - path to FASTA file with viral reference genome
`target_depth` - number of reads to retain from each amplicon. For example, '1000' would retain a maximum of 1,000 reads from each amplicon specified in `amplicon_bed`
`output_prefix` - string to use as the prefix for all output files
`consensus_min_depth` - when generating a consensus sequence, the minimum coverage required for a consensus basecall. For example, if this is set to '20' sites that have coverage of '12' would be masked with an 'N'
`haplotype_min_frequency` - frequency of deduplicated reads to report in haplotypes. For example, setting this to '0.01' and 'target_depth' to 1000 would retain sequences that are found in at least 10 identical downsampled reads (1%). Note that this expects that there are 'target_depth' reads for each amplicon; amplicons with unexpectedly low coverage will still report sequences found in at least 10 identical reads.

## Configuration parameters set in the workflow

`reads` - 'data/'
`amplicon_bed` - 'resources/amplicons/SARS-CoV-2.MJ-RBD-NTD.amplicon.bed'
`primer_bed` - 'resources/primers/SARS-CoV-2.MJ-RBD-NTD.primer.bed'
`ref_fasta` - 'resources/genomes/MN908947.3.fa'
`target_depth` - '1000'
`consensus_min_depth` - '20'
`haplotype_min_frequency` - '0.001'

## Output files

- `[output_prefix].[target_depth]X.fastq.gz` - FASTQ files after downsampling to target depth
- `[output_prefix].[target_depth]X.primer_removed.sorted.bam` - BAM file after mapping downsampled FASTQ to `ref_fasta` and removing primers specified in `primer_bed`
- `[output_prefix].[target_depth]X.primer_removed.sorted.bam.bai` - BAM index of the mapped file
- `[output_prefix].[target_depth]X..primer_removed.consensus.fa` - FASTA file with consensus sequence derived from BAM file
- `[output_prefix].[target_depth]X..primer_removed.consensus.qual.txt` - Text file with consensus sequence quality scores derived from BAM file
- `[output_prefix].[target_depth]X.haplotypes.fastq` - FASTQ file where identical reads at a frequency greater than `haplotype_min_frequency` are retained. Useful for mapping to a reference, extracting all reads corresponding to a single amplicon, and identifying sequences corresponding to haplotypes. 

## Example data

The archive includes interleaved sequencing data for 9 air samples from Ramuta et al. 2022 in the `data` folder.

## Analysis

Run snakemake with: `snakemake --use-conda -- cores 8`

## For more information

Contact Dr. David O'Connor (dhoconno@wisc.edu), Dr. Shelby O'Connor (slfeinberg@wisc.edu), or Mitchell Ramuta (ramuta2@wisc.edu)