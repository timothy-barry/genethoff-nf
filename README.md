---
editor_options: 
  markdown: 
    wrap: 72
---

# Donor-seq pipeline

This is a pipeline for processing donor-seq data. Donor-seq is an assay
for profiling the off-target editing activity of genome editors based on
prime assembly (introduced in [Levesque et al
2025](https://www.biorxiv.org/content/10.1101/2025.06.16.659926v1)).
This pipeline is a reworked, modified version of the
[GENETHOFF](https://github.com/gcorre/GENETHOFF/releases) GUIDE-seq
pipeline.

## Requirements and installation

Install `Nextflow` using the instructions
[here](https://www.nextflow.io/docs/latest/install.html). Next, set up
Nextflow for your cluster using the instructions
[here](https://timothy-barry.github.io/sceptre-book/pipeline-args.html).
Install `conda` using the instructions
[here](https://docs.conda.io/projects/conda/en/stable/index.html).
Installing `Nextflow` and `conda` typically takes less than 20 minutes.

The pipeline depends on a `conda` environment called
`guideseq-pipeline`. The `.yml` file for `guideseq-pipeline` is pasted
below (and is also present in this repository).

```         
name: guideseq-pipeline
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  # Bioinformatics
  - cutadapt=5.1
  - bowtie2=2.5.4
  - samtools=1.22.1
  - bedtools=2.31.1
  # R
  - r-base=4.3
  - r-dplyr=1.1.4
  # Python
  - python=3.11
```

Nextflow will try to automatically create this `conda` environment using
the `.yml` file; thus, you should not need to manually download the
software packages listed above.

Finally, ensure that you have a reference genome (typically GRCh38)
installed and indexed via `bowtie2`. For the best results, include only
the main chromosomes from the primary assembly. The directory containing
your reference genome should have the following files. (The file base
name need not be `hg38_main_chroms`.)

```         
hg38_main_chroms.1.bt2     hg38_main_chroms.3.bt2     hg38_main_chroms.fa        hg38_main_chroms.rev.1.bt2
hg38_main_chroms.2.bt2     hg38_main_chroms.4.bt2     hg38_main_chroms.fa.fai    hg38_main_chroms.rev.2.bt2
```

**Tip (optional)**: To accelerate installation of the `conda`
environment, install `mamba`, a more efficient version of `conda`.
`mamba` is included as part of most modern `conda` installations. Within
your global `Nextflow` config file (typically stored in
`~/.nextflow/config`), add the line `conda.useMamba = true` to use
`mamba` rather than `conda` within Nextflow.

## Software testing

The pipeline has been tested on a SLURM cluster running Red Hat
Enterprise Linux version 9.6 and a Macbook Pro running macOS Sonoma
version 14.5 using `Nextflow` version `25.10.0` and `conda` version
`24.11.3`.

## Demo

A demo dataset is contained within this repository. To access the demo
dataset, first clone the repository.

```         
git clone git@github.com:timothy-barry/genethoff-nf.git
```

Next, `cd` into the `genethoff-nf/demo` directory. This directory
contains several data files and scripts. The data files are as follows:

```         
Jing_AAVS1_n1-10_Donor-Seq_AAVS1_GSP_plus_1_I2_first10000.fastq
Jing_AAVS1_n1-10_Donor-Seq_AAVS1_GSP_plus_1_R1_first10000.fastq
Jing_AAVS1_n1-10_Donor-Seq_AAVS1_GSP_plus_1_R2_first10000.fastq
Jing_AAVS1_n1-10_Donor-Seq_IL2RG_GSP_minus_1_I2_first10000.fastq
Jing_AAVS1_n1-10_Donor-Seq_IL2RG_GSP_minus_1_R1_first10000.fastq
Jing_AAVS1_n1-10_Donor-Seq_IL2RG_GSP_minus_1_R2_first10000.fastq
```

These files correspond to two donor-seq samples:
`Jing_AAVS1_n1-10_Donor-Seq_AAVS1_GSP_plus` and
`Jing_AAVS1_n1-10_Donor-Seq_IL2RG_GSP_minus`. Each sample has an
associated R1, R2, and I2 file; R1 and R2 contain the paired-end reads,
while I2 contains the UMI information. The pipeline assumes that the UMI
is contained in the first nine bases of the I2 reads (although this
default behavior can be modified).

Next, `datasheet.tsv` is a tab-separated file storing sample metadata,
with each row corresponding to a different sample. Columns `R1`, `R2`,
and `I2` store the file paths to the R1, R2, and I2 files, respectively.
(*Be sure to update these files paths*!) Further, `sample_id` is a
string uniquely identifying the sample. Next, `negative_R2_leading` is
the technical oligo at the beginning of the R2 read when using the
negative primer; `positive_R2_leading` is the technical oligo at the
beginning of the R2 read when using the positive primer; and
`negative_R1_trailing` and `positive_R1_trailing` are the reverse
complements of `negative_R2_leading` and `positive_R2_leading`,
respectively.

Next, `nextflow.config` is a config file storing the pipeline
parameters. Briefly, `samplesheet`, `outdir`, and `index` store file
paths to the samplesheet, results directory, and genome index directory,
respectively. (*Be sure to update these files paths*!) Next,
`min_length`, `min_frag_length`, `max_frag_length`, `min_mapq`, and
`keep_multimapped_reads` are alignment parameters. Furthermore,
`umi_length` and `umi_side` indicate where the pipeline is to look for
the UMI within the I2 file. (`umi_length` is the length of the UMI,
while `umi_side` is either `5` or `3`, indicating whether the UMI is on
the 5' or 3' end of the read.) Finally, `R2_trailing` is the i5 Illumina
adapter sequence. (In most GUIDE-seq protocols, the i5 Illumina adapter
sequence `CTGTCTCTTATACAC`.)

Finally, `launch.sh` is the launch script. (*Be sure to update the file
path of the `.yml`* *file!*). When running the pipeline locally, launch
the pipeline via `bash launch_nf_pipeline.sh`; when running the pipeline
on a SLURM cluster, launch the pipeline via
`sbatch launch_nf_pipeline.sh`, etc. Resource allocation requests for
the `Nextflow` driver should be issued at the top of this script.

## Expected output

After running the pipeline, outputs should be written to the directory
`genethoff-nf/demo/results`. The results directory contains two
subdirectories --- namely, `Jing_AAVS1_n1-10_Donor-Seq_AAVS1_GSP_plus_1`
and `Jing_AAVS1_n1-10_Donor-Seq_IL2RG_GSP_minus_1` --- corresponding to
the two samples in the samplesheet. Each subdirectory contains `*.log`
files containing plain-text information about how the pipeline ran and
an `.rds` file storing the counts data frame. The `.rds` file can be
opened in R via `readRDS()`. Each row corresponds to an "occupied base",
or a base to which at least one donor-seq read mapped. The columns are
as follows:

-   `chr`: chromosome

-   `coord`: coordinate of the chromosome

-   `strand`: whether the donor-seq read mapped to the + or - strand

-   `umi_count`: UMI count of the base

-    `total_read_count`: total read count across the UMIs

-   `mean_mapq`: mean mapq score of the reads (before collapsing by UMI)

## Expected runtime

When running this pipeline for the first time, Nextflow automatically
will try to create a `guideseq-pipeline` `conda` environment. When using
`mamba` (see "Requirements and installation"), the environment creation
step should take under 10 minutes. The pipeline itself should run in
under two minutes (depending on how long it takes for jobs to get
scheduled).
