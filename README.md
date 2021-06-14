# SARS-service_BU-ISCIII


[![GitHub Actions CI Status](https://github.com/nf-core/viralrecon/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/viralrecon/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/viralrecon/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/viralrecon/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3901628.svg)](https://doi.org/10.5281/zenodo.3901628)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/viralrecon.svg)](https://hub.docker.com/r/nfcore/viralrecon)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23viralrecon-4A154B?logo=slack)](https://nfcore.slack.com/channels/viralrecon)

## Introduction

**viralrecon_BU-ISCIII** is a bioinformatics analysis pipeline derived from nf-core/viralrecon one, which is used to perform assembly and intra-host/low-frequency variant calling for viral samples. This pipeline, similarly to what developed from nf-core/viralrecon, supports short-read Illumina sequencing data from both shotgun (e.g. sequencing directly from clinical samples) and enrichment-based library preparation methods (e.g. amplicon-based: [ARTIC SARS-CoV-2 enrichment protocol](https://artic.network/ncov-2019); or probe-capture-based). 

In order to make reproducible the analysis run by BU_ISCIII, here we describe step by step how the viralrecon_BU-ISCIII can be run. ??

## Viene de DSL2??# 
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with Docker containers making installation trivial and results highly reproducible. Furthermore, automated continuous integration tests that run the pipeline on a full-sized dataset using AWS cloud ensure that the code is stable.

## Local data organization summary

Once every service is ordered, a template of data structure is created, starting with a new folder with the name of the service, the date, the type of project (i.e. SARS-Cov) and the name of the researcher, typically "SRVxxxx_date_SARS_name-researcher". 
In this folder there will be other directories: 

ANALYSIS  DOC  RAW  REFERENCES  RESULTS  TMP

In ANALYSIS folder is where all the analysis process will run. In RAW, inside a folder named with the RUN  there are the R1/2 fastq.gz files, where they will be taken from in order to perform the analysis. DOC, REFERENCES, RESULTS and TMP are empty. 

**ANALYSIS_folder**

00-reads. 
Here there is a Lablog file which has the commands to be executed to make symbolic links to every RAW/RUN_xxx/R1/2.fastq.gz files to be analysed.

Lablog
It contains the commands necesarry to make the sample_id.txt file, which is derived from all the samples located in the RAW folder.
It also has the commands to create the other two folders in whicc the analysis will be splited:

1. Analysis_01 -> Lablog 
2. Where the lablog of the previous services has to be located. This lablog makes symbolic links with 00-reads folder and samples_id.txt file, which contain the raw data (Fastqc format) and sample names respectively, which are going to be analyzed.
3. 
  A. First, the lablog create a csv file format which is the samplesheet file, to be taken by viralrecon pipeline. 
  B. 

2. Analysis_02_Met. 


date_ANALYSIS01_AMPLICONS_HUMAN. 

    Lablog -> 1) Make symbolic links eith 00-reads folder and samples_id.txt file. 2º) In addition, it contains the command necessary to perform the samplesheet.csv file;3º) the main nextflow command to run the variant analysis of SARS using viralrecon and kraken2 (wich will be located in the _01_viral_recon_mapping.sh bash file. 4º) Copy the bash script "create_summary_report.sh" and the python script percentajeNs.py from the "proccessin_Data/Bioinformatics/services_and_colaborations" folder. *) There is a commented command in this lablog for the assembly, but we dont run it. 
    
    variants_table 
        Lablog -> It makes symbolic links to vcf.gz from variants/varsca2/SAMPLE_NAME.AF0.75.vcf.gz; merge all VCF files using bcftools and extract data using POS/REF/ALT from VCF files, in order to make a table for variants analysis. 
    date_viralrecon_mapping
.
├── assembly
│   ├── cutadapt
│   │   ├── fastqc
│   │   │   └── zips
│   │   └── log
│   └── kraken2
├── multiqc
│   └── multiqc_data
├── pipeline_info
├── preprocess
│   ├── fastp
│   │   ├── fastqc
│   │   │   └── zips
│   │   └── log
│   └── fastqc
│       └── zips
└── variants
    ├── bam
    │   ├── log
    │   ├── mosdepth
    │   │   ├── amplicon
    │   │   │   └── plots
    │   │   └── genome
    │   │       └── plots
    │   ├── mpileup
    │   ├── picard_metrics
    │   └── samtools_stats
    ├── bcftools
    │   ├── bcftools_stats
    │   ├── consensus
    │   │   └── base_qc
    │   ├── quast
    │   │   ├── aligned_stats
    │   │   ├── basic_stats
    │   │   ├── contigs_reports
    │   │   │   └── minimap_output
    │   │   ├── genome_stats
    │   │   └── icarus_viewers
    │   └── snpeff
    ├── intersect
    │   ├── 214704
    │   ├── 214721
    │   ├── 214724
    │   ├── 214760
    │   ├── 214763
    │   ├── 214765
    │   ├── 214766
    │   ├── 214780
    │   └── 214807
    ├── ivar
    │   ├── bcftools_stats
    │   ├── consensus
    │   │   └── base_qc
    │   ├── log
    │   ├── quast
    │   │   └── AF0.75
    │   │       ├── aligned_stats
    │   │       ├── basic_stats
    │   │       ├── contigs_reports
    │   │       │   └── minimap_output
    │   │       ├── genome_stats
    │   │       └── icarus_viewers
    │   └── snpeff
    └── varscan2
        ├── bcftools_stats
        ├── consensus
        │   └── base_qc
        ├── log
        ├── quast
        │   └── AF0.75
        │       ├── aligned_stats
        │       ├── basic_stats
        │       ├── contigs_reports
        │       │   └── minimap_output
        │       ├── genome_stats
        │       └── icarus_viewers
        └── snpeff

├── assembly
│   ├── cutadapt
│   │   ├── fastqc
│   │   │   ├── 214704_1.ptrim_fastqc.html
│   │   │   ├── 214704_2.ptrim_fastqc.html
│   │   │   ├── 214721_1.ptrim_fastqc.html
│   │   │   ├── 214721_2.ptrim_fastqc.html
│   │   │   ├── 214724_1.ptrim_fastqc.html
│   │   │   ├── 214724_2.ptrim_fastqc.html
│   │   │   ├── 214760_1.ptrim_fastqc.html
│   │   │   ├── 214760_2.ptrim_fastqc.html
│   │   │   ├── 214763_1.ptrim_fastqc.html
│   │   │   ├── 214763_2.ptrim_fastqc.html
│   │   │   ├── 214765_1.ptrim_fastqc.html
│   │   │   ├── 214765_2.ptrim_fastqc.html
│   │   │   ├── 214766_1.ptrim_fastqc.html
│   │   │   ├── 214766_2.ptrim_fastqc.html
│   │   │   ├── 214780_1.ptrim_fastqc.html
│   │   │   ├── 214780_2.ptrim_fastqc.html
│   │   │   ├── 214807_1.ptrim_fastqc.html
│   │   │   ├── 214807_2.ptrim_fastqc.html
│   │   │   └── zips
│   │   │       ├── 214704_1.ptrim_fastqc.zip
│   │   │       ├── 214704_2.ptrim_fastqc.zip
│   │   │       ├── 214721_1.ptrim_fastqc.zip
│   │   │       ├── 214721_2.ptrim_fastqc.zip
│   │   │       ├── 214724_1.ptrim_fastqc.zip
│   │   │       ├── 214724_2.ptrim_fastqc.zip
│   │   │       ├── 214760_1.ptrim_fastqc.zip
│   │   │       ├── 214760_2.ptrim_fastqc.zip
│   │   │       ├── 214763_1.ptrim_fastqc.zip
│   │   │       ├── 214763_2.ptrim_fastqc.zip
│   │   │       ├── 214765_1.ptrim_fastqc.zip
│   │   │       ├── 214765_2.ptrim_fastqc.zip
│   │   │       ├── 214766_1.ptrim_fastqc.zip
│   │   │       ├── 214766_2.ptrim_fastqc.zip
│   │   │       ├── 214780_1.ptrim_fastqc.zip
│   │   │       ├── 214780_2.ptrim_fastqc.zip
│   │   │       ├── 214807_1.ptrim_fastqc.zip
│   │   │       └── 214807_2.ptrim_fastqc.zip
│   │   └── log
│   │       ├── 214704.cutadapt.log
│   │       ├── 214721.cutadapt.log
│   │       ├── 214724.cutadapt.log
│   │       ├── 214760.cutadapt.log
│   │       ├── 214763.cutadapt.log
│   │       ├── 214765.cutadapt.log
│   │       ├── 214766.cutadapt.log
│   │       ├── 214780.cutadapt.log
│   │       └── 214807.cutadapt.log
│   ├── kraken2
│   │   ├── 214704.kraken2.report.txt
│   │   ├── 214721.kraken2.report.txt
│   │   ├── 214724.kraken2.report.txt
│   │   ├── 214760.kraken2.report.txt
│   │   ├── 214763.kraken2.report.txt
│   │   ├── 214765.kraken2.report.txt
│   │   ├── 214766.kraken2.report.txt
│   │   ├── 214780.kraken2.report.txt
│   │   └── 214807.kraken2.report.txt
│   └── summary_assembly_metrics_mqc.tsv
├── multiqc
│   ├── multiqc_data
│   │   ├── mqc_cutadapt_trimmed_sequences_plot_Counts.yaml
│   │   ├── mqc_cutadapt_trimmed_sequences_plot_Obs_Exp.yaml
│   │   ├── multiqc_bcftools_stats_bcftools_bcftools.yaml
│   │   ├── multiqc_bcftools_stats_bcftools_ivar.yaml
│   │   ├── multiqc_bcftools_stats_bcftools_varscan2.yaml
│   │   ├── multiqc_bowtie2.yaml
│   │   ├── multiqc_cutadapt.yaml
│   │   ├── multiqc_data.json
│   │   ├── multiqc_de_novo_assembly_metrics.yaml
│   │   ├── multiqc_fastp.yaml
│   │   ├── multiqc_fastqc_fastqc_cutadapt.yaml
│   │   ├── multiqc_fastqc_fastqc_fastp.yaml
│   │   ├── multiqc_fastqc_fastqc_raw.yaml
│   │   ├── multiqc_general_stats.yaml
│   │   ├── multiqc_ivar_primers.yaml
│   │   ├── multiqc_ivar_summary.yaml
│   │   ├── multiqc.log
│   │   ├── multiqc_picard_AlignmentSummaryMetrics.yaml
│   │   ├── multiqc_picard_insertSize.yaml
│   │   ├── multiqc_picard_wgsmetrics.yaml
│   │   ├── multiqc_quast_quast_bcftools.yaml
│   │   ├── multiqc_quast_quast_ivar.yaml
│   │   ├── multiqc_quast_quast_varscan2.yaml
│   │   ├── multiqc_samtools_flagstat_samtools_bowtie2.yaml
│   │   ├── multiqc_samtools_flagstat_samtools_ivar.yaml
│   │   ├── multiqc_samtools_idxstats_samtools_bowtie2.yaml
│   │   ├── multiqc_samtools_idxstats_samtools_ivar.yaml
│   │   ├── multiqc_samtools_stats_samtools_bowtie2.yaml
│   │   ├── multiqc_samtools_stats_samtools_ivar.yaml
│   │   ├── multiqc_snpeff_snpeff_bcftools.yaml
│   │   ├── multiqc_snpeff_snpeff_ivar.yaml
│   │   ├── multiqc_snpeff_snpeff_varscan2.yaml
│   │   ├── multiqc_sources.yaml
│   │   ├── multiqc_variant_calling_metrics.yaml
│   │   └── multiqc_varscan2_summary.yaml
│   └── multiqc_report.html
├── pipeline_info
│   ├── execution_report.html
│   ├── execution_timeline.html
│   ├── execution_trace.txt
│   ├── pipeline_dag.svg
│   ├── pipeline_report.html
│   ├── pipeline_report.txt
│   ├── results_description.html
│   ├── samplesheet.valid.csv
│   └── software_versions.csv
├── preprocess
│   ├── fastp
│   │   ├── 214704.fastp.html
│   │   ├── 214704.fastp.json
│   │   ├── 214721.fastp.html
│   │   ├── 214721.fastp.json
│   │   ├── 214724.fastp.html
│   │   ├── 214724.fastp.json
│   │   ├── 214760.fastp.html
│   │   ├── 214760.fastp.json
│   │   ├── 214763.fastp.html
│   │   ├── 214763.fastp.json
│   │   ├── 214765.fastp.html
│   │   ├── 214765.fastp.json
│   │   ├── 214766.fastp.html
│   │   ├── 214766.fastp.json
│   │   ├── 214780.fastp.html
│   │   ├── 214780.fastp.json
│   │   ├── 214807.fastp.html
│   │   ├── 214807.fastp.json
│   │   ├── fastqc
│   │   │   ├── 214704_1.trim_fastqc.html
│   │   │   ├── 214704_2.trim_fastqc.html
│   │   │   ├── 214721_1.trim_fastqc.html
│   │   │   ├── 214721_2.trim_fastqc.html
│   │   │   ├── 214724_1.trim_fastqc.html
│   │   │   ├── 214724_2.trim_fastqc.html
│   │   │   ├── 214760_1.trim_fastqc.html
│   │   │   ├── 214760_2.trim_fastqc.html
│   │   │   ├── 214763_1.trim_fastqc.html
│   │   │   ├── 214763_2.trim_fastqc.html
│   │   │   ├── 214765_1.trim_fastqc.html
│   │   │   ├── 214765_2.trim_fastqc.html
│   │   │   ├── 214766_1.trim_fastqc.html
│   │   │   ├── 214766_2.trim_fastqc.html
│   │   │   ├── 214780_1.trim_fastqc.html
│   │   │   ├── 214780_2.trim_fastqc.html
│   │   │   ├── 214807_1.trim_fastqc.html
│   │   │   ├── 214807_2.trim_fastqc.html
│   │   │   └── zips
│   │   │       ├── 214704_1.trim_fastqc.zip
│   │   │       ├── 214704_2.trim_fastqc.zip
│   │   │       ├── 214721_1.trim_fastqc.zip
│   │   │       ├── 214721_2.trim_fastqc.zip
│   │   │       ├── 214724_1.trim_fastqc.zip
│   │   │       ├── 214724_2.trim_fastqc.zip
│   │   │       ├── 214760_1.trim_fastqc.zip
│   │   │       ├── 214760_2.trim_fastqc.zip
│   │   │       ├── 214763_1.trim_fastqc.zip
│   │   │       ├── 214763_2.trim_fastqc.zip
│   │   │       ├── 214765_1.trim_fastqc.zip
│   │   │       ├── 214765_2.trim_fastqc.zip
│   │   │       ├── 214766_1.trim_fastqc.zip
│   │   │       ├── 214766_2.trim_fastqc.zip
│   │   │       ├── 214780_1.trim_fastqc.zip
│   │   │       ├── 214780_2.trim_fastqc.zip
│   │   │       ├── 214807_1.trim_fastqc.zip
│   │   │       └── 214807_2.trim_fastqc.zip
│   │   └── log
│   │       ├── 214704.fastp.log
│   │       ├── 214721.fastp.log
│   │       ├── 214724.fastp.log
│   │       ├── 214760.fastp.log
│   │       ├── 214763.fastp.log
│   │       ├── 214765.fastp.log
│   │       ├── 214766.fastp.log
│   │       ├── 214780.fastp.log
│   │       └── 214807.fastp.log
│   └── fastqc
│       ├── 214704_1.merged_fastqc.html
│       ├── 214704_2.merged_fastqc.html
│       ├── 214721_1.merged_fastqc.html
│       ├── 214721_2.merged_fastqc.html
│       ├── 214724_1.merged_fastqc.html
│       ├── 214724_2.merged_fastqc.html
│       ├── 214760_1.merged_fastqc.html
│       ├── 214760_2.merged_fastqc.html
│       ├── 214763_1.merged_fastqc.html
│       ├── 214763_2.merged_fastqc.html
│       ├── 214765_1.merged_fastqc.html
│       ├── 214765_2.merged_fastqc.html
│       ├── 214766_1.merged_fastqc.html
│       ├── 214766_2.merged_fastqc.html
│       ├── 214780_1.merged_fastqc.html
│       ├── 214780_2.merged_fastqc.html
│       ├── 214807_1.merged_fastqc.html
│       ├── 214807_2.merged_fastqc.html
│       └── zips
│           ├── 214704_1.merged_fastqc.zip
│           ├── 214704_2.merged_fastqc.zip
│           ├── 214721_1.merged_fastqc.zip
│           ├── 214721_2.merged_fastqc.zip
│           ├── 214724_1.merged_fastqc.zip
│           ├── 214724_2.merged_fastqc.zip
│           ├── 214760_1.merged_fastqc.zip
│           ├── 214760_2.merged_fastqc.zip
│           ├── 214763_1.merged_fastqc.zip
│           ├── 214763_2.merged_fastqc.zip
│           ├── 214765_1.merged_fastqc.zip
│           ├── 214765_2.merged_fastqc.zip
│           ├── 214766_1.merged_fastqc.zip
│           ├── 214766_2.merged_fastqc.zip
│           ├── 214780_1.merged_fastqc.zip
│           ├── 214780_2.merged_fastqc.zip
│           ├── 214807_1.merged_fastqc.zip
│           └── 214807_2.merged_fastqc.zip
└── variants
    ├── bam
    │   ├── 214704.bam
    │   ├── 214704.sorted.bam
    │   ├── 214704.sorted.bam.bai
    │   ├── 214704.trim.sorted.bam
    │   ├── 214704.trim.sorted.bam.bai
    │   ├── 214721.bam
    │   ├── 214721.sorted.bam
    │   ├── 214721.sorted.bam.bai
    │   ├── 214721.trim.sorted.bam
    │   ├── 214721.trim.sorted.bam.bai
    │   ├── 214724.bam
    │   ├── 214724.sorted.bam
    │   ├── 214724.sorted.bam.bai
    │   ├── 214724.trim.sorted.bam
    │   ├── 214724.trim.sorted.bam.bai
    │   ├── 214760.bam
    │   ├── 214760.sorted.bam
    │   ├── 214760.sorted.bam.bai
    │   ├── 214760.trim.sorted.bam
    │   ├── 214760.trim.sorted.bam.bai
    │   ├── 214763.bam
    │   ├── 214763.sorted.bam
    │   ├── 214763.sorted.bam.bai
    │   ├── 214763.trim.sorted.bam
    │   ├── 214763.trim.sorted.bam.bai
    │   ├── 214765.bam
    │   ├── 214765.sorted.bam
    │   ├── 214765.sorted.bam.bai
    │   ├── 214765.trim.sorted.bam
    │   ├── 214765.trim.sorted.bam.bai
    │   ├── 214766.bam
    │   ├── 214766.sorted.bam
    │   ├── 214766.sorted.bam.bai
    │   ├── 214766.trim.sorted.bam
    │   ├── 214766.trim.sorted.bam.bai
    │   ├── 214780.bam
    │   ├── 214780.sorted.bam
    │   ├── 214780.sorted.bam.bai
    │   ├── 214780.trim.sorted.bam
    │   ├── 214780.trim.sorted.bam.bai
    │   ├── 214807.bam
    │   ├── 214807.sorted.bam
    │   ├── 214807.sorted.bam.bai
    │   ├── 214807.trim.sorted.bam
    │   ├── 214807.trim.sorted.bam.bai
    │   ├── log
    │   │   ├── 214704.bowtie2.log
    │   │   ├── 214704.trim.ivar.log
    │   │   ├── 214721.bowtie2.log
    │   │   ├── 214721.trim.ivar.log
    │   │   ├── 214724.bowtie2.log
    │   │   ├── 214724.trim.ivar.log
    │   │   ├── 214760.bowtie2.log
    │   │   ├── 214760.trim.ivar.log
    │   │   ├── 214763.bowtie2.log
    │   │   ├── 214763.trim.ivar.log
    │   │   ├── 214765.bowtie2.log
    │   │   ├── 214765.trim.ivar.log
    │   │   ├── 214766.bowtie2.log
    │   │   ├── 214766.trim.ivar.log
    │   │   ├── 214780.bowtie2.log
    │   │   ├── 214780.trim.ivar.log
    │   │   ├── 214807.bowtie2.log
    │   │   └── 214807.trim.ivar.log
    │   ├── mosdepth
    │   │   ├── amplicon
    │   │   │   ├── 214704.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214704.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214704.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214704.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214704.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214704.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214704.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214704.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214704.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214721.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214721.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214721.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214721.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214721.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214721.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214721.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214721.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214721.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214724.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214724.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214724.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214724.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214724.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214724.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214724.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214724.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214724.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214760.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214760.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214760.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214760.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214760.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214760.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214760.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214760.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214760.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214763.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214763.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214763.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214763.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214763.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214763.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214763.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214763.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214763.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214765.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214765.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214765.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214765.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214765.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214765.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214765.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214765.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214765.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214766.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214766.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214766.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214766.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214766.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214766.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214766.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214766.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214766.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214780.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214780.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214780.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214780.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214780.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214780.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214780.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214780.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214780.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214807.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214807.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214807.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214807.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214807.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214807.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214807.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214807.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214807.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   └── plots
    │   │   │       ├── 214704.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214704.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214721.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214721.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214724.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214724.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214760.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214760.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214763.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214763.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214765.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214765.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214766.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214766.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214780.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214780.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214807.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214807.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── all_samples.trim.amplicon.regions.coverage.tsv
    │   │   │       └── all_samples.trim.amplicon.regions.heatmap.pdf
    │   │   └── genome
    │   │       ├── 214704.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214704.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214704.trim.genome.mosdepth.summary.txt
    │   │       ├── 214704.trim.genome.per-base.bed.gz
    │   │       ├── 214704.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214704.trim.genome.regions.bed.gz
    │   │       ├── 214704.trim.genome.regions.bed.gz.csi
    │   │       ├── 214721.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214721.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214721.trim.genome.mosdepth.summary.txt
    │   │       ├── 214721.trim.genome.per-base.bed.gz
    │   │       ├── 214721.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214721.trim.genome.regions.bed.gz
    │   │       ├── 214721.trim.genome.regions.bed.gz.csi
    │   │       ├── 214724.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214724.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214724.trim.genome.mosdepth.summary.txt
    │   │       ├── 214724.trim.genome.per-base.bed.gz
    │   │       ├── 214724.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214724.trim.genome.regions.bed.gz
    │   │       ├── 214724.trim.genome.regions.bed.gz.csi
    │   │       ├── 214760.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214760.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214760.trim.genome.mosdepth.summary.txt
    │   │       ├── 214760.trim.genome.per-base.bed.gz
    │   │       ├── 214760.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214760.trim.genome.regions.bed.gz
    │   │       ├── 214760.trim.genome.regions.bed.gz.csi
    │   │       ├── 214763.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214763.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214763.trim.genome.mosdepth.summary.txt
    │   │       ├── 214763.trim.genome.per-base.bed.gz
    │   │       ├── 214763.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214763.trim.genome.regions.bed.gz
    │   │       ├── 214763.trim.genome.regions.bed.gz.csi
    │   │       ├── 214765.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214765.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214765.trim.genome.mosdepth.summary.txt
    │   │       ├── 214765.trim.genome.per-base.bed.gz
    │   │       ├── 214765.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214765.trim.genome.regions.bed.gz
    │   │       ├── 214765.trim.genome.regions.bed.gz.csi
    │   │       ├── 214766.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214766.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214766.trim.genome.mosdepth.summary.txt
    │   │       ├── 214766.trim.genome.per-base.bed.gz
    │   │       ├── 214766.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214766.trim.genome.regions.bed.gz
    │   │       ├── 214766.trim.genome.regions.bed.gz.csi
    │   │       ├── 214780.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214780.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214780.trim.genome.mosdepth.summary.txt
    │   │       ├── 214780.trim.genome.per-base.bed.gz
    │   │       ├── 214780.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214780.trim.genome.regions.bed.gz
    │   │       ├── 214780.trim.genome.regions.bed.gz.csi
    │   │       ├── 214807.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214807.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214807.trim.genome.mosdepth.summary.txt
    │   │       ├── 214807.trim.genome.per-base.bed.gz
    │   │       ├── 214807.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214807.trim.genome.regions.bed.gz
    │   │       ├── 214807.trim.genome.regions.bed.gz.csi
    │   │       └── plots
    │   │           ├── 214704.trim.genome.regions.coverage.pdf
    │   │           ├── 214704.trim.genome.regions.coverage.tsv
    │   │           ├── 214721.trim.genome.regions.coverage.pdf
    │   │           ├── 214721.trim.genome.regions.coverage.tsv
    │   │           ├── 214724.trim.genome.regions.coverage.pdf
    │   │           ├── 214724.trim.genome.regions.coverage.tsv
    │   │           ├── 214760.trim.genome.regions.coverage.pdf
    │   │           ├── 214760.trim.genome.regions.coverage.tsv
    │   │           ├── 214763.trim.genome.regions.coverage.pdf
    │   │           ├── 214763.trim.genome.regions.coverage.tsv
    │   │           ├── 214765.trim.genome.regions.coverage.pdf
    │   │           ├── 214765.trim.genome.regions.coverage.tsv
    │   │           ├── 214766.trim.genome.regions.coverage.pdf
    │   │           ├── 214766.trim.genome.regions.coverage.tsv
    │   │           ├── 214780.trim.genome.regions.coverage.pdf
    │   │           ├── 214780.trim.genome.regions.coverage.tsv
    │   │           ├── 214807.trim.genome.regions.coverage.pdf
    │   │           ├── 214807.trim.genome.regions.coverage.tsv
    │   │           └── all_samples.trim.genome.regions.coverage.tsv
    │   ├── mpileup
    │   │   ├── 214704.trim.mpileup
    │   │   ├── 214721.trim.mpileup
    │   │   ├── 214724.trim.mpileup
    │   │   ├── 214760.trim.mpileup
    │   │   ├── 214763.trim.mpileup
    │   │   ├── 214765.trim.mpileup
    │   │   ├── 214766.trim.mpileup
    │   │   ├── 214780.trim.mpileup
    │   │   └── 214807.trim.mpileup
    │   ├── picard_metrics
    │   │   ├── 214704.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214704.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214704.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214704.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214704.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214704.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214704.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214704.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214704.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214704.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214721.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214721.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214721.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214721.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214721.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214721.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214721.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214721.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214721.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214721.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214724.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214724.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214724.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214724.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214724.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214724.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214724.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214724.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214724.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214724.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214760.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214760.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214760.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214760.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214760.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214760.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214760.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214760.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214760.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214760.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214763.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214763.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214763.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214763.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214763.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214763.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214763.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214763.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214763.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214763.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214765.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214765.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214765.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214765.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214765.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214765.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214765.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214765.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214765.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214765.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214766.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214766.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214766.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214766.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214766.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214766.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214766.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214766.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214766.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214766.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214780.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214780.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214780.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214780.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214780.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214780.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214780.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214780.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214780.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214780.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214807.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214807.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214807.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214807.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214807.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214807.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214807.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214807.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214807.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   └── 214807.trim.CollectWgsMetrics.coverage_metrics
    │   └── samtools_stats
    │       ├── 214704.sorted.bam.flagstat
    │       ├── 214704.sorted.bam.idxstats
    │       ├── 214704.sorted.bam.stats
    │       ├── 214704.trim.sorted.bam.flagstat
    │       ├── 214704.trim.sorted.bam.idxstats
    │       ├── 214704.trim.sorted.bam.stats
    │       ├── 214721.sorted.bam.flagstat
    │       ├── 214721.sorted.bam.idxstats
    │       ├── 214721.sorted.bam.stats
    │       ├── 214721.trim.sorted.bam.flagstat
    │       ├── 214721.trim.sorted.bam.idxstats
    │       ├── 214721.trim.sorted.bam.stats
    │       ├── 214724.sorted.bam.flagstat
    │       ├── 214724.sorted.bam.idxstats
    │       ├── 214724.sorted.bam.stats
    │       ├── 214724.trim.sorted.bam.flagstat
    │       ├── 214724.trim.sorted.bam.idxstats
    │       ├── 214724.trim.sorted.bam.stats
    │       ├── 214760.sorted.bam.flagstat
    │       ├── 214760.sorted.bam.idxstats
    │       ├── 214760.sorted.bam.stats
    │       ├── 214760.trim.sorted.bam.flagstat
    │       ├── 214760.trim.sorted.bam.idxstats
    │       ├── 214760.trim.sorted.bam.stats
    │       ├── 214763.sorted.bam.flagstat
    │       ├── 214763.sorted.bam.idxstats
    │       ├── 214763.sorted.bam.stats
    │       ├── 214763.trim.sorted.bam.flagstat
    │       ├── 214763.trim.sorted.bam.idxstats
    │       ├── 214763.trim.sorted.bam.stats
    │       ├── 214765.sorted.bam.flagstat
    │       ├── 214765.sorted.bam.idxstats
    │       ├── 214765.sorted.bam.stats
    │       ├── 214765.trim.sorted.bam.flagstat
    │       ├── 214765.trim.sorted.bam.idxstats
    │       ├── 214765.trim.sorted.bam.stats
    │       ├── 214766.sorted.bam.flagstat
    │       ├── 214766.sorted.bam.idxstats
    │       ├── 214766.sorted.bam.stats
    │       ├── 214766.trim.sorted.bam.flagstat
    │       ├── 214766.trim.sorted.bam.idxstats
    │       ├── 214766.trim.sorted.bam.stats
    │       ├── 214780.sorted.bam.flagstat
    │       ├── 214780.sorted.bam.idxstats
    │       ├── 214780.sorted.bam.stats
    │       ├── 214780.trim.sorted.bam.flagstat
    │       ├── 214780.trim.sorted.bam.idxstats
    │       ├── 214780.trim.sorted.bam.stats
    │       ├── 214807.sorted.bam.flagstat
    │       ├── 214807.sorted.bam.idxstats
    │       ├── 214807.sorted.bam.stats
    │       ├── 214807.trim.sorted.bam.flagstat
    │       ├── 214807.trim.sorted.bam.idxstats
    │       └── 214807.trim.sorted.bam.stats
    ├── bcftools
    │   ├── 214704.vcf.gz
    │   ├── 214704.vcf.gz.tbi
    │   ├── 214721.vcf.gz
    │   ├── 214721.vcf.gz.tbi
    │   ├── 214724.vcf.gz
    │   ├── 214724.vcf.gz.tbi
    │   ├── 214760.vcf.gz
    │   ├── 214760.vcf.gz.tbi
    │   ├── 214763.vcf.gz
    │   ├── 214763.vcf.gz.tbi
    │   ├── 214765.vcf.gz
    │   ├── 214765.vcf.gz.tbi
    │   ├── 214766.vcf.gz
    │   ├── 214766.vcf.gz.tbi
    │   ├── 214780.vcf.gz
    │   ├── 214780.vcf.gz.tbi
    │   ├── 214807.vcf.gz
    │   ├── 214807.vcf.gz.tbi
    │   ├── bcftools_stats
    │   │   ├── 214704.bcftools_stats.txt
    │   │   ├── 214721.bcftools_stats.txt
    │   │   ├── 214724.bcftools_stats.txt
    │   │   ├── 214760.bcftools_stats.txt
    │   │   ├── 214763.bcftools_stats.txt
    │   │   ├── 214765.bcftools_stats.txt
    │   │   ├── 214766.bcftools_stats.txt
    │   │   ├── 214780.bcftools_stats.txt
    │   │   └── 214807.bcftools_stats.txt
    │   ├── consensus
    │   │   ├── 214704.consensus.masked.fa
    │   │   ├── 214721.consensus.masked.fa
    │   │   ├── 214724.consensus.masked.fa
    │   │   ├── 214760.consensus.masked.fa
    │   │   ├── 214763.consensus.masked.fa
    │   │   ├── 214765.consensus.masked.fa
    │   │   ├── 214766.consensus.masked.fa
    │   │   ├── 214780.consensus.masked.fa
    │   │   ├── 214807.consensus.masked.fa
    │   │   └── base_qc
    │   │       ├── 214704.ACTG_density.pdf
    │   │       ├── 214704.base_counts.pdf
    │   │       ├── 214704.base_counts.tsv
    │   │       ├── 214704.N_density.pdf
    │   │       ├── 214704.N_run.tsv
    │   │       ├── 214721.ACTG_density.pdf
    │   │       ├── 214721.base_counts.pdf
    │   │       ├── 214721.base_counts.tsv
    │   │       ├── 214721.N_density.pdf
    │   │       ├── 214721.N_run.tsv
    │   │       ├── 214724.ACTG_density.pdf
    │   │       ├── 214724.base_counts.pdf
    │   │       ├── 214724.base_counts.tsv
    │   │       ├── 214724.N_density.pdf
    │   │       ├── 214724.N_run.tsv
    │   │       ├── 214760.ACTG_density.pdf
    │   │       ├── 214760.base_counts.pdf
    │   │       ├── 214760.base_counts.tsv
    │   │       ├── 214760.N_density.pdf
    │   │       ├── 214760.N_run.tsv
    │   │       ├── 214763.ACTG_density.pdf
    │   │       ├── 214763.base_counts.pdf
    │   │       ├── 214763.base_counts.tsv
    │   │       ├── 214763.N_density.pdf
    │   │       ├── 214763.N_run.tsv
    │   │       ├── 214765.ACTG_density.pdf
    │   │       ├── 214765.base_counts.pdf
    │   │       ├── 214765.base_counts.tsv
    │   │       ├── 214765.N_density.pdf
    │   │       ├── 214765.N_run.tsv
    │   │       ├── 214766.ACTG_density.pdf
    │   │       ├── 214766.base_counts.pdf
    │   │       ├── 214766.base_counts.tsv
    │   │       ├── 214766.N_density.pdf
    │   │       ├── 214766.N_run.tsv
    │   │       ├── 214780.ACTG_density.pdf
    │   │       ├── 214780.base_counts.pdf
    │   │       ├── 214780.base_counts.tsv
    │   │       ├── 214780.N_density.pdf
    │   │       ├── 214780.N_run.tsv
    │   │       ├── 214807.ACTG_density.pdf
    │   │       ├── 214807.base_counts.pdf
    │   │       ├── 214807.base_counts.tsv
    │   │       ├── 214807.N_density.pdf
    │   │       └── 214807.N_run.tsv
    │   ├── quast
    │   │   ├── aligned_stats
    │   │   │   ├── cumulative_plot.pdf
    │   │   │   ├── NAx_plot.pdf
    │   │   │   └── NGAx_plot.pdf
    │   │   ├── basic_stats
    │   │   │   ├── 214704.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214721.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214724.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214760.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214763.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214765.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214766.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214780.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214807.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── cumulative_plot.pdf
    │   │   │   ├── GC_content_plot.pdf
    │   │   │   ├── gc.icarus.txt
    │   │   │   ├── NGx_plot.pdf
    │   │   │   └── Nx_plot.pdf
    │   │   ├── contigs_reports
    │   │   │   ├── 214704_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214721_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214724_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214760_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214763_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214765_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214766_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214780_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214807_consensus_masked.mis_contigs.fa
    │   │   │   ├── all_alignments_214704-consensus-masked.tsv
    │   │   │   ├── all_alignments_214721-consensus-masked.tsv
    │   │   │   ├── all_alignments_214724-consensus-masked.tsv
    │   │   │   ├── all_alignments_214760-consensus-masked.tsv
    │   │   │   ├── all_alignments_214763-consensus-masked.tsv
    │   │   │   ├── all_alignments_214765-consensus-masked.tsv
    │   │   │   ├── all_alignments_214766-consensus-masked.tsv
    │   │   │   ├── all_alignments_214780-consensus-masked.tsv
    │   │   │   ├── all_alignments_214807-consensus-masked.tsv
    │   │   │   ├── contigs_report_214704-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214704-consensus-masked.stderr
    │   │   │   ├── contigs_report_214704-consensus-masked.stdout
    │   │   │   ├── contigs_report_214704-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214721-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214721-consensus-masked.stderr
    │   │   │   ├── contigs_report_214721-consensus-masked.stdout
    │   │   │   ├── contigs_report_214721-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214724-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214724-consensus-masked.stderr
    │   │   │   ├── contigs_report_214724-consensus-masked.stdout
    │   │   │   ├── contigs_report_214724-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214760-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214760-consensus-masked.stderr
    │   │   │   ├── contigs_report_214760-consensus-masked.stdout
    │   │   │   ├── contigs_report_214760-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214763-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214763-consensus-masked.stderr
    │   │   │   ├── contigs_report_214763-consensus-masked.stdout
    │   │   │   ├── contigs_report_214763-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214765-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214765-consensus-masked.stderr
    │   │   │   ├── contigs_report_214765-consensus-masked.stdout
    │   │   │   ├── contigs_report_214765-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214766-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214766-consensus-masked.stderr
    │   │   │   ├── contigs_report_214766-consensus-masked.stdout
    │   │   │   ├── contigs_report_214766-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214780-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214780-consensus-masked.stderr
    │   │   │   ├── contigs_report_214780-consensus-masked.stdout
    │   │   │   ├── contigs_report_214780-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214807-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214807-consensus-masked.stderr
    │   │   │   ├── contigs_report_214807-consensus-masked.stdout
    │   │   │   ├── contigs_report_214807-consensus-masked.unaligned.info
    │   │   │   ├── minimap_output
    │   │   │   │   ├── 214704-consensus-masked.coords
    │   │   │   │   ├── 214704-consensus-masked.coords.filtered
    │   │   │   │   ├── 214704-consensus-masked.coords_tmp
    │   │   │   │   ├── 214704-consensus-masked.sf
    │   │   │   │   ├── 214704-consensus-masked.unaligned
    │   │   │   │   ├── 214704-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214721-consensus-masked.coords
    │   │   │   │   ├── 214721-consensus-masked.coords.filtered
    │   │   │   │   ├── 214721-consensus-masked.coords_tmp
    │   │   │   │   ├── 214721-consensus-masked.sf
    │   │   │   │   ├── 214721-consensus-masked.unaligned
    │   │   │   │   ├── 214721-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214724-consensus-masked.coords
    │   │   │   │   ├── 214724-consensus-masked.coords.filtered
    │   │   │   │   ├── 214724-consensus-masked.coords_tmp
    │   │   │   │   ├── 214724-consensus-masked.sf
    │   │   │   │   ├── 214724-consensus-masked.unaligned
    │   │   │   │   ├── 214724-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214760-consensus-masked.coords
    │   │   │   │   ├── 214760-consensus-masked.coords.filtered
    │   │   │   │   ├── 214760-consensus-masked.coords_tmp
    │   │   │   │   ├── 214760-consensus-masked.sf
    │   │   │   │   ├── 214760-consensus-masked.unaligned
    │   │   │   │   ├── 214760-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214763-consensus-masked.coords
    │   │   │   │   ├── 214763-consensus-masked.coords.filtered
    │   │   │   │   ├── 214763-consensus-masked.coords_tmp
    │   │   │   │   ├── 214763-consensus-masked.sf
    │   │   │   │   ├── 214763-consensus-masked.unaligned
    │   │   │   │   ├── 214763-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214765-consensus-masked.coords
    │   │   │   │   ├── 214765-consensus-masked.coords.filtered
    │   │   │   │   ├── 214765-consensus-masked.coords_tmp
    │   │   │   │   ├── 214765-consensus-masked.sf
    │   │   │   │   ├── 214765-consensus-masked.unaligned
    │   │   │   │   ├── 214765-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214766-consensus-masked.coords
    │   │   │   │   ├── 214766-consensus-masked.coords.filtered
    │   │   │   │   ├── 214766-consensus-masked.coords_tmp
    │   │   │   │   ├── 214766-consensus-masked.sf
    │   │   │   │   ├── 214766-consensus-masked.unaligned
    │   │   │   │   ├── 214766-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214780-consensus-masked.coords
    │   │   │   │   ├── 214780-consensus-masked.coords.filtered
    │   │   │   │   ├── 214780-consensus-masked.coords_tmp
    │   │   │   │   ├── 214780-consensus-masked.sf
    │   │   │   │   ├── 214780-consensus-masked.unaligned
    │   │   │   │   ├── 214780-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214807-consensus-masked.coords
    │   │   │   │   ├── 214807-consensus-masked.coords.filtered
    │   │   │   │   ├── 214807-consensus-masked.coords_tmp
    │   │   │   │   ├── 214807-consensus-masked.sf
    │   │   │   │   ├── 214807-consensus-masked.unaligned
    │   │   │   │   └── 214807-consensus-masked.used_snps.gz
    │   │   │   ├── misassemblies_frcurve_plot.pdf
    │   │   │   ├── misassemblies_plot.pdf
    │   │   │   ├── misassemblies_report.tex
    │   │   │   ├── misassemblies_report.tsv
    │   │   │   ├── misassemblies_report.txt
    │   │   │   ├── transposed_report_misassemblies.tex
    │   │   │   ├── transposed_report_misassemblies.tsv
    │   │   │   ├── transposed_report_misassemblies.txt
    │   │   │   ├── unaligned_report.tex
    │   │   │   ├── unaligned_report.tsv
    │   │   │   └── unaligned_report.txt
    │   │   ├── genome_stats
    │   │   │   ├── 214704-consensus-masked_gaps.txt
    │   │   │   ├── 214704-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214721-consensus-masked_gaps.txt
    │   │   │   ├── 214721-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214724-consensus-masked_gaps.txt
    │   │   │   ├── 214724-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214760-consensus-masked_gaps.txt
    │   │   │   ├── 214760-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214763-consensus-masked_gaps.txt
    │   │   │   ├── 214763-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214765-consensus-masked_gaps.txt
    │   │   │   ├── 214765-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214766-consensus-masked_gaps.txt
    │   │   │   ├── 214766-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214780-consensus-masked_gaps.txt
    │   │   │   ├── 214780-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214807-consensus-masked_gaps.txt
    │   │   │   ├── 214807-consensus-masked_genomic_features_any.txt
    │   │   │   ├── complete_features_histogram.pdf
    │   │   │   ├── features_cumulative_plot.pdf
    │   │   │   ├── features_frcurve_plot.pdf
    │   │   │   ├── genome_fraction_histogram.pdf
    │   │   │   └── genome_info.txt
    │   │   ├── icarus.html
    │   │   ├── icarus_viewers
    │   │   │   ├── alignment_viewer.html
    │   │   │   └── contig_size_viewer.html
    │   │   ├── quast.log
    │   │   ├── report.html
    │   │   ├── report.pdf
    │   │   ├── report.tex
    │   │   ├── report.tsv
    │   │   ├── report.txt
    │   │   ├── transposed_report.tex
    │   │   ├── transposed_report.tsv
    │   │   └── transposed_report.txt
    │   └── snpeff
    │       ├── 214704.snpEff.csv
    │       ├── 214704.snpEff.genes.txt
    │       ├── 214704.snpEff.summary.html
    │       ├── 214704.snpEff.vcf.gz
    │       ├── 214704.snpEff.vcf.gz.tbi
    │       ├── 214704.snpSift.table.txt
    │       ├── 214721.snpEff.csv
    │       ├── 214721.snpEff.genes.txt
    │       ├── 214721.snpEff.summary.html
    │       ├── 214721.snpEff.vcf.gz
    │       ├── 214721.snpEff.vcf.gz.tbi
    │       ├── 214721.snpSift.table.txt
    │       ├── 214724.snpEff.csv
    │       ├── 214724.snpEff.genes.txt
    │       ├── 214724.snpEff.summary.html
    │       ├── 214724.snpEff.vcf.gz
    │       ├── 214724.snpEff.vcf.gz.tbi
    │       ├── 214724.snpSift.table.txt
    │       ├── 214760.snpEff.csv
    │       ├── 214760.snpEff.genes.txt
    │       ├── 214760.snpEff.summary.html
    │       ├── 214760.snpEff.vcf.gz
    │       ├── 214760.snpEff.vcf.gz.tbi
    │       ├── 214760.snpSift.table.txt
    │       ├── 214763.snpEff.csv
    │       ├── 214763.snpEff.genes.txt
    │       ├── 214763.snpEff.summary.html
    │       ├── 214763.snpEff.vcf.gz
    │       ├── 214763.snpEff.vcf.gz.tbi
    │       ├── 214763.snpSift.table.txt
    │       ├── 214765.snpEff.csv
    │       ├── 214765.snpEff.genes.txt
    │       ├── 214765.snpEff.summary.html
    │       ├── 214765.snpEff.vcf.gz
    │       ├── 214765.snpEff.vcf.gz.tbi
    │       ├── 214765.snpSift.table.txt
    │       ├── 214766.snpEff.csv
    │       ├── 214766.snpEff.genes.txt
    │       ├── 214766.snpEff.summary.html
    │       ├── 214766.snpEff.vcf.gz
    │       ├── 214766.snpEff.vcf.gz.tbi
    │       ├── 214766.snpSift.table.txt
    │       ├── 214780.snpEff.csv
    │       ├── 214780.snpEff.genes.txt
    │       ├── 214780.snpEff.summary.html
    │       ├── 214780.snpEff.vcf.gz
    │       ├── 214780.snpEff.vcf.gz.tbi
    │       ├── 214780.snpSift.table.txt
    │       ├── 214807.snpEff.csv
    │       ├── 214807.snpEff.genes.txt
    │       ├── 214807.snpEff.summary.html
    │       ├── 214807.snpEff.vcf.gz
    │       ├── 214807.snpEff.vcf.gz.tbi
    │       └── 214807.snpSift.table.txt
    ├── intersect
    │   ├── 214704
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214721
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214724
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214760
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214763
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214765
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214766
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214780
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   └── 214807
    │       ├── 0000.vcf.gz
    │       ├── 0000.vcf.gz.tbi
    │       ├── 0001.vcf.gz
    │       ├── 0001.vcf.gz.tbi
    │       ├── 0002.vcf.gz
    │       ├── 0002.vcf.gz.tbi
    │       ├── README.txt
    │       └── sites.txt
    ├── ivar
    │   ├── 214704.AF0.75.vcf.gz
    │   ├── 214704.AF0.75.vcf.gz.tbi
    │   ├── 214704.tsv
    │   ├── 214704.vcf.gz
    │   ├── 214704.vcf.gz.tbi
    │   ├── 214721.AF0.75.vcf.gz
    │   ├── 214721.AF0.75.vcf.gz.tbi
    │   ├── 214721.tsv
    │   ├── 214721.vcf.gz
    │   ├── 214721.vcf.gz.tbi
    │   ├── 214724.AF0.75.vcf.gz
    │   ├── 214724.AF0.75.vcf.gz.tbi
    │   ├── 214724.tsv
    │   ├── 214724.vcf.gz
    │   ├── 214724.vcf.gz.tbi
    │   ├── 214760.AF0.75.vcf.gz
    │   ├── 214760.AF0.75.vcf.gz.tbi
    │   ├── 214760.tsv
    │   ├── 214760.vcf.gz
    │   ├── 214760.vcf.gz.tbi
    │   ├── 214763.AF0.75.vcf.gz
    │   ├── 214763.AF0.75.vcf.gz.tbi
    │   ├── 214763.tsv
    │   ├── 214763.vcf.gz
    │   ├── 214763.vcf.gz.tbi
    │   ├── 214765.AF0.75.vcf.gz
    │   ├── 214765.AF0.75.vcf.gz.tbi
    │   ├── 214765.tsv
    │   ├── 214765.vcf.gz
    │   ├── 214765.vcf.gz.tbi
    │   ├── 214766.AF0.75.vcf.gz
    │   ├── 214766.AF0.75.vcf.gz.tbi
    │   ├── 214766.tsv
    │   ├── 214766.vcf.gz
    │   ├── 214766.vcf.gz.tbi
    │   ├── 214780.AF0.75.vcf.gz
    │   ├── 214780.AF0.75.vcf.gz.tbi
    │   ├── 214780.tsv
    │   ├── 214780.vcf.gz
    │   ├── 214780.vcf.gz.tbi
    │   ├── 214807.AF0.75.vcf.gz
    │   ├── 214807.AF0.75.vcf.gz.tbi
    │   ├── 214807.tsv
    │   ├── 214807.vcf.gz
    │   ├── 214807.vcf.gz.tbi
    │   ├── bcftools_stats
    │   │   ├── 214704.AF0.75.bcftools_stats.txt
    │   │   ├── 214704.bcftools_stats.txt
    │   │   ├── 214721.AF0.75.bcftools_stats.txt
    │   │   ├── 214721.bcftools_stats.txt
    │   │   ├── 214724.AF0.75.bcftools_stats.txt
    │   │   ├── 214724.bcftools_stats.txt
    │   │   ├── 214760.AF0.75.bcftools_stats.txt
    │   │   ├── 214760.bcftools_stats.txt
    │   │   ├── 214763.AF0.75.bcftools_stats.txt
    │   │   ├── 214763.bcftools_stats.txt
    │   │   ├── 214765.AF0.75.bcftools_stats.txt
    │   │   ├── 214765.bcftools_stats.txt
    │   │   ├── 214766.AF0.75.bcftools_stats.txt
    │   │   ├── 214766.bcftools_stats.txt
    │   │   ├── 214780.AF0.75.bcftools_stats.txt
    │   │   ├── 214780.bcftools_stats.txt
    │   │   ├── 214807.AF0.75.bcftools_stats.txt
    │   │   └── 214807.bcftools_stats.txt
    │   ├── consensus
    │   │   ├── 214704.AF0.75.consensus.fa
    │   │   ├── 214704.AF0.75.consensus.qual.txt
    │   │   ├── 214721.AF0.75.consensus.fa
    │   │   ├── 214721.AF0.75.consensus.qual.txt
    │   │   ├── 214724.AF0.75.consensus.fa
    │   │   ├── 214724.AF0.75.consensus.qual.txt
    │   │   ├── 214760.AF0.75.consensus.fa
    │   │   ├── 214760.AF0.75.consensus.qual.txt
    │   │   ├── 214763.AF0.75.consensus.fa
    │   │   ├── 214763.AF0.75.consensus.qual.txt
    │   │   ├── 214765.AF0.75.consensus.fa
    │   │   ├── 214765.AF0.75.consensus.qual.txt
    │   │   ├── 214766.AF0.75.consensus.fa
    │   │   ├── 214766.AF0.75.consensus.qual.txt
    │   │   ├── 214780.AF0.75.consensus.fa
    │   │   ├── 214780.AF0.75.consensus.qual.txt
    │   │   ├── 214807.AF0.75.consensus.fa
    │   │   ├── 214807.AF0.75.consensus.qual.txt
    │   │   └── base_qc
    │   │       ├── 214704.AF0.75.ACTG_density.pdf
    │   │       ├── 214704.AF0.75.base_counts.pdf
    │   │       ├── 214704.AF0.75.base_counts.tsv
    │   │       ├── 214704.AF0.75.N_density.pdf
    │   │       ├── 214704.AF0.75.N_run.tsv
    │   │       ├── 214704.AF0.75.R_density.pdf
    │   │       ├── 214704.AF0.75.Y_density.pdf
    │   │       ├── 214721.AF0.75.ACTG_density.pdf
    │   │       ├── 214721.AF0.75.base_counts.pdf
    │   │       ├── 214721.AF0.75.base_counts.tsv
    │   │       ├── 214721.AF0.75.N_density.pdf
    │   │       ├── 214721.AF0.75.N_run.tsv
    │   │       ├── 214721.AF0.75.W_density.pdf
    │   │       ├── 214721.AF0.75.Y_density.pdf
    │   │       ├── 214724.AF0.75.ACTG_density.pdf
    │   │       ├── 214724.AF0.75.base_counts.pdf
    │   │       ├── 214724.AF0.75.base_counts.tsv
    │   │       ├── 214724.AF0.75.M_density.pdf
    │   │       ├── 214724.AF0.75.N_density.pdf
    │   │       ├── 214724.AF0.75.N_run.tsv
    │   │       ├── 214724.AF0.75.R_density.pdf
    │   │       ├── 214760.AF0.75.ACTG_density.pdf
    │   │       ├── 214760.AF0.75.base_counts.pdf
    │   │       ├── 214760.AF0.75.base_counts.tsv
    │   │       ├── 214760.AF0.75.N_density.pdf
    │   │       ├── 214760.AF0.75.N_run.tsv
    │   │       ├── 214760.AF0.75.R_density.pdf
    │   │       ├── 214760.AF0.75.W_density.pdf
    │   │       ├── 214760.AF0.75.Y_density.pdf
    │   │       ├── 214763.AF0.75.ACTG_density.pdf
    │   │       ├── 214763.AF0.75.base_counts.pdf
    │   │       ├── 214763.AF0.75.base_counts.tsv
    │   │       ├── 214763.AF0.75.N_density.pdf
    │   │       ├── 214763.AF0.75.N_run.tsv
    │   │       ├── 214763.AF0.75.R_density.pdf
    │   │       ├── 214763.AF0.75.Y_density.pdf
    │   │       ├── 214765.AF0.75.ACTG_density.pdf
    │   │       ├── 214765.AF0.75.base_counts.pdf
    │   │       ├── 214765.AF0.75.base_counts.tsv
    │   │       ├── 214765.AF0.75.N_density.pdf
    │   │       ├── 214765.AF0.75.N_run.tsv
    │   │       ├── 214765.AF0.75.R_density.pdf
    │   │       ├── 214766.AF0.75.ACTG_density.pdf
    │   │       ├── 214766.AF0.75.base_counts.pdf
    │   │       ├── 214766.AF0.75.base_counts.tsv
    │   │       ├── 214766.AF0.75.N_density.pdf
    │   │       ├── 214766.AF0.75.N_run.tsv
    │   │       ├── 214766.AF0.75.R_density.pdf
    │   │       ├── 214780.AF0.75.ACTG_density.pdf
    │   │       ├── 214780.AF0.75.base_counts.pdf
    │   │       ├── 214780.AF0.75.base_counts.tsv
    │   │       ├── 214780.AF0.75.N_density.pdf
    │   │       ├── 214780.AF0.75.N_run.tsv
    │   │       ├── 214780.AF0.75.R_density.pdf
    │   │       ├── 214807.AF0.75.ACTG_density.pdf
    │   │       ├── 214807.AF0.75.base_counts.pdf
    │   │       ├── 214807.AF0.75.base_counts.tsv
    │   │       ├── 214807.AF0.75.M_density.pdf
    │   │       ├── 214807.AF0.75.N_density.pdf
    │   │       ├── 214807.AF0.75.N_run.tsv
    │   │       ├── 214807.AF0.75.R_density.pdf
    │   │       └── 214807.AF0.75.W_density.pdf
    │   ├── log
    │   │   ├── 214704.AF0.75.variant.counts.log
    │   │   ├── 214704.variant.counts.log
    │   │   ├── 214721.AF0.75.variant.counts.log
    │   │   ├── 214721.variant.counts.log
    │   │   ├── 214724.AF0.75.variant.counts.log
    │   │   ├── 214724.variant.counts.log
    │   │   ├── 214760.AF0.75.variant.counts.log
    │   │   ├── 214760.variant.counts.log
    │   │   ├── 214763.AF0.75.variant.counts.log
    │   │   ├── 214763.variant.counts.log
    │   │   ├── 214765.AF0.75.variant.counts.log
    │   │   ├── 214765.variant.counts.log
    │   │   ├── 214766.AF0.75.variant.counts.log
    │   │   ├── 214766.variant.counts.log
    │   │   ├── 214780.AF0.75.variant.counts.log
    │   │   ├── 214780.variant.counts.log
    │   │   ├── 214807.AF0.75.variant.counts.log
    │   │   └── 214807.variant.counts.log
    │   ├── quast
    │   │   └── AF0.75
    │   │       ├── aligned_stats
    │   │       │   ├── cumulative_plot.pdf
    │   │       │   ├── NAx_plot.pdf
    │   │       │   └── NGAx_plot.pdf
    │   │       ├── basic_stats
    │   │       │   ├── 214704.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214721.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214724.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214760.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214763.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214765.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214766.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214780.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214807.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── cumulative_plot.pdf
    │   │       │   ├── GC_content_plot.pdf
    │   │       │   ├── gc.icarus.txt
    │   │       │   ├── NGx_plot.pdf
    │   │       │   └── Nx_plot.pdf
    │   │       ├── contigs_reports
    │   │       │   ├── 214704_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214721_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214724_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214760_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214763_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214765_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214766_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214780_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214807_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── all_alignments_214704-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214721-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214724-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214760-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214763-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214765-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214766-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214780-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214807-AF0-75-consensus.tsv
    │   │       │   ├── contigs_report_214704-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214704-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214704-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214704-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214721-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214721-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214721-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214721-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214724-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214724-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214724-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214724-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214760-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214760-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214760-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214760-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214763-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214763-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214763-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214763-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214765-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214765-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214765-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214765-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214766-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214766-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214766-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214766-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214780-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214780-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214780-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214780-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214807-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214807-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214807-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214807-AF0-75-consensus.unaligned.info
    │   │       │   ├── minimap_output
    │   │       │   │   ├── 214704-AF0-75-consensus.coords
    │   │       │   │   ├── 214704-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214704-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214704-AF0-75-consensus.sf
    │   │       │   │   ├── 214704-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214704-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214721-AF0-75-consensus.coords
    │   │       │   │   ├── 214721-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214721-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214721-AF0-75-consensus.sf
    │   │       │   │   ├── 214721-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214721-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214724-AF0-75-consensus.coords
    │   │       │   │   ├── 214724-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214724-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214724-AF0-75-consensus.sf
    │   │       │   │   ├── 214724-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214724-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214760-AF0-75-consensus.coords
    │   │       │   │   ├── 214760-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214760-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214760-AF0-75-consensus.sf
    │   │       │   │   ├── 214760-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214760-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214763-AF0-75-consensus.coords
    │   │       │   │   ├── 214763-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214763-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214763-AF0-75-consensus.sf
    │   │       │   │   ├── 214763-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214763-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214765-AF0-75-consensus.coords
    │   │       │   │   ├── 214765-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214765-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214765-AF0-75-consensus.sf
    │   │       │   │   ├── 214765-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214765-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214766-AF0-75-consensus.coords
    │   │       │   │   ├── 214766-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214766-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214766-AF0-75-consensus.sf
    │   │       │   │   ├── 214766-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214766-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214780-AF0-75-consensus.coords
    │   │       │   │   ├── 214780-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214780-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214780-AF0-75-consensus.sf
    │   │       │   │   ├── 214780-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214780-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214807-AF0-75-consensus.coords
    │   │       │   │   ├── 214807-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214807-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214807-AF0-75-consensus.sf
    │   │       │   │   ├── 214807-AF0-75-consensus.unaligned
    │   │       │   │   └── 214807-AF0-75-consensus.used_snps.gz
    │   │       │   ├── misassemblies_frcurve_plot.pdf
    │   │       │   ├── misassemblies_plot.pdf
    │   │       │   ├── misassemblies_report.tex
    │   │       │   ├── misassemblies_report.tsv
    │   │       │   ├── misassemblies_report.txt
    │   │       │   ├── transposed_report_misassemblies.tex
    │   │       │   ├── transposed_report_misassemblies.tsv
    │   │       │   ├── transposed_report_misassemblies.txt
    │   │       │   ├── unaligned_report.tex
    │   │       │   ├── unaligned_report.tsv
    │   │       │   └── unaligned_report.txt
    │   │       ├── genome_stats
    │   │       │   ├── 214704-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214704-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214721-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214721-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214724-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214724-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214760-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214760-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214763-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214763-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214765-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214765-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214766-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214766-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214780-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214780-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214807-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214807-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── complete_features_histogram.pdf
    │   │       │   ├── features_cumulative_plot.pdf
    │   │       │   ├── features_frcurve_plot.pdf
    │   │       │   ├── genome_fraction_histogram.pdf
    │   │       │   └── genome_info.txt
    │   │       ├── icarus.html
    │   │       ├── icarus_viewers
    │   │       │   ├── alignment_viewer.html
    │   │       │   └── contig_size_viewer.html
    │   │       ├── quast.log
    │   │       ├── report.html
    │   │       ├── report.pdf
    │   │       ├── report.tex
    │   │       ├── report.tsv
    │   │       ├── report.txt
    │   │       ├── transposed_report.tex
    │   │       ├── transposed_report.tsv
    │   │       └── transposed_report.txt
    │   └── snpeff
    │       ├── 214704.AF0.75.snpEff.csv
    │       ├── 214704.AF0.75.snpEff.genes.txt
    │       ├── 214704.AF0.75.snpEff.summary.html
    │       ├── 214704.AF0.75.snpEff.vcf.gz
    │       ├── 214704.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214704.AF0.75.snpSift.table.txt
    │       ├── 214704.snpEff.csv
    │       ├── 214704.snpEff.genes.txt
    │       ├── 214704.snpEff.summary.html
    │       ├── 214704.snpEff.vcf.gz
    │       ├── 214704.snpEff.vcf.gz.tbi
    │       ├── 214704.snpSift.table.txt
    │       ├── 214721.AF0.75.snpEff.csv
    │       ├── 214721.AF0.75.snpEff.genes.txt
    │       ├── 214721.AF0.75.snpEff.summary.html
    │       ├── 214721.AF0.75.snpEff.vcf.gz
    │       ├── 214721.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214721.AF0.75.snpSift.table.txt
    │       ├── 214721.snpEff.csv
    │       ├── 214721.snpEff.genes.txt
    │       ├── 214721.snpEff.summary.html
    │       ├── 214721.snpEff.vcf.gz
    │       ├── 214721.snpEff.vcf.gz.tbi
    │       ├── 214721.snpSift.table.txt
    │       ├── 214724.AF0.75.snpEff.csv
    │       ├── 214724.AF0.75.snpEff.genes.txt
    │       ├── 214724.AF0.75.snpEff.summary.html
    │       ├── 214724.AF0.75.snpEff.vcf.gz
    │       ├── 214724.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214724.AF0.75.snpSift.table.txt
    │       ├── 214724.snpEff.csv
    │       ├── 214724.snpEff.genes.txt
    │       ├── 214724.snpEff.summary.html
    │       ├── 214724.snpEff.vcf.gz
    │       ├── 214724.snpEff.vcf.gz.tbi
    │       ├── 214724.snpSift.table.txt
    │       ├── 214760.AF0.75.snpEff.csv
    │       ├── 214760.AF0.75.snpEff.genes.txt
    │       ├── 214760.AF0.75.snpEff.summary.html
    │       ├── 214760.AF0.75.snpEff.vcf.gz
    │       ├── 214760.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214760.AF0.75.snpSift.table.txt
    │       ├── 214760.snpEff.csv
    │       ├── 214760.snpEff.genes.txt
    │       ├── 214760.snpEff.summary.html
    │       ├── 214760.snpEff.vcf.gz
    │       ├── 214760.snpEff.vcf.gz.tbi
    │       ├── 214760.snpSift.table.txt
    │       ├── 214763.AF0.75.snpEff.csv
    │       ├── 214763.AF0.75.snpEff.genes.txt
    │       ├── 214763.AF0.75.snpEff.summary.html
    │       ├── 214763.AF0.75.snpEff.vcf.gz
    │       ├── 214763.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214763.AF0.75.snpSift.table.txt
    │       ├── 214763.snpEff.csv
    │       ├── 214763.snpEff.genes.txt
    │       ├── 214763.snpEff.summary.html
    │       ├── 214763.snpEff.vcf.gz
    │       ├── 214763.snpEff.vcf.gz.tbi
    │       ├── 214763.snpSift.table.txt
    │       ├── 214765.AF0.75.snpEff.csv
    │       ├── 214765.AF0.75.snpEff.genes.txt
    │       ├── 214765.AF0.75.snpEff.summary.html
    │       ├── 214765.AF0.75.snpEff.vcf.gz
    │       ├── 214765.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214765.AF0.75.snpSift.table.txt
    │       ├── 214765.snpEff.csv
    │       ├── 214765.snpEff.genes.txt
    │       ├── 214765.snpEff.summary.html
    │       ├── 214765.snpEff.vcf.gz
    │       ├── 214765.snpEff.vcf.gz.tbi
    │       ├── 214765.snpSift.table.txt
    │       ├── 214766.AF0.75.snpEff.csv
    │       ├── 214766.AF0.75.snpEff.genes.txt
    │       ├── 214766.AF0.75.snpEff.summary.html
    │       ├── 214766.AF0.75.snpEff.vcf.gz
    │       ├── 214766.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214766.AF0.75.snpSift.table.txt
    │       ├── 214766.snpEff.csv
    │       ├── 214766.snpEff.genes.txt
    │       ├── 214766.snpEff.summary.html
    │       ├── 214766.snpEff.vcf.gz
    │       ├── 214766.snpEff.vcf.gz.tbi
    │       ├── 214766.snpSift.table.txt
    │       ├── 214780.AF0.75.snpEff.csv
    │       ├── 214780.AF0.75.snpEff.genes.txt
    │       ├── 214780.AF0.75.snpEff.summary.html
    │       ├── 214780.AF0.75.snpEff.vcf.gz
    │       ├── 214780.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214780.AF0.75.snpSift.table.txt
    │       ├── 214780.snpEff.csv
    │       ├── 214780.snpEff.genes.txt
    │       ├── 214780.snpEff.summary.html
    │       ├── 214780.snpEff.vcf.gz
    │       ├── 214780.snpEff.vcf.gz.tbi
    │       ├── 214780.snpSift.table.txt
    │       ├── 214807.AF0.75.snpEff.csv
    │       ├── 214807.AF0.75.snpEff.genes.txt
    │       ├── 214807.AF0.75.snpEff.summary.html
    │       ├── 214807.AF0.75.snpEff.vcf.gz
    │       ├── 214807.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214807.AF0.75.snpSift.table.txt
    │       ├── 214807.snpEff.csv
    │       ├── 214807.snpEff.genes.txt
    │       ├── 214807.snpEff.summary.html
    │       ├── 214807.snpEff.vcf.gz
    │       ├── 214807.snpEff.vcf.gz.tbi
    │       ├── 214807.snpSift.table.txt
    │       ├── snpSift_template2.txt
    │       ├── snpSift_template_filtered.txt
    │       └── snpSift_template.txt
    ├── summary_variants_metrics_mqc.tsv
    └── varscan2
        ├── 214704.AF0.75.vcf.gz
        ├── 214704.AF0.75.vcf.gz.tbi
        ├── 214704.vcf.gz
        ├── 214704.vcf.gz.tbi
        ├── 214721.AF0.75.vcf.gz
        ├── 214721.AF0.75.vcf.gz.tbi
        ├── 214721.vcf.gz
        ├── 214721.vcf.gz.tbi
        ├── 214724.AF0.75.vcf.gz
        ├── 214724.AF0.75.vcf.gz.tbi
        ├── 214724.vcf.gz
        ├── 214724.vcf.gz.tbi
        ├── 214760.AF0.75.vcf.gz
        ├── 214760.AF0.75.vcf.gz.tbi
        ├── 214760.vcf.gz
        ├── 214760.vcf.gz.tbi
        ├── 214763.AF0.75.vcf.gz
        ├── 214763.AF0.75.vcf.gz.tbi
        ├── 214763.vcf.gz
        ├── 214763.vcf.gz.tbi
        ├── 214765.AF0.75.vcf.gz
        ├── 214765.AF0.75.vcf.gz.tbi
        ├── 214765.vcf.gz
        ├── 214765.vcf.gz.tbi
        ├── 214766.AF0.75.vcf.gz
        ├── 214766.AF0.75.vcf.gz.tbi
        ├── 214766.vcf.gz
        ├── 214766.vcf.gz.tbi
        ├── 214780.AF0.75.vcf.gz
        ├── 214780.AF0.75.vcf.gz.tbi
        ├── 214780.vcf.gz
        ├── 214780.vcf.gz.tbi
        ├── 214807.AF0.75.vcf.gz
        ├── 214807.AF0.75.vcf.gz.tbi
        ├── 214807.vcf.gz
        ├── 214807.vcf.gz.tbi
        ├── bcftools_stats
        │   ├── 214704.AF0.75.bcftools_stats.txt
        │   ├── 214704.bcftools_stats.txt
        │   ├── 214721.AF0.75.bcftools_stats.txt
        │   ├── 214721.bcftools_stats.txt
        │   ├── 214724.AF0.75.bcftools_stats.txt
        │   ├── 214724.bcftools_stats.txt
        │   ├── 214760.AF0.75.bcftools_stats.txt
        │   ├── 214760.bcftools_stats.txt
        │   ├── 214763.AF0.75.bcftools_stats.txt
        │   ├── 214763.bcftools_stats.txt
        │   ├── 214765.AF0.75.bcftools_stats.txt
        │   ├── 214765.bcftools_stats.txt
        │   ├── 214766.AF0.75.bcftools_stats.txt
        │   ├── 214766.bcftools_stats.txt
        │   ├── 214780.AF0.75.bcftools_stats.txt
        │   ├── 214780.bcftools_stats.txt
        │   ├── 214807.AF0.75.bcftools_stats.txt
        │   └── 214807.bcftools_stats.txt
        ├── consensus
        │   ├── 214704.AF0.75.consensus.masked.fa
        │   ├── 214721.AF0.75.consensus.masked.fa
        │   ├── 214724.AF0.75.consensus.masked.fa
        │   ├── 214760.AF0.75.consensus.masked.fa
        │   ├── 214763.AF0.75.consensus.masked.fa
        │   ├── 214765.AF0.75.consensus.masked.fa
        │   ├── 214766.AF0.75.consensus.masked.fa
        │   ├── 214780.AF0.75.consensus.masked.fa
        │   ├── 214807.AF0.75.consensus.masked.fa
        │   ├── base_qc
        │   │   ├── 214704.AF0.75.ACTG_density.pdf
        │   │   ├── 214704.AF0.75.base_counts.pdf
        │   │   ├── 214704.AF0.75.base_counts.tsv
        │   │   ├── 214704.AF0.75.N_density.pdf
        │   │   ├── 214704.AF0.75.N_run.tsv
        │   │   ├── 214721.AF0.75.ACTG_density.pdf
        │   │   ├── 214721.AF0.75.base_counts.pdf
        │   │   ├── 214721.AF0.75.base_counts.tsv
        │   │   ├── 214721.AF0.75.N_density.pdf
        │   │   ├── 214721.AF0.75.N_run.tsv
        │   │   ├── 214724.AF0.75.ACTG_density.pdf
        │   │   ├── 214724.AF0.75.base_counts.pdf
        │   │   ├── 214724.AF0.75.base_counts.tsv
        │   │   ├── 214724.AF0.75.N_density.pdf
        │   │   ├── 214724.AF0.75.N_run.tsv
        │   │   ├── 214760.AF0.75.ACTG_density.pdf
        │   │   ├── 214760.AF0.75.base_counts.pdf
        │   │   ├── 214760.AF0.75.base_counts.tsv
        │   │   ├── 214760.AF0.75.N_density.pdf
        │   │   ├── 214760.AF0.75.N_run.tsv
        │   │   ├── 214763.AF0.75.ACTG_density.pdf
        │   │   ├── 214763.AF0.75.base_counts.pdf
        │   │   ├── 214763.AF0.75.base_counts.tsv
        │   │   ├── 214763.AF0.75.N_density.pdf
        │   │   ├── 214763.AF0.75.N_run.tsv
        │   │   ├── 214765.AF0.75.ACTG_density.pdf
        │   │   ├── 214765.AF0.75.base_counts.pdf
        │   │   ├── 214765.AF0.75.base_counts.tsv
        │   │   ├── 214765.AF0.75.N_density.pdf
        │   │   ├── 214765.AF0.75.N_run.tsv
        │   │   ├── 214766.AF0.75.ACTG_density.pdf
        │   │   ├── 214766.AF0.75.base_counts.pdf
        │   │   ├── 214766.AF0.75.base_counts.tsv
        │   │   ├── 214766.AF0.75.N_density.pdf
        │   │   ├── 214766.AF0.75.N_run.tsv
        │   │   ├── 214780.AF0.75.ACTG_density.pdf
        │   │   ├── 214780.AF0.75.base_counts.pdf
        │   │   ├── 214780.AF0.75.base_counts.tsv
        │   │   ├── 214780.AF0.75.N_density.pdf
        │   │   ├── 214780.AF0.75.N_run.tsv
        │   │   ├── 214807.AF0.75.ACTG_density.pdf
        │   │   ├── 214807.AF0.75.base_counts.pdf
        │   │   ├── 214807.AF0.75.base_counts.tsv
        │   │   ├── 214807.AF0.75.N_density.pdf
        │   │   └── 214807.AF0.75.N_run.tsv
        │   └── snpSift_template.txt
        ├── log
        │   ├── 214704.varscan2.log
        │   ├── 214721.varscan2.log
        │   ├── 214724.varscan2.log
        │   ├── 214760.varscan2.log
        │   ├── 214763.varscan2.log
        │   ├── 214765.varscan2.log
        │   ├── 214766.varscan2.log
        │   ├── 214780.varscan2.log
        │   └── 214807.varscan2.log
        ├── quast
        │   └── AF0.75
        │       ├── aligned_stats
        │       │   ├── cumulative_plot.pdf
        │       │   ├── NAx_plot.pdf
        │       │   └── NGAx_plot.pdf
        │       ├── basic_stats
        │       │   ├── 214704.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214721.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214724.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214760.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214763.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214765.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214766.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214780.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214807.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── cumulative_plot.pdf
        │       │   ├── GC_content_plot.pdf
        │       │   ├── gc.icarus.txt
        │       │   ├── NGx_plot.pdf
        │       │   └── Nx_plot.pdf
        │       ├── contigs_reports
        │       │   ├── 214704_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214721_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214724_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214760_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214763_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214765_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214766_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214780_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214807_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── all_alignments_214704-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214721-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214724-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214760-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214763-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214765-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214766-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214780-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214807-AF0-75-consensus-masked.tsv
        │       │   ├── contigs_report_214704-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214704-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214704-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214704-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214721-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214721-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214721-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214721-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214724-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214724-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214724-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214724-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214760-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214760-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214760-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214760-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214763-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214763-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214763-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214763-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214765-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214765-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214765-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214765-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214766-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214766-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214766-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214766-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214780-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214780-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214780-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214780-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214807-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214807-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214807-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214807-AF0-75-consensus-masked.unaligned.info
        │       │   ├── minimap_output
        │       │   │   ├── 214704-AF0-75-consensus-masked.coords
        │       │   │   ├── 214704-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214704-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214704-AF0-75-consensus-masked.sf
        │       │   │   ├── 214704-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214704-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214721-AF0-75-consensus-masked.coords
        │       │   │   ├── 214721-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214721-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214721-AF0-75-consensus-masked.sf
        │       │   │   ├── 214721-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214721-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214724-AF0-75-consensus-masked.coords
        │       │   │   ├── 214724-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214724-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214724-AF0-75-consensus-masked.sf
        │       │   │   ├── 214724-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214724-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214760-AF0-75-consensus-masked.coords
        │       │   │   ├── 214760-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214760-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214760-AF0-75-consensus-masked.sf
        │       │   │   ├── 214760-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214760-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214763-AF0-75-consensus-masked.coords
        │       │   │   ├── 214763-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214763-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214763-AF0-75-consensus-masked.sf
        │       │   │   ├── 214763-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214763-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214765-AF0-75-consensus-masked.coords
        │       │   │   ├── 214765-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214765-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214765-AF0-75-consensus-masked.sf
        │       │   │   ├── 214765-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214765-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214766-AF0-75-consensus-masked.coords
        │       │   │   ├── 214766-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214766-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214766-AF0-75-consensus-masked.sf
        │       │   │   ├── 214766-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214766-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214780-AF0-75-consensus-masked.coords
        │       │   │   ├── 214780-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214780-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214780-AF0-75-consensus-masked.sf
        │       │   │   ├── 214780-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214780-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214807-AF0-75-consensus-masked.coords
        │       │   │   ├── 214807-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214807-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214807-AF0-75-consensus-masked.sf
        │       │   │   ├── 214807-AF0-75-consensus-masked.unaligned
        │       │   │   └── 214807-AF0-75-consensus-masked.used_snps.gz
        │       │   ├── misassemblies_frcurve_plot.pdf
        │       │   ├── misassemblies_plot.pdf
        │       │   ├── misassemblies_report.tex
        │       │   ├── misassemblies_report.tsv
        │       │   ├── misassemblies_report.txt
        │       │   ├── transposed_report_misassemblies.tex
        │       │   ├── transposed_report_misassemblies.tsv
        │       │   ├── transposed_report_misassemblies.txt
        │       │   ├── unaligned_report.tex
        │       │   ├── unaligned_report.tsv
        │       │   └── unaligned_report.txt
        │       ├── genome_stats
        │       │   ├── 214704-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214704-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214721-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214721-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214724-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214724-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214760-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214760-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214763-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214763-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214765-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214765-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214766-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214766-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214780-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214780-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214807-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214807-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── complete_features_histogram.pdf
        │       │   ├── features_cumulative_plot.pdf
        │       │   ├── features_frcurve_plot.pdf
        │       │   ├── genome_fraction_histogram.pdf
        │       │   └── genome_info.txt
        │       ├── icarus.html
        │       ├── icarus_viewers
        │       │   ├── alignment_viewer.html
        │       │   └── contig_size_viewer.html
        │       ├── quast.log
        │       ├── report.html
        │       ├── report.pdf
        │       ├── report.tex
        │       ├── report.tsv
        │       ├── report.txt
        │       ├── transposed_report.tex
        │       ├── transposed_report.tsv
        │       └── transposed_report.txt
        └── snpeff
            ├── 214704.AF0.75.snpEff.csv
            ├── 214704.AF0.75.snpEff.genes.txt
            ├── 214704.AF0.75.snpEff.summary.html
            ├── 214704.AF0.75.snpEff.vcf.gz
            ├── 214704.AF0.75.snpEff.vcf.gz.tbi
            ├── 214704.AF0.75.snpSift.table.txt
            ├── 214704.snpEff.csv
            ├── 214704.snpEff.genes.txt
            ├── 214704.snpEff.summary.html
            ├── 214704.snpEff.vcf.gz
            ├── 214704.snpEff.vcf.gz.tbi
            ├── 214704.snpSift.table.txt
            ├── 214721.AF0.75.snpEff.csv
            ├── 214721.AF0.75.snpEff.genes.txt
            ├── 214721.AF0.75.snpEff.summary.html
            ├── 214721.AF0.75.snpEff.vcf.gz
            ├── 214721.AF0.75.snpEff.vcf.gz.tbi
            ├── 214721.AF0.75.snpSift.table.txt
            ├── 214721.snpEff.csv
            ├── 214721.snpEff.genes.txt
            ├── 214721.snpEff.summary.html
            ├── 214721.snpEff.vcf.gz
            ├── 214721.snpEff.vcf.gz.tbi
            ├── 214721.snpSift.table.txt
            ├── 214724.AF0.75.snpEff.csv
            ├── 214724.AF0.75.snpEff.genes.txt
            ├── 214724.AF0.75.snpEff.summary.html
            ├── 214724.AF0.75.snpEff.vcf.gz
            ├── 214724.AF0.75.snpEff.vcf.gz.tbi
            ├── 214724.AF0.75.snpSift.table.txt
            ├── 214724.snpEff.csv
            ├── 214724.snpEff.genes.txt
            ├── 214724.snpEff.summary.html
            ├── 214724.snpEff.vcf.gz
            ├── 214724.snpEff.vcf.gz.tbi
            ├── 214724.snpSift.table.txt
            ├── 214760.AF0.75.snpEff.csv
            ├── 214760.AF0.75.snpEff.genes.txt
            ├── 214760.AF0.75.snpEff.summary.html
            ├── 214760.AF0.75.snpEff.vcf.gz
            ├── 214760.AF0.75.snpEff.vcf.gz.tbi
            ├── 214760.AF0.75.snpSift.table.txt
            ├── 214760.snpEff.csv
            ├── 214760.snpEff.genes.txt
            ├── 214760.snpEff.summary.html
            ├── 214760.snpEff.vcf.gz
            ├── 214760.snpEff.vcf.gz.tbi
            ├── 214760.snpSift.table.txt
            ├── 214763.AF0.75.snpEff.csv
            ├── 214763.AF0.75.snpEff.genes.txt
            ├── 214763.AF0.75.snpEff.summary.html
            ├── 214763.AF0.75.snpEff.vcf.gz
            ├── 214763.AF0.75.snpEff.vcf.gz.tbi
            ├── 214763.AF0.75.snpSift.table.txt
            ├── 214763.snpEff.csv
            ├── 214763.snpEff.genes.txt
            ├── 214763.snpEff.summary.html
            ├── 214763.snpEff.vcf.gz
            ├── 214763.snpEff.vcf.gz.tbi
            ├── 214763.snpSift.table.txt
            ├── 214765.AF0.75.snpEff.csv
            ├── 214765.AF0.75.snpEff.genes.txt
            ├── 214765.AF0.75.snpEff.summary.html
            ├── 214765.AF0.75.snpEff.vcf.gz
            ├── 214765.AF0.75.snpEff.vcf.gz.tbi
            ├── 214765.AF0.75.snpSift.table.txt
            ├── 214765.snpEff.csv
            ├── 214765.snpEff.genes.txt
            ├── 214765.snpEff.summary.html
            ├── 214765.snpEff.vcf.gz
            ├── 214765.snpEff.vcf.gz.tbi
            ├── 214765.snpSift.table.txt
            ├── 214766.AF0.75.snpEff.csv
            ├── 214766.AF0.75.snpEff.genes.txt
            ├── 214766.AF0.75.snpEff.summary.html
            ├── 214766.AF0.75.snpEff.vcf.gz
            ├── 214766.AF0.75.snpEff.vcf.gz.tbi
            ├── 214766.AF0.75.snpSift.table.txt
            ├── 214766.snpEff.csv
            ├── 214766.snpEff.genes.txt
            ├── 214766.snpEff.summary.html
            ├── 214766.snpEff.vcf.gz
            ├── 214766.snpEff.vcf.gz.tbi
            ├── 214766.snpSift.table.txt
            ├── 214780.AF0.75.snpEff.csv
            ├── 214780.AF0.75.snpEff.genes.txt
            ├── 214780.AF0.75.snpEff.summary.html
            ├── 214780.AF0.75.snpEff.vcf.gz
            ├── 214780.AF0.75.snpEff.vcf.gz.tbi
            ├── 214780.AF0.75.snpSift.table.txt
            ├── 214780.snpEff.csv
            ├── 214780.snpEff.genes.txt
            ├── 214780.snpEff.summary.html
            ├── 214780.snpEff.vcf.gz
            ├── 214780.snpEff.vcf.gz.tbi
            ├── 214780.snpSift.table.txt
            ├── 214807.AF0.75.snpEff.csv
            ├── 214807.AF0.75.snpEff.genes.txt
            ├── 214807.AF0.75.snpEff.summary.html
            ├── 214807.AF0.75.snpEff.vcf.gz
            ├── 214807.AF0.75.snpEff.vcf.gz.tbi
            ├── 214807.AF0.75.snpSift.table.txt
            ├── 214807.snpEff.csv
            ├── 214807.snpEff.genes.txt
            ├── 214807.snpEff.summary.html
            ├── 214807.snpEff.vcf.gz
            ├── 214807.snpEff.vcf.gz.tbi
            └── 214807.snpSift.table.txt


├── assembly
│   ├── cutadapt
│   │   ├── fastqc
│   │   │   ├── 214704_1.ptrim_fastqc.html
│   │   │   ├── 214704_2.ptrim_fastqc.html
│   │   │   ├── 214721_1.ptrim_fastqc.html
│   │   │   ├── 214721_2.ptrim_fastqc.html
│   │   │   ├── 214724_1.ptrim_fastqc.html
│   │   │   ├── 214724_2.ptrim_fastqc.html
│   │   │   ├── 214760_1.ptrim_fastqc.html
│   │   │   ├── 214760_2.ptrim_fastqc.html
│   │   │   ├── 214763_1.ptrim_fastqc.html
│   │   │   ├── 214763_2.ptrim_fastqc.html
│   │   │   ├── 214765_1.ptrim_fastqc.html
│   │   │   ├── 214765_2.ptrim_fastqc.html
│   │   │   ├── 214766_1.ptrim_fastqc.html
│   │   │   ├── 214766_2.ptrim_fastqc.html
│   │   │   ├── 214780_1.ptrim_fastqc.html
│   │   │   ├── 214780_2.ptrim_fastqc.html
│   │   │   ├── 214807_1.ptrim_fastqc.html
│   │   │   ├── 214807_2.ptrim_fastqc.html
│   │   │   └── zips
│   │   │       ├── 214704_1.ptrim_fastqc.zip
│   │   │       ├── 214704_2.ptrim_fastqc.zip
│   │   │       ├── 214721_1.ptrim_fastqc.zip
│   │   │       ├── 214721_2.ptrim_fastqc.zip
│   │   │       ├── 214724_1.ptrim_fastqc.zip
│   │   │       ├── 214724_2.ptrim_fastqc.zip
│   │   │       ├── 214760_1.ptrim_fastqc.zip
│   │   │       ├── 214760_2.ptrim_fastqc.zip
│   │   │       ├── 214763_1.ptrim_fastqc.zip
│   │   │       ├── 214763_2.ptrim_fastqc.zip
│   │   │       ├── 214765_1.ptrim_fastqc.zip
│   │   │       ├── 214765_2.ptrim_fastqc.zip
│   │   │       ├── 214766_1.ptrim_fastqc.zip
│   │   │       ├── 214766_2.ptrim_fastqc.zip
│   │   │       ├── 214780_1.ptrim_fastqc.zip
│   │   │       ├── 214780_2.ptrim_fastqc.zip
│   │   │       ├── 214807_1.ptrim_fastqc.zip
│   │   │       └── 214807_2.ptrim_fastqc.zip
│   │   └── log
│   │       ├── 214704.cutadapt.log
│   │       ├── 214721.cutadapt.log
│   │       ├── 214724.cutadapt.log
│   │       ├── 214760.cutadapt.log
│   │       ├── 214763.cutadapt.log
│   │       ├── 214765.cutadapt.log
│   │       ├── 214766.cutadapt.log
│   │       ├── 214780.cutadapt.log
│   │       └── 214807.cutadapt.log
│   ├── kraken2
│   │   ├── 214704.kraken2.report.txt
│   │   ├── 214721.kraken2.report.txt
│   │   ├── 214724.kraken2.report.txt
│   │   ├── 214760.kraken2.report.txt
│   │   ├── 214763.kraken2.report.txt
│   │   ├── 214765.kraken2.report.txt
│   │   ├── 214766.kraken2.report.txt
│   │   ├── 214780.kraken2.report.txt
│   │   └── 214807.kraken2.report.txt
│   └── summary_assembly_metrics_mqc.tsv
├── multiqc
│   ├── multiqc_data
│   │   ├── mqc_cutadapt_trimmed_sequences_plot_Counts.yaml
│   │   ├── mqc_cutadapt_trimmed_sequences_plot_Obs_Exp.yaml
│   │   ├── multiqc_bcftools_stats_bcftools_bcftools.yaml
│   │   ├── multiqc_bcftools_stats_bcftools_ivar.yaml
│   │   ├── multiqc_bcftools_stats_bcftools_varscan2.yaml
│   │   ├── multiqc_bowtie2.yaml
│   │   ├── multiqc_cutadapt.yaml
│   │   ├── multiqc_data.json
│   │   ├── multiqc_de_novo_assembly_metrics.yaml
│   │   ├── multiqc_fastp.yaml
│   │   ├── multiqc_fastqc_fastqc_cutadapt.yaml
│   │   ├── multiqc_fastqc_fastqc_fastp.yaml
│   │   ├── multiqc_fastqc_fastqc_raw.yaml
│   │   ├── multiqc_general_stats.yaml
│   │   ├── multiqc_ivar_primers.yaml
│   │   ├── multiqc_ivar_summary.yaml
│   │   ├── multiqc.log
│   │   ├── multiqc_picard_AlignmentSummaryMetrics.yaml
│   │   ├── multiqc_picard_insertSize.yaml
│   │   ├── multiqc_picard_wgsmetrics.yaml
│   │   ├── multiqc_quast_quast_bcftools.yaml
│   │   ├── multiqc_quast_quast_ivar.yaml
│   │   ├── multiqc_quast_quast_varscan2.yaml
│   │   ├── multiqc_samtools_flagstat_samtools_bowtie2.yaml
│   │   ├── multiqc_samtools_flagstat_samtools_ivar.yaml
│   │   ├── multiqc_samtools_idxstats_samtools_bowtie2.yaml
│   │   ├── multiqc_samtools_idxstats_samtools_ivar.yaml
│   │   ├── multiqc_samtools_stats_samtools_bowtie2.yaml
│   │   ├── multiqc_samtools_stats_samtools_ivar.yaml
│   │   ├── multiqc_snpeff_snpeff_bcftools.yaml
│   │   ├── multiqc_snpeff_snpeff_ivar.yaml
│   │   ├── multiqc_snpeff_snpeff_varscan2.yaml
│   │   ├── multiqc_sources.yaml
│   │   ├── multiqc_variant_calling_metrics.yaml
│   │   └── multiqc_varscan2_summary.yaml
│   └── multiqc_report.html
├── pipeline_info
│   ├── execution_report.html
│   ├── execution_timeline.html
│   ├── execution_trace.txt
│   ├── pipeline_dag.svg
│   ├── pipeline_report.html
│   ├── pipeline_report.txt
│   ├── results_description.html
│   ├── samplesheet.valid.csv
│   └── software_versions.csv
├── preprocess
│   ├── fastp
│   │   ├── 214704.fastp.html
│   │   ├── 214704.fastp.json
│   │   ├── 214721.fastp.html
│   │   ├── 214721.fastp.json
│   │   ├── 214724.fastp.html
│   │   ├── 214724.fastp.json
│   │   ├── 214760.fastp.html
│   │   ├── 214760.fastp.json
│   │   ├── 214763.fastp.html
│   │   ├── 214763.fastp.json
│   │   ├── 214765.fastp.html
│   │   ├── 214765.fastp.json
│   │   ├── 214766.fastp.html
│   │   ├── 214766.fastp.json
│   │   ├── 214780.fastp.html
│   │   ├── 214780.fastp.json
│   │   ├── 214807.fastp.html
│   │   ├── 214807.fastp.json
│   │   ├── fastqc
│   │   │   ├── 214704_1.trim_fastqc.html
│   │   │   ├── 214704_2.trim_fastqc.html
│   │   │   ├── 214721_1.trim_fastqc.html
│   │   │   ├── 214721_2.trim_fastqc.html
│   │   │   ├── 214724_1.trim_fastqc.html
│   │   │   ├── 214724_2.trim_fastqc.html
│   │   │   ├── 214760_1.trim_fastqc.html
│   │   │   ├── 214760_2.trim_fastqc.html
│   │   │   ├── 214763_1.trim_fastqc.html
│   │   │   ├── 214763_2.trim_fastqc.html
│   │   │   ├── 214765_1.trim_fastqc.html
│   │   │   ├── 214765_2.trim_fastqc.html
│   │   │   ├── 214766_1.trim_fastqc.html
│   │   │   ├── 214766_2.trim_fastqc.html
│   │   │   ├── 214780_1.trim_fastqc.html
│   │   │   ├── 214780_2.trim_fastqc.html
│   │   │   ├── 214807_1.trim_fastqc.html
│   │   │   ├── 214807_2.trim_fastqc.html
│   │   │   └── zips
│   │   │       ├── 214704_1.trim_fastqc.zip
│   │   │       ├── 214704_2.trim_fastqc.zip
│   │   │       ├── 214721_1.trim_fastqc.zip
│   │   │       ├── 214721_2.trim_fastqc.zip
│   │   │       ├── 214724_1.trim_fastqc.zip
│   │   │       ├── 214724_2.trim_fastqc.zip
│   │   │       ├── 214760_1.trim_fastqc.zip
│   │   │       ├── 214760_2.trim_fastqc.zip
│   │   │       ├── 214763_1.trim_fastqc.zip
│   │   │       ├── 214763_2.trim_fastqc.zip
│   │   │       ├── 214765_1.trim_fastqc.zip
│   │   │       ├── 214765_2.trim_fastqc.zip
│   │   │       ├── 214766_1.trim_fastqc.zip
│   │   │       ├── 214766_2.trim_fastqc.zip
│   │   │       ├── 214780_1.trim_fastqc.zip
│   │   │       ├── 214780_2.trim_fastqc.zip
│   │   │       ├── 214807_1.trim_fastqc.zip
│   │   │       └── 214807_2.trim_fastqc.zip
│   │   └── log
│   │       ├── 214704.fastp.log
│   │       ├── 214721.fastp.log
│   │       ├── 214724.fastp.log
│   │       ├── 214760.fastp.log
│   │       ├── 214763.fastp.log
│   │       ├── 214765.fastp.log
│   │       ├── 214766.fastp.log
│   │       ├── 214780.fastp.log
│   │       └── 214807.fastp.log
│   └── fastqc
│       ├── 214704_1.merged_fastqc.html
│       ├── 214704_2.merged_fastqc.html
│       ├── 214721_1.merged_fastqc.html
│       ├── 214721_2.merged_fastqc.html
│       ├── 214724_1.merged_fastqc.html
│       ├── 214724_2.merged_fastqc.html
│       ├── 214760_1.merged_fastqc.html
│       ├── 214760_2.merged_fastqc.html
│       ├── 214763_1.merged_fastqc.html
│       ├── 214763_2.merged_fastqc.html
│       ├── 214765_1.merged_fastqc.html
│       ├── 214765_2.merged_fastqc.html
│       ├── 214766_1.merged_fastqc.html
│       ├── 214766_2.merged_fastqc.html
│       ├── 214780_1.merged_fastqc.html
│       ├── 214780_2.merged_fastqc.html
│       ├── 214807_1.merged_fastqc.html
│       ├── 214807_2.merged_fastqc.html
│       └── zips
│           ├── 214704_1.merged_fastqc.zip
│           ├── 214704_2.merged_fastqc.zip
│           ├── 214721_1.merged_fastqc.zip
│           ├── 214721_2.merged_fastqc.zip
│           ├── 214724_1.merged_fastqc.zip
│           ├── 214724_2.merged_fastqc.zip
│           ├── 214760_1.merged_fastqc.zip
│           ├── 214760_2.merged_fastqc.zip
│           ├── 214763_1.merged_fastqc.zip
│           ├── 214763_2.merged_fastqc.zip
│           ├── 214765_1.merged_fastqc.zip
│           ├── 214765_2.merged_fastqc.zip
│           ├── 214766_1.merged_fastqc.zip
│           ├── 214766_2.merged_fastqc.zip
│           ├── 214780_1.merged_fastqc.zip
│           ├── 214780_2.merged_fastqc.zip
│           ├── 214807_1.merged_fastqc.zip
│           └── 214807_2.merged_fastqc.zip
└── variants
    ├── bam
    │   ├── 214704.bam
    │   ├── 214704.sorted.bam
    │   ├── 214704.sorted.bam.bai
    │   ├── 214704.trim.sorted.bam
    │   ├── 214704.trim.sorted.bam.bai
    │   ├── 214721.bam
    │   ├── 214721.sorted.bam
    │   ├── 214721.sorted.bam.bai
    │   ├── 214721.trim.sorted.bam
    │   ├── 214721.trim.sorted.bam.bai
    │   ├── 214724.bam
    │   ├── 214724.sorted.bam
    │   ├── 214724.sorted.bam.bai
    │   ├── 214724.trim.sorted.bam
    │   ├── 214724.trim.sorted.bam.bai
    │   ├── 214760.bam
    │   ├── 214760.sorted.bam
    │   ├── 214760.sorted.bam.bai
    │   ├── 214760.trim.sorted.bam
    │   ├── 214760.trim.sorted.bam.bai
    │   ├── 214763.bam
    │   ├── 214763.sorted.bam
    │   ├── 214763.sorted.bam.bai
    │   ├── 214763.trim.sorted.bam
    │   ├── 214763.trim.sorted.bam.bai
    │   ├── 214765.bam
    │   ├── 214765.sorted.bam
    │   ├── 214765.sorted.bam.bai
    │   ├── 214765.trim.sorted.bam
    │   ├── 214765.trim.sorted.bam.bai
    │   ├── 214766.bam
    │   ├── 214766.sorted.bam
    │   ├── 214766.sorted.bam.bai
    │   ├── 214766.trim.sorted.bam
    │   ├── 214766.trim.sorted.bam.bai
    │   ├── 214780.bam
    │   ├── 214780.sorted.bam
    │   ├── 214780.sorted.bam.bai
    │   ├── 214780.trim.sorted.bam
    │   ├── 214780.trim.sorted.bam.bai
    │   ├── 214807.bam
    │   ├── 214807.sorted.bam
    │   ├── 214807.sorted.bam.bai
    │   ├── 214807.trim.sorted.bam
    │   ├── 214807.trim.sorted.bam.bai
    │   ├── log
    │   │   ├── 214704.bowtie2.log
    │   │   ├── 214704.trim.ivar.log
    │   │   ├── 214721.bowtie2.log
    │   │   ├── 214721.trim.ivar.log
    │   │   ├── 214724.bowtie2.log
    │   │   ├── 214724.trim.ivar.log
    │   │   ├── 214760.bowtie2.log
    │   │   ├── 214760.trim.ivar.log
    │   │   ├── 214763.bowtie2.log
    │   │   ├── 214763.trim.ivar.log
    │   │   ├── 214765.bowtie2.log
    │   │   ├── 214765.trim.ivar.log
    │   │   ├── 214766.bowtie2.log
    │   │   ├── 214766.trim.ivar.log
    │   │   ├── 214780.bowtie2.log
    │   │   ├── 214780.trim.ivar.log
    │   │   ├── 214807.bowtie2.log
    │   │   └── 214807.trim.ivar.log
    │   ├── mosdepth
    │   │   ├── amplicon
    │   │   │   ├── 214704.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214704.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214704.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214704.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214704.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214704.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214704.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214704.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214704.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214721.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214721.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214721.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214721.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214721.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214721.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214721.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214721.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214721.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214724.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214724.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214724.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214724.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214724.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214724.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214724.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214724.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214724.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214760.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214760.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214760.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214760.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214760.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214760.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214760.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214760.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214760.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214763.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214763.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214763.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214763.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214763.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214763.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214763.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214763.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214763.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214765.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214765.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214765.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214765.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214765.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214765.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214765.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214765.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214765.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214766.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214766.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214766.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214766.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214766.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214766.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214766.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214766.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214766.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214780.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214780.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214780.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214780.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214780.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214780.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214780.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214780.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214780.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   ├── 214807.trim.amplicon.mosdepth.global.dist.txt
    │   │   │   ├── 214807.trim.amplicon.mosdepth.region.dist.txt
    │   │   │   ├── 214807.trim.amplicon.mosdepth.summary.txt
    │   │   │   ├── 214807.trim.amplicon.per-base.bed.gz
    │   │   │   ├── 214807.trim.amplicon.per-base.bed.gz.csi
    │   │   │   ├── 214807.trim.amplicon.regions.bed.gz
    │   │   │   ├── 214807.trim.amplicon.regions.bed.gz.csi
    │   │   │   ├── 214807.trim.amplicon.thresholds.bed.gz
    │   │   │   ├── 214807.trim.amplicon.thresholds.bed.gz.csi
    │   │   │   └── plots
    │   │   │       ├── 214704.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214704.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214721.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214721.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214724.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214724.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214760.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214760.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214763.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214763.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214765.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214765.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214766.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214766.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214780.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214780.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── 214807.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214807.trim.amplicon.regions.coverage.tsv
    │   │   │       ├── all_samples.trim.amplicon.regions.coverage.tsv
    │   │   │       └── all_samples.trim.amplicon.regions.heatmap.pdf
    │   │   └── genome
    │   │       ├── 214704.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214704.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214704.trim.genome.mosdepth.summary.txt
    │   │       ├── 214704.trim.genome.per-base.bed.gz
    │   │       ├── 214704.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214704.trim.genome.regions.bed.gz
    │   │       ├── 214704.trim.genome.regions.bed.gz.csi
    │   │       ├── 214721.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214721.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214721.trim.genome.mosdepth.summary.txt
    │   │       ├── 214721.trim.genome.per-base.bed.gz
    │   │       ├── 214721.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214721.trim.genome.regions.bed.gz
    │   │       ├── 214721.trim.genome.regions.bed.gz.csi
    │   │       ├── 214724.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214724.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214724.trim.genome.mosdepth.summary.txt
    │   │       ├── 214724.trim.genome.per-base.bed.gz
    │   │       ├── 214724.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214724.trim.genome.regions.bed.gz
    │   │       ├── 214724.trim.genome.regions.bed.gz.csi
    │   │       ├── 214760.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214760.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214760.trim.genome.mosdepth.summary.txt
    │   │       ├── 214760.trim.genome.per-base.bed.gz
    │   │       ├── 214760.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214760.trim.genome.regions.bed.gz
    │   │       ├── 214760.trim.genome.regions.bed.gz.csi
    │   │       ├── 214763.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214763.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214763.trim.genome.mosdepth.summary.txt
    │   │       ├── 214763.trim.genome.per-base.bed.gz
    │   │       ├── 214763.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214763.trim.genome.regions.bed.gz
    │   │       ├── 214763.trim.genome.regions.bed.gz.csi
    │   │       ├── 214765.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214765.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214765.trim.genome.mosdepth.summary.txt
    │   │       ├── 214765.trim.genome.per-base.bed.gz
    │   │       ├── 214765.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214765.trim.genome.regions.bed.gz
    │   │       ├── 214765.trim.genome.regions.bed.gz.csi
    │   │       ├── 214766.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214766.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214766.trim.genome.mosdepth.summary.txt
    │   │       ├── 214766.trim.genome.per-base.bed.gz
    │   │       ├── 214766.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214766.trim.genome.regions.bed.gz
    │   │       ├── 214766.trim.genome.regions.bed.gz.csi
    │   │       ├── 214780.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214780.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214780.trim.genome.mosdepth.summary.txt
    │   │       ├── 214780.trim.genome.per-base.bed.gz
    │   │       ├── 214780.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214780.trim.genome.regions.bed.gz
    │   │       ├── 214780.trim.genome.regions.bed.gz.csi
    │   │       ├── 214807.trim.genome.mosdepth.global.dist.txt
    │   │       ├── 214807.trim.genome.mosdepth.region.dist.txt
    │   │       ├── 214807.trim.genome.mosdepth.summary.txt
    │   │       ├── 214807.trim.genome.per-base.bed.gz
    │   │       ├── 214807.trim.genome.per-base.bed.gz.csi
    │   │       ├── 214807.trim.genome.regions.bed.gz
    │   │       ├── 214807.trim.genome.regions.bed.gz.csi
    │   │       └── plots
    │   │           ├── 214704.trim.genome.regions.coverage.pdf
    │   │           ├── 214704.trim.genome.regions.coverage.tsv
    │   │           ├── 214721.trim.genome.regions.coverage.pdf
    │   │           ├── 214721.trim.genome.regions.coverage.tsv
    │   │           ├── 214724.trim.genome.regions.coverage.pdf
    │   │           ├── 214724.trim.genome.regions.coverage.tsv
    │   │           ├── 214760.trim.genome.regions.coverage.pdf
    │   │           ├── 214760.trim.genome.regions.coverage.tsv
    │   │           ├── 214763.trim.genome.regions.coverage.pdf
    │   │           ├── 214763.trim.genome.regions.coverage.tsv
    │   │           ├── 214765.trim.genome.regions.coverage.pdf
    │   │           ├── 214765.trim.genome.regions.coverage.tsv
    │   │           ├── 214766.trim.genome.regions.coverage.pdf
    │   │           ├── 214766.trim.genome.regions.coverage.tsv
    │   │           ├── 214780.trim.genome.regions.coverage.pdf
    │   │           ├── 214780.trim.genome.regions.coverage.tsv
    │   │           ├── 214807.trim.genome.regions.coverage.pdf
    │   │           ├── 214807.trim.genome.regions.coverage.tsv
    │   │           └── all_samples.trim.genome.regions.coverage.tsv
    │   ├── mpileup
    │   │   ├── 214704.trim.mpileup
    │   │   ├── 214721.trim.mpileup
    │   │   ├── 214724.trim.mpileup
    │   │   ├── 214760.trim.mpileup
    │   │   ├── 214763.trim.mpileup
    │   │   ├── 214765.trim.mpileup
    │   │   ├── 214766.trim.mpileup
    │   │   ├── 214780.trim.mpileup
    │   │   └── 214807.trim.mpileup
    │   ├── picard_metrics
    │   │   ├── 214704.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214704.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214704.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214704.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214704.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214704.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214704.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214704.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214704.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214704.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214721.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214721.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214721.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214721.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214721.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214721.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214721.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214721.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214721.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214721.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214724.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214724.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214724.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214724.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214724.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214724.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214724.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214724.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214724.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214724.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214760.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214760.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214760.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214760.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214760.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214760.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214760.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214760.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214760.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214760.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214763.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214763.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214763.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214763.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214763.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214763.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214763.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214763.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214763.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214763.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214765.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214765.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214765.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214765.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214765.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214765.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214765.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214765.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214765.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214765.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214766.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214766.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214766.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214766.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214766.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214766.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214766.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214766.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214766.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214766.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214780.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214780.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214780.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214780.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214780.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214780.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214780.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214780.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214780.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   ├── 214780.trim.CollectWgsMetrics.coverage_metrics
    │   │   ├── 214807.trim.CollectMultipleMetrics.alignment_summary_metrics
    │   │   ├── 214807.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
    │   │   ├── 214807.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
    │   │   ├── 214807.trim.CollectMultipleMetrics.insert_size_histogram.pdf
    │   │   ├── 214807.trim.CollectMultipleMetrics.insert_size_metrics
    │   │   ├── 214807.trim.CollectMultipleMetrics.quality_by_cycle_metrics
    │   │   ├── 214807.trim.CollectMultipleMetrics.quality_by_cycle.pdf
    │   │   ├── 214807.trim.CollectMultipleMetrics.quality_distribution_metrics
    │   │   ├── 214807.trim.CollectMultipleMetrics.quality_distribution.pdf
    │   │   └── 214807.trim.CollectWgsMetrics.coverage_metrics
    │   └── samtools_stats
    │       ├── 214704.sorted.bam.flagstat
    │       ├── 214704.sorted.bam.idxstats
    │       ├── 214704.sorted.bam.stats
    │       ├── 214704.trim.sorted.bam.flagstat
    │       ├── 214704.trim.sorted.bam.idxstats
    │       ├── 214704.trim.sorted.bam.stats
    │       ├── 214721.sorted.bam.flagstat
    │       ├── 214721.sorted.bam.idxstats
    │       ├── 214721.sorted.bam.stats
    │       ├── 214721.trim.sorted.bam.flagstat
    │       ├── 214721.trim.sorted.bam.idxstats
    │       ├── 214721.trim.sorted.bam.stats
    │       ├── 214724.sorted.bam.flagstat
    │       ├── 214724.sorted.bam.idxstats
    │       ├── 214724.sorted.bam.stats
    │       ├── 214724.trim.sorted.bam.flagstat
    │       ├── 214724.trim.sorted.bam.idxstats
    │       ├── 214724.trim.sorted.bam.stats
    │       ├── 214760.sorted.bam.flagstat
    │       ├── 214760.sorted.bam.idxstats
    │       ├── 214760.sorted.bam.stats
    │       ├── 214760.trim.sorted.bam.flagstat
    │       ├── 214760.trim.sorted.bam.idxstats
    │       ├── 214760.trim.sorted.bam.stats
    │       ├── 214763.sorted.bam.flagstat
    │       ├── 214763.sorted.bam.idxstats
    │       ├── 214763.sorted.bam.stats
    │       ├── 214763.trim.sorted.bam.flagstat
    │       ├── 214763.trim.sorted.bam.idxstats
    │       ├── 214763.trim.sorted.bam.stats
    │       ├── 214765.sorted.bam.flagstat
    │       ├── 214765.sorted.bam.idxstats
    │       ├── 214765.sorted.bam.stats
    │       ├── 214765.trim.sorted.bam.flagstat
    │       ├── 214765.trim.sorted.bam.idxstats
    │       ├── 214765.trim.sorted.bam.stats
    │       ├── 214766.sorted.bam.flagstat
    │       ├── 214766.sorted.bam.idxstats
    │       ├── 214766.sorted.bam.stats
    │       ├── 214766.trim.sorted.bam.flagstat
    │       ├── 214766.trim.sorted.bam.idxstats
    │       ├── 214766.trim.sorted.bam.stats
    │       ├── 214780.sorted.bam.flagstat
    │       ├── 214780.sorted.bam.idxstats
    │       ├── 214780.sorted.bam.stats
    │       ├── 214780.trim.sorted.bam.flagstat
    │       ├── 214780.trim.sorted.bam.idxstats
    │       ├── 214780.trim.sorted.bam.stats
    │       ├── 214807.sorted.bam.flagstat
    │       ├── 214807.sorted.bam.idxstats
    │       ├── 214807.sorted.bam.stats
    │       ├── 214807.trim.sorted.bam.flagstat
    │       ├── 214807.trim.sorted.bam.idxstats
    │       └── 214807.trim.sorted.bam.stats
    ├── bcftools
    │   ├── 214704.vcf.gz
    │   ├── 214704.vcf.gz.tbi
    │   ├── 214721.vcf.gz
    │   ├── 214721.vcf.gz.tbi
    │   ├── 214724.vcf.gz
    │   ├── 214724.vcf.gz.tbi
    │   ├── 214760.vcf.gz
    │   ├── 214760.vcf.gz.tbi
    │   ├── 214763.vcf.gz
    │   ├── 214763.vcf.gz.tbi
    │   ├── 214765.vcf.gz
    │   ├── 214765.vcf.gz.tbi
    │   ├── 214766.vcf.gz
    │   ├── 214766.vcf.gz.tbi
    │   ├── 214780.vcf.gz
    │   ├── 214780.vcf.gz.tbi
    │   ├── 214807.vcf.gz
    │   ├── 214807.vcf.gz.tbi
    │   ├── bcftools_stats
    │   │   ├── 214704.bcftools_stats.txt
    │   │   ├── 214721.bcftools_stats.txt
    │   │   ├── 214724.bcftools_stats.txt
    │   │   ├── 214760.bcftools_stats.txt
    │   │   ├── 214763.bcftools_stats.txt
    │   │   ├── 214765.bcftools_stats.txt
    │   │   ├── 214766.bcftools_stats.txt
    │   │   ├── 214780.bcftools_stats.txt
    │   │   └── 214807.bcftools_stats.txt
    │   ├── consensus
    │   │   ├── 214704.consensus.masked.fa
    │   │   ├── 214721.consensus.masked.fa
    │   │   ├── 214724.consensus.masked.fa
    │   │   ├── 214760.consensus.masked.fa
    │   │   ├── 214763.consensus.masked.fa
    │   │   ├── 214765.consensus.masked.fa
    │   │   ├── 214766.consensus.masked.fa
    │   │   ├── 214780.consensus.masked.fa
    │   │   ├── 214807.consensus.masked.fa
    │   │   └── base_qc
    │   │       ├── 214704.ACTG_density.pdf
    │   │       ├── 214704.base_counts.pdf
    │   │       ├── 214704.base_counts.tsv
    │   │       ├── 214704.N_density.pdf
    │   │       ├── 214704.N_run.tsv
    │   │       ├── 214721.ACTG_density.pdf
    │   │       ├── 214721.base_counts.pdf
    │   │       ├── 214721.base_counts.tsv
    │   │       ├── 214721.N_density.pdf
    │   │       ├── 214721.N_run.tsv
    │   │       ├── 214724.ACTG_density.pdf
    │   │       ├── 214724.base_counts.pdf
    │   │       ├── 214724.base_counts.tsv
    │   │       ├── 214724.N_density.pdf
    │   │       ├── 214724.N_run.tsv
    │   │       ├── 214760.ACTG_density.pdf
    │   │       ├── 214760.base_counts.pdf
    │   │       ├── 214760.base_counts.tsv
    │   │       ├── 214760.N_density.pdf
    │   │       ├── 214760.N_run.tsv
    │   │       ├── 214763.ACTG_density.pdf
    │   │       ├── 214763.base_counts.pdf
    │   │       ├── 214763.base_counts.tsv
    │   │       ├── 214763.N_density.pdf
    │   │       ├── 214763.N_run.tsv
    │   │       ├── 214765.ACTG_density.pdf
    │   │       ├── 214765.base_counts.pdf
    │   │       ├── 214765.base_counts.tsv
    │   │       ├── 214765.N_density.pdf
    │   │       ├── 214765.N_run.tsv
    │   │       ├── 214766.ACTG_density.pdf
    │   │       ├── 214766.base_counts.pdf
    │   │       ├── 214766.base_counts.tsv
    │   │       ├── 214766.N_density.pdf
    │   │       ├── 214766.N_run.tsv
    │   │       ├── 214780.ACTG_density.pdf
    │   │       ├── 214780.base_counts.pdf
    │   │       ├── 214780.base_counts.tsv
    │   │       ├── 214780.N_density.pdf
    │   │       ├── 214780.N_run.tsv
    │   │       ├── 214807.ACTG_density.pdf
    │   │       ├── 214807.base_counts.pdf
    │   │       ├── 214807.base_counts.tsv
    │   │       ├── 214807.N_density.pdf
    │   │       └── 214807.N_run.tsv
    │   ├── quast
    │   │   ├── aligned_stats
    │   │   │   ├── cumulative_plot.pdf
    │   │   │   ├── NAx_plot.pdf
    │   │   │   └── NGAx_plot.pdf
    │   │   ├── basic_stats
    │   │   │   ├── 214704.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214721.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214724.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214760.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214763.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214765.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214766.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214780.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── 214807.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── cumulative_plot.pdf
    │   │   │   ├── GC_content_plot.pdf
    │   │   │   ├── gc.icarus.txt
    │   │   │   ├── NGx_plot.pdf
    │   │   │   └── Nx_plot.pdf
    │   │   ├── contigs_reports
    │   │   │   ├── 214704_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214721_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214724_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214760_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214763_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214765_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214766_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214780_consensus_masked.mis_contigs.fa
    │   │   │   ├── 214807_consensus_masked.mis_contigs.fa
    │   │   │   ├── all_alignments_214704-consensus-masked.tsv
    │   │   │   ├── all_alignments_214721-consensus-masked.tsv
    │   │   │   ├── all_alignments_214724-consensus-masked.tsv
    │   │   │   ├── all_alignments_214760-consensus-masked.tsv
    │   │   │   ├── all_alignments_214763-consensus-masked.tsv
    │   │   │   ├── all_alignments_214765-consensus-masked.tsv
    │   │   │   ├── all_alignments_214766-consensus-masked.tsv
    │   │   │   ├── all_alignments_214780-consensus-masked.tsv
    │   │   │   ├── all_alignments_214807-consensus-masked.tsv
    │   │   │   ├── contigs_report_214704-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214704-consensus-masked.stderr
    │   │   │   ├── contigs_report_214704-consensus-masked.stdout
    │   │   │   ├── contigs_report_214704-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214721-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214721-consensus-masked.stderr
    │   │   │   ├── contigs_report_214721-consensus-masked.stdout
    │   │   │   ├── contigs_report_214721-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214724-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214724-consensus-masked.stderr
    │   │   │   ├── contigs_report_214724-consensus-masked.stdout
    │   │   │   ├── contigs_report_214724-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214760-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214760-consensus-masked.stderr
    │   │   │   ├── contigs_report_214760-consensus-masked.stdout
    │   │   │   ├── contigs_report_214760-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214763-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214763-consensus-masked.stderr
    │   │   │   ├── contigs_report_214763-consensus-masked.stdout
    │   │   │   ├── contigs_report_214763-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214765-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214765-consensus-masked.stderr
    │   │   │   ├── contigs_report_214765-consensus-masked.stdout
    │   │   │   ├── contigs_report_214765-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214766-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214766-consensus-masked.stderr
    │   │   │   ├── contigs_report_214766-consensus-masked.stdout
    │   │   │   ├── contigs_report_214766-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214780-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214780-consensus-masked.stderr
    │   │   │   ├── contigs_report_214780-consensus-masked.stdout
    │   │   │   ├── contigs_report_214780-consensus-masked.unaligned.info
    │   │   │   ├── contigs_report_214807-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214807-consensus-masked.stderr
    │   │   │   ├── contigs_report_214807-consensus-masked.stdout
    │   │   │   ├── contigs_report_214807-consensus-masked.unaligned.info
    │   │   │   ├── minimap_output
    │   │   │   │   ├── 214704-consensus-masked.coords
    │   │   │   │   ├── 214704-consensus-masked.coords.filtered
    │   │   │   │   ├── 214704-consensus-masked.coords_tmp
    │   │   │   │   ├── 214704-consensus-masked.sf
    │   │   │   │   ├── 214704-consensus-masked.unaligned
    │   │   │   │   ├── 214704-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214721-consensus-masked.coords
    │   │   │   │   ├── 214721-consensus-masked.coords.filtered
    │   │   │   │   ├── 214721-consensus-masked.coords_tmp
    │   │   │   │   ├── 214721-consensus-masked.sf
    │   │   │   │   ├── 214721-consensus-masked.unaligned
    │   │   │   │   ├── 214721-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214724-consensus-masked.coords
    │   │   │   │   ├── 214724-consensus-masked.coords.filtered
    │   │   │   │   ├── 214724-consensus-masked.coords_tmp
    │   │   │   │   ├── 214724-consensus-masked.sf
    │   │   │   │   ├── 214724-consensus-masked.unaligned
    │   │   │   │   ├── 214724-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214760-consensus-masked.coords
    │   │   │   │   ├── 214760-consensus-masked.coords.filtered
    │   │   │   │   ├── 214760-consensus-masked.coords_tmp
    │   │   │   │   ├── 214760-consensus-masked.sf
    │   │   │   │   ├── 214760-consensus-masked.unaligned
    │   │   │   │   ├── 214760-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214763-consensus-masked.coords
    │   │   │   │   ├── 214763-consensus-masked.coords.filtered
    │   │   │   │   ├── 214763-consensus-masked.coords_tmp
    │   │   │   │   ├── 214763-consensus-masked.sf
    │   │   │   │   ├── 214763-consensus-masked.unaligned
    │   │   │   │   ├── 214763-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214765-consensus-masked.coords
    │   │   │   │   ├── 214765-consensus-masked.coords.filtered
    │   │   │   │   ├── 214765-consensus-masked.coords_tmp
    │   │   │   │   ├── 214765-consensus-masked.sf
    │   │   │   │   ├── 214765-consensus-masked.unaligned
    │   │   │   │   ├── 214765-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214766-consensus-masked.coords
    │   │   │   │   ├── 214766-consensus-masked.coords.filtered
    │   │   │   │   ├── 214766-consensus-masked.coords_tmp
    │   │   │   │   ├── 214766-consensus-masked.sf
    │   │   │   │   ├── 214766-consensus-masked.unaligned
    │   │   │   │   ├── 214766-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214780-consensus-masked.coords
    │   │   │   │   ├── 214780-consensus-masked.coords.filtered
    │   │   │   │   ├── 214780-consensus-masked.coords_tmp
    │   │   │   │   ├── 214780-consensus-masked.sf
    │   │   │   │   ├── 214780-consensus-masked.unaligned
    │   │   │   │   ├── 214780-consensus-masked.used_snps.gz
    │   │   │   │   ├── 214807-consensus-masked.coords
    │   │   │   │   ├── 214807-consensus-masked.coords.filtered
    │   │   │   │   ├── 214807-consensus-masked.coords_tmp
    │   │   │   │   ├── 214807-consensus-masked.sf
    │   │   │   │   ├── 214807-consensus-masked.unaligned
    │   │   │   │   └── 214807-consensus-masked.used_snps.gz
    │   │   │   ├── misassemblies_frcurve_plot.pdf
    │   │   │   ├── misassemblies_plot.pdf
    │   │   │   ├── misassemblies_report.tex
    │   │   │   ├── misassemblies_report.tsv
    │   │   │   ├── misassemblies_report.txt
    │   │   │   ├── transposed_report_misassemblies.tex
    │   │   │   ├── transposed_report_misassemblies.tsv
    │   │   │   ├── transposed_report_misassemblies.txt
    │   │   │   ├── unaligned_report.tex
    │   │   │   ├── unaligned_report.tsv
    │   │   │   └── unaligned_report.txt
    │   │   ├── genome_stats
    │   │   │   ├── 214704-consensus-masked_gaps.txt
    │   │   │   ├── 214704-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214721-consensus-masked_gaps.txt
    │   │   │   ├── 214721-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214724-consensus-masked_gaps.txt
    │   │   │   ├── 214724-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214760-consensus-masked_gaps.txt
    │   │   │   ├── 214760-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214763-consensus-masked_gaps.txt
    │   │   │   ├── 214763-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214765-consensus-masked_gaps.txt
    │   │   │   ├── 214765-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214766-consensus-masked_gaps.txt
    │   │   │   ├── 214766-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214780-consensus-masked_gaps.txt
    │   │   │   ├── 214780-consensus-masked_genomic_features_any.txt
    │   │   │   ├── 214807-consensus-masked_gaps.txt
    │   │   │   ├── 214807-consensus-masked_genomic_features_any.txt
    │   │   │   ├── complete_features_histogram.pdf
    │   │   │   ├── features_cumulative_plot.pdf
    │   │   │   ├── features_frcurve_plot.pdf
    │   │   │   ├── genome_fraction_histogram.pdf
    │   │   │   └── genome_info.txt
    │   │   ├── icarus.html
    │   │   ├── icarus_viewers
    │   │   │   ├── alignment_viewer.html
    │   │   │   └── contig_size_viewer.html
    │   │   ├── quast.log
    │   │   ├── report.html
    │   │   ├── report.pdf
    │   │   ├── report.tex
    │   │   ├── report.tsv
    │   │   ├── report.txt
    │   │   ├── transposed_report.tex
    │   │   ├── transposed_report.tsv
    │   │   └── transposed_report.txt
    │   └── snpeff
    │       ├── 214704.snpEff.csv
    │       ├── 214704.snpEff.genes.txt
    │       ├── 214704.snpEff.summary.html
    │       ├── 214704.snpEff.vcf.gz
    │       ├── 214704.snpEff.vcf.gz.tbi
    │       ├── 214704.snpSift.table.txt
    │       ├── 214721.snpEff.csv
    │       ├── 214721.snpEff.genes.txt
    │       ├── 214721.snpEff.summary.html
    │       ├── 214721.snpEff.vcf.gz
    │       ├── 214721.snpEff.vcf.gz.tbi
    │       ├── 214721.snpSift.table.txt
    │       ├── 214724.snpEff.csv
    │       ├── 214724.snpEff.genes.txt
    │       ├── 214724.snpEff.summary.html
    │       ├── 214724.snpEff.vcf.gz
    │       ├── 214724.snpEff.vcf.gz.tbi
    │       ├── 214724.snpSift.table.txt
    │       ├── 214760.snpEff.csv
    │       ├── 214760.snpEff.genes.txt
    │       ├── 214760.snpEff.summary.html
    │       ├── 214760.snpEff.vcf.gz
    │       ├── 214760.snpEff.vcf.gz.tbi
    │       ├── 214760.snpSift.table.txt
    │       ├── 214763.snpEff.csv
    │       ├── 214763.snpEff.genes.txt
    │       ├── 214763.snpEff.summary.html
    │       ├── 214763.snpEff.vcf.gz
    │       ├── 214763.snpEff.vcf.gz.tbi
    │       ├── 214763.snpSift.table.txt
    │       ├── 214765.snpEff.csv
    │       ├── 214765.snpEff.genes.txt
    │       ├── 214765.snpEff.summary.html
    │       ├── 214765.snpEff.vcf.gz
    │       ├── 214765.snpEff.vcf.gz.tbi
    │       ├── 214765.snpSift.table.txt
    │       ├── 214766.snpEff.csv
    │       ├── 214766.snpEff.genes.txt
    │       ├── 214766.snpEff.summary.html
    │       ├── 214766.snpEff.vcf.gz
    │       ├── 214766.snpEff.vcf.gz.tbi
    │       ├── 214766.snpSift.table.txt
    │       ├── 214780.snpEff.csv
    │       ├── 214780.snpEff.genes.txt
    │       ├── 214780.snpEff.summary.html
    │       ├── 214780.snpEff.vcf.gz
    │       ├── 214780.snpEff.vcf.gz.tbi
    │       ├── 214780.snpSift.table.txt
    │       ├── 214807.snpEff.csv
    │       ├── 214807.snpEff.genes.txt
    │       ├── 214807.snpEff.summary.html
    │       ├── 214807.snpEff.vcf.gz
    │       ├── 214807.snpEff.vcf.gz.tbi
    │       └── 214807.snpSift.table.txt
    ├── intersect
    │   ├── 214704
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214721
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214724
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214760
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214763
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214765
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214766
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   ├── 214780
    │   │   ├── 0000.vcf.gz
    │   │   ├── 0000.vcf.gz.tbi
    │   │   ├── 0001.vcf.gz
    │   │   ├── 0001.vcf.gz.tbi
    │   │   ├── 0002.vcf.gz
    │   │   ├── 0002.vcf.gz.tbi
    │   │   ├── README.txt
    │   │   └── sites.txt
    │   └── 214807
    │       ├── 0000.vcf.gz
    │       ├── 0000.vcf.gz.tbi
    │       ├── 0001.vcf.gz
    │       ├── 0001.vcf.gz.tbi
    │       ├── 0002.vcf.gz
    │       ├── 0002.vcf.gz.tbi
    │       ├── README.txt
    │       └── sites.txt
    ├── ivar
    │   ├── 214704.AF0.75.vcf.gz
    │   ├── 214704.AF0.75.vcf.gz.tbi
    │   ├── 214704.tsv
    │   ├── 214704.vcf.gz
    │   ├── 214704.vcf.gz.tbi
    │   ├── 214721.AF0.75.vcf.gz
    │   ├── 214721.AF0.75.vcf.gz.tbi
    │   ├── 214721.tsv
    │   ├── 214721.vcf.gz
    │   ├── 214721.vcf.gz.tbi
    │   ├── 214724.AF0.75.vcf.gz
    │   ├── 214724.AF0.75.vcf.gz.tbi
    │   ├── 214724.tsv
    │   ├── 214724.vcf.gz
    │   ├── 214724.vcf.gz.tbi
    │   ├── 214760.AF0.75.vcf.gz
    │   ├── 214760.AF0.75.vcf.gz.tbi
    │   ├── 214760.tsv
    │   ├── 214760.vcf.gz
    │   ├── 214760.vcf.gz.tbi
    │   ├── 214763.AF0.75.vcf.gz
    │   ├── 214763.AF0.75.vcf.gz.tbi
    │   ├── 214763.tsv
    │   ├── 214763.vcf.gz
    │   ├── 214763.vcf.gz.tbi
    │   ├── 214765.AF0.75.vcf.gz
    │   ├── 214765.AF0.75.vcf.gz.tbi
    │   ├── 214765.tsv
    │   ├── 214765.vcf.gz
    │   ├── 214765.vcf.gz.tbi
    │   ├── 214766.AF0.75.vcf.gz
    │   ├── 214766.AF0.75.vcf.gz.tbi
    │   ├── 214766.tsv
    │   ├── 214766.vcf.gz
    │   ├── 214766.vcf.gz.tbi
    │   ├── 214780.AF0.75.vcf.gz
    │   ├── 214780.AF0.75.vcf.gz.tbi
    │   ├── 214780.tsv
    │   ├── 214780.vcf.gz
    │   ├── 214780.vcf.gz.tbi
    │   ├── 214807.AF0.75.vcf.gz
    │   ├── 214807.AF0.75.vcf.gz.tbi
    │   ├── 214807.tsv
    │   ├── 214807.vcf.gz
    │   ├── 214807.vcf.gz.tbi
    │   ├── bcftools_stats
    │   │   ├── 214704.AF0.75.bcftools_stats.txt
    │   │   ├── 214704.bcftools_stats.txt
    │   │   ├── 214721.AF0.75.bcftools_stats.txt
    │   │   ├── 214721.bcftools_stats.txt
    │   │   ├── 214724.AF0.75.bcftools_stats.txt
    │   │   ├── 214724.bcftools_stats.txt
    │   │   ├── 214760.AF0.75.bcftools_stats.txt
    │   │   ├── 214760.bcftools_stats.txt
    │   │   ├── 214763.AF0.75.bcftools_stats.txt
    │   │   ├── 214763.bcftools_stats.txt
    │   │   ├── 214765.AF0.75.bcftools_stats.txt
    │   │   ├── 214765.bcftools_stats.txt
    │   │   ├── 214766.AF0.75.bcftools_stats.txt
    │   │   ├── 214766.bcftools_stats.txt
    │   │   ├── 214780.AF0.75.bcftools_stats.txt
    │   │   ├── 214780.bcftools_stats.txt
    │   │   ├── 214807.AF0.75.bcftools_stats.txt
    │   │   └── 214807.bcftools_stats.txt
    │   ├── consensus
    │   │   ├── 214704.AF0.75.consensus.fa
    │   │   ├── 214704.AF0.75.consensus.qual.txt
    │   │   ├── 214721.AF0.75.consensus.fa
    │   │   ├── 214721.AF0.75.consensus.qual.txt
    │   │   ├── 214724.AF0.75.consensus.fa
    │   │   ├── 214724.AF0.75.consensus.qual.txt
    │   │   ├── 214760.AF0.75.consensus.fa
    │   │   ├── 214760.AF0.75.consensus.qual.txt
    │   │   ├── 214763.AF0.75.consensus.fa
    │   │   ├── 214763.AF0.75.consensus.qual.txt
    │   │   ├── 214765.AF0.75.consensus.fa
    │   │   ├── 214765.AF0.75.consensus.qual.txt
    │   │   ├── 214766.AF0.75.consensus.fa
    │   │   ├── 214766.AF0.75.consensus.qual.txt
    │   │   ├── 214780.AF0.75.consensus.fa
    │   │   ├── 214780.AF0.75.consensus.qual.txt
    │   │   ├── 214807.AF0.75.consensus.fa
    │   │   ├── 214807.AF0.75.consensus.qual.txt
    │   │   └── base_qc
    │   │       ├── 214704.AF0.75.ACTG_density.pdf
    │   │       ├── 214704.AF0.75.base_counts.pdf
    │   │       ├── 214704.AF0.75.base_counts.tsv
    │   │       ├── 214704.AF0.75.N_density.pdf
    │   │       ├── 214704.AF0.75.N_run.tsv
    │   │       ├── 214704.AF0.75.R_density.pdf
    │   │       ├── 214704.AF0.75.Y_density.pdf
    │   │       ├── 214721.AF0.75.ACTG_density.pdf
    │   │       ├── 214721.AF0.75.base_counts.pdf
    │   │       ├── 214721.AF0.75.base_counts.tsv
    │   │       ├── 214721.AF0.75.N_density.pdf
    │   │       ├── 214721.AF0.75.N_run.tsv
    │   │       ├── 214721.AF0.75.W_density.pdf
    │   │       ├── 214721.AF0.75.Y_density.pdf
    │   │       ├── 214724.AF0.75.ACTG_density.pdf
    │   │       ├── 214724.AF0.75.base_counts.pdf
    │   │       ├── 214724.AF0.75.base_counts.tsv
    │   │       ├── 214724.AF0.75.M_density.pdf
    │   │       ├── 214724.AF0.75.N_density.pdf
    │   │       ├── 214724.AF0.75.N_run.tsv
    │   │       ├── 214724.AF0.75.R_density.pdf
    │   │       ├── 214760.AF0.75.ACTG_density.pdf
    │   │       ├── 214760.AF0.75.base_counts.pdf
    │   │       ├── 214760.AF0.75.base_counts.tsv
    │   │       ├── 214760.AF0.75.N_density.pdf
    │   │       ├── 214760.AF0.75.N_run.tsv
    │   │       ├── 214760.AF0.75.R_density.pdf
    │   │       ├── 214760.AF0.75.W_density.pdf
    │   │       ├── 214760.AF0.75.Y_density.pdf
    │   │       ├── 214763.AF0.75.ACTG_density.pdf
    │   │       ├── 214763.AF0.75.base_counts.pdf
    │   │       ├── 214763.AF0.75.base_counts.tsv
    │   │       ├── 214763.AF0.75.N_density.pdf
    │   │       ├── 214763.AF0.75.N_run.tsv
    │   │       ├── 214763.AF0.75.R_density.pdf
    │   │       ├── 214763.AF0.75.Y_density.pdf
    │   │       ├── 214765.AF0.75.ACTG_density.pdf
    │   │       ├── 214765.AF0.75.base_counts.pdf
    │   │       ├── 214765.AF0.75.base_counts.tsv
    │   │       ├── 214765.AF0.75.N_density.pdf
    │   │       ├── 214765.AF0.75.N_run.tsv
    │   │       ├── 214765.AF0.75.R_density.pdf
    │   │       ├── 214766.AF0.75.ACTG_density.pdf
    │   │       ├── 214766.AF0.75.base_counts.pdf
    │   │       ├── 214766.AF0.75.base_counts.tsv
    │   │       ├── 214766.AF0.75.N_density.pdf
    │   │       ├── 214766.AF0.75.N_run.tsv
    │   │       ├── 214766.AF0.75.R_density.pdf
    │   │       ├── 214780.AF0.75.ACTG_density.pdf
    │   │       ├── 214780.AF0.75.base_counts.pdf
    │   │       ├── 214780.AF0.75.base_counts.tsv
    │   │       ├── 214780.AF0.75.N_density.pdf
    │   │       ├── 214780.AF0.75.N_run.tsv
    │   │       ├── 214780.AF0.75.R_density.pdf
    │   │       ├── 214807.AF0.75.ACTG_density.pdf
    │   │       ├── 214807.AF0.75.base_counts.pdf
    │   │       ├── 214807.AF0.75.base_counts.tsv
    │   │       ├── 214807.AF0.75.M_density.pdf
    │   │       ├── 214807.AF0.75.N_density.pdf
    │   │       ├── 214807.AF0.75.N_run.tsv
    │   │       ├── 214807.AF0.75.R_density.pdf
    │   │       └── 214807.AF0.75.W_density.pdf
    │   ├── log
    │   │   ├── 214704.AF0.75.variant.counts.log
    │   │   ├── 214704.variant.counts.log
    │   │   ├── 214721.AF0.75.variant.counts.log
    │   │   ├── 214721.variant.counts.log
    │   │   ├── 214724.AF0.75.variant.counts.log
    │   │   ├── 214724.variant.counts.log
    │   │   ├── 214760.AF0.75.variant.counts.log
    │   │   ├── 214760.variant.counts.log
    │   │   ├── 214763.AF0.75.variant.counts.log
    │   │   ├── 214763.variant.counts.log
    │   │   ├── 214765.AF0.75.variant.counts.log
    │   │   ├── 214765.variant.counts.log
    │   │   ├── 214766.AF0.75.variant.counts.log
    │   │   ├── 214766.variant.counts.log
    │   │   ├── 214780.AF0.75.variant.counts.log
    │   │   ├── 214780.variant.counts.log
    │   │   ├── 214807.AF0.75.variant.counts.log
    │   │   └── 214807.variant.counts.log
    │   ├── quast
    │   │   └── AF0.75
    │   │       ├── aligned_stats
    │   │       │   ├── cumulative_plot.pdf
    │   │       │   ├── NAx_plot.pdf
    │   │       │   └── NGAx_plot.pdf
    │   │       ├── basic_stats
    │   │       │   ├── 214704.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214721.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214724.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214760.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214763.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214765.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214766.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214780.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── 214807.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── cumulative_plot.pdf
    │   │       │   ├── GC_content_plot.pdf
    │   │       │   ├── gc.icarus.txt
    │   │       │   ├── NGx_plot.pdf
    │   │       │   └── Nx_plot.pdf
    │   │       ├── contigs_reports
    │   │       │   ├── 214704_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214721_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214724_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214760_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214763_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214765_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214766_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214780_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── 214807_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── all_alignments_214704-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214721-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214724-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214760-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214763-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214765-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214766-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214780-AF0-75-consensus.tsv
    │   │       │   ├── all_alignments_214807-AF0-75-consensus.tsv
    │   │       │   ├── contigs_report_214704-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214704-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214704-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214704-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214721-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214721-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214721-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214721-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214724-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214724-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214724-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214724-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214760-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214760-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214760-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214760-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214763-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214763-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214763-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214763-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214765-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214765-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214765-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214765-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214766-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214766-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214766-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214766-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214780-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214780-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214780-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214780-AF0-75-consensus.unaligned.info
    │   │       │   ├── contigs_report_214807-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214807-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214807-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214807-AF0-75-consensus.unaligned.info
    │   │       │   ├── minimap_output
    │   │       │   │   ├── 214704-AF0-75-consensus.coords
    │   │       │   │   ├── 214704-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214704-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214704-AF0-75-consensus.sf
    │   │       │   │   ├── 214704-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214704-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214721-AF0-75-consensus.coords
    │   │       │   │   ├── 214721-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214721-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214721-AF0-75-consensus.sf
    │   │       │   │   ├── 214721-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214721-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214724-AF0-75-consensus.coords
    │   │       │   │   ├── 214724-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214724-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214724-AF0-75-consensus.sf
    │   │       │   │   ├── 214724-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214724-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214760-AF0-75-consensus.coords
    │   │       │   │   ├── 214760-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214760-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214760-AF0-75-consensus.sf
    │   │       │   │   ├── 214760-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214760-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214763-AF0-75-consensus.coords
    │   │       │   │   ├── 214763-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214763-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214763-AF0-75-consensus.sf
    │   │       │   │   ├── 214763-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214763-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214765-AF0-75-consensus.coords
    │   │       │   │   ├── 214765-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214765-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214765-AF0-75-consensus.sf
    │   │       │   │   ├── 214765-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214765-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214766-AF0-75-consensus.coords
    │   │       │   │   ├── 214766-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214766-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214766-AF0-75-consensus.sf
    │   │       │   │   ├── 214766-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214766-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214780-AF0-75-consensus.coords
    │   │       │   │   ├── 214780-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214780-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214780-AF0-75-consensus.sf
    │   │       │   │   ├── 214780-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214780-AF0-75-consensus.used_snps.gz
    │   │       │   │   ├── 214807-AF0-75-consensus.coords
    │   │       │   │   ├── 214807-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214807-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214807-AF0-75-consensus.sf
    │   │       │   │   ├── 214807-AF0-75-consensus.unaligned
    │   │       │   │   └── 214807-AF0-75-consensus.used_snps.gz
    │   │       │   ├── misassemblies_frcurve_plot.pdf
    │   │       │   ├── misassemblies_plot.pdf
    │   │       │   ├── misassemblies_report.tex
    │   │       │   ├── misassemblies_report.tsv
    │   │       │   ├── misassemblies_report.txt
    │   │       │   ├── transposed_report_misassemblies.tex
    │   │       │   ├── transposed_report_misassemblies.tsv
    │   │       │   ├── transposed_report_misassemblies.txt
    │   │       │   ├── unaligned_report.tex
    │   │       │   ├── unaligned_report.tsv
    │   │       │   └── unaligned_report.txt
    │   │       ├── genome_stats
    │   │       │   ├── 214704-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214704-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214721-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214721-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214724-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214724-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214760-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214760-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214763-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214763-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214765-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214765-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214766-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214766-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214780-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214780-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── 214807-AF0-75-consensus_gaps.txt
    │   │       │   ├── 214807-AF0-75-consensus_genomic_features_any.txt
    │   │       │   ├── complete_features_histogram.pdf
    │   │       │   ├── features_cumulative_plot.pdf
    │   │       │   ├── features_frcurve_plot.pdf
    │   │       │   ├── genome_fraction_histogram.pdf
    │   │       │   └── genome_info.txt
    │   │       ├── icarus.html
    │   │       ├── icarus_viewers
    │   │       │   ├── alignment_viewer.html
    │   │       │   └── contig_size_viewer.html
    │   │       ├── quast.log
    │   │       ├── report.html
    │   │       ├── report.pdf
    │   │       ├── report.tex
    │   │       ├── report.tsv
    │   │       ├── report.txt
    │   │       ├── transposed_report.tex
    │   │       ├── transposed_report.tsv
    │   │       └── transposed_report.txt
    │   └── snpeff
    │       ├── 214704.AF0.75.snpEff.csv
    │       ├── 214704.AF0.75.snpEff.genes.txt
    │       ├── 214704.AF0.75.snpEff.summary.html
    │       ├── 214704.AF0.75.snpEff.vcf.gz
    │       ├── 214704.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214704.AF0.75.snpSift.table.txt
    │       ├── 214704.snpEff.csv
    │       ├── 214704.snpEff.genes.txt
    │       ├── 214704.snpEff.summary.html
    │       ├── 214704.snpEff.vcf.gz
    │       ├── 214704.snpEff.vcf.gz.tbi
    │       ├── 214704.snpSift.table.txt
    │       ├── 214721.AF0.75.snpEff.csv
    │       ├── 214721.AF0.75.snpEff.genes.txt
    │       ├── 214721.AF0.75.snpEff.summary.html
    │       ├── 214721.AF0.75.snpEff.vcf.gz
    │       ├── 214721.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214721.AF0.75.snpSift.table.txt
    │       ├── 214721.snpEff.csv
    │       ├── 214721.snpEff.genes.txt
    │       ├── 214721.snpEff.summary.html
    │       ├── 214721.snpEff.vcf.gz
    │       ├── 214721.snpEff.vcf.gz.tbi
    │       ├── 214721.snpSift.table.txt
    │       ├── 214724.AF0.75.snpEff.csv
    │       ├── 214724.AF0.75.snpEff.genes.txt
    │       ├── 214724.AF0.75.snpEff.summary.html
    │       ├── 214724.AF0.75.snpEff.vcf.gz
    │       ├── 214724.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214724.AF0.75.snpSift.table.txt
    │       ├── 214724.snpEff.csv
    │       ├── 214724.snpEff.genes.txt
    │       ├── 214724.snpEff.summary.html
    │       ├── 214724.snpEff.vcf.gz
    │       ├── 214724.snpEff.vcf.gz.tbi
    │       ├── 214724.snpSift.table.txt
    │       ├── 214760.AF0.75.snpEff.csv
    │       ├── 214760.AF0.75.snpEff.genes.txt
    │       ├── 214760.AF0.75.snpEff.summary.html
    │       ├── 214760.AF0.75.snpEff.vcf.gz
    │       ├── 214760.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214760.AF0.75.snpSift.table.txt
    │       ├── 214760.snpEff.csv
    │       ├── 214760.snpEff.genes.txt
    │       ├── 214760.snpEff.summary.html
    │       ├── 214760.snpEff.vcf.gz
    │       ├── 214760.snpEff.vcf.gz.tbi
    │       ├── 214760.snpSift.table.txt
    │       ├── 214763.AF0.75.snpEff.csv
    │       ├── 214763.AF0.75.snpEff.genes.txt
    │       ├── 214763.AF0.75.snpEff.summary.html
    │       ├── 214763.AF0.75.snpEff.vcf.gz
    │       ├── 214763.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214763.AF0.75.snpSift.table.txt
    │       ├── 214763.snpEff.csv
    │       ├── 214763.snpEff.genes.txt
    │       ├── 214763.snpEff.summary.html
    │       ├── 214763.snpEff.vcf.gz
    │       ├── 214763.snpEff.vcf.gz.tbi
    │       ├── 214763.snpSift.table.txt
    │       ├── 214765.AF0.75.snpEff.csv
    │       ├── 214765.AF0.75.snpEff.genes.txt
    │       ├── 214765.AF0.75.snpEff.summary.html
    │       ├── 214765.AF0.75.snpEff.vcf.gz
    │       ├── 214765.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214765.AF0.75.snpSift.table.txt
    │       ├── 214765.snpEff.csv
    │       ├── 214765.snpEff.genes.txt
    │       ├── 214765.snpEff.summary.html
    │       ├── 214765.snpEff.vcf.gz
    │       ├── 214765.snpEff.vcf.gz.tbi
    │       ├── 214765.snpSift.table.txt
    │       ├── 214766.AF0.75.snpEff.csv
    │       ├── 214766.AF0.75.snpEff.genes.txt
    │       ├── 214766.AF0.75.snpEff.summary.html
    │       ├── 214766.AF0.75.snpEff.vcf.gz
    │       ├── 214766.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214766.AF0.75.snpSift.table.txt
    │       ├── 214766.snpEff.csv
    │       ├── 214766.snpEff.genes.txt
    │       ├── 214766.snpEff.summary.html
    │       ├── 214766.snpEff.vcf.gz
    │       ├── 214766.snpEff.vcf.gz.tbi
    │       ├── 214766.snpSift.table.txt
    │       ├── 214780.AF0.75.snpEff.csv
    │       ├── 214780.AF0.75.snpEff.genes.txt
    │       ├── 214780.AF0.75.snpEff.summary.html
    │       ├── 214780.AF0.75.snpEff.vcf.gz
    │       ├── 214780.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214780.AF0.75.snpSift.table.txt
    │       ├── 214780.snpEff.csv
    │       ├── 214780.snpEff.genes.txt
    │       ├── 214780.snpEff.summary.html
    │       ├── 214780.snpEff.vcf.gz
    │       ├── 214780.snpEff.vcf.gz.tbi
    │       ├── 214780.snpSift.table.txt
    │       ├── 214807.AF0.75.snpEff.csv
    │       ├── 214807.AF0.75.snpEff.genes.txt
    │       ├── 214807.AF0.75.snpEff.summary.html
    │       ├── 214807.AF0.75.snpEff.vcf.gz
    │       ├── 214807.AF0.75.snpEff.vcf.gz.tbi
    │       ├── 214807.AF0.75.snpSift.table.txt
    │       ├── 214807.snpEff.csv
    │       ├── 214807.snpEff.genes.txt
    │       ├── 214807.snpEff.summary.html
    │       ├── 214807.snpEff.vcf.gz
    │       ├── 214807.snpEff.vcf.gz.tbi
    │       ├── 214807.snpSift.table.txt
    │       ├── snpSift_template2.txt
    │       ├── snpSift_template_filtered.txt
    │       └── snpSift_template.txt
    ├── summary_variants_metrics_mqc.tsv
    └── varscan2
        ├── 214704.AF0.75.vcf.gz
        ├── 214704.AF0.75.vcf.gz.tbi
        ├── 214704.vcf.gz
        ├── 214704.vcf.gz.tbi
        ├── 214721.AF0.75.vcf.gz
        ├── 214721.AF0.75.vcf.gz.tbi
        ├── 214721.vcf.gz
        ├── 214721.vcf.gz.tbi
        ├── 214724.AF0.75.vcf.gz
        ├── 214724.AF0.75.vcf.gz.tbi
        ├── 214724.vcf.gz
        ├── 214724.vcf.gz.tbi
        ├── 214760.AF0.75.vcf.gz
        ├── 214760.AF0.75.vcf.gz.tbi
        ├── 214760.vcf.gz
        ├── 214760.vcf.gz.tbi
        ├── 214763.AF0.75.vcf.gz
        ├── 214763.AF0.75.vcf.gz.tbi
        ├── 214763.vcf.gz
        ├── 214763.vcf.gz.tbi
        ├── 214765.AF0.75.vcf.gz
        ├── 214765.AF0.75.vcf.gz.tbi
        ├── 214765.vcf.gz
        ├── 214765.vcf.gz.tbi
        ├── 214766.AF0.75.vcf.gz
        ├── 214766.AF0.75.vcf.gz.tbi
        ├── 214766.vcf.gz
        ├── 214766.vcf.gz.tbi
        ├── 214780.AF0.75.vcf.gz
        ├── 214780.AF0.75.vcf.gz.tbi
        ├── 214780.vcf.gz
        ├── 214780.vcf.gz.tbi
        ├── 214807.AF0.75.vcf.gz
        ├── 214807.AF0.75.vcf.gz.tbi
        ├── 214807.vcf.gz
        ├── 214807.vcf.gz.tbi
        ├── bcftools_stats
        │   ├── 214704.AF0.75.bcftools_stats.txt
        │   ├── 214704.bcftools_stats.txt
        │   ├── 214721.AF0.75.bcftools_stats.txt
        │   ├── 214721.bcftools_stats.txt
        │   ├── 214724.AF0.75.bcftools_stats.txt
        │   ├── 214724.bcftools_stats.txt
        │   ├── 214760.AF0.75.bcftools_stats.txt
        │   ├── 214760.bcftools_stats.txt
        │   ├── 214763.AF0.75.bcftools_stats.txt
        │   ├── 214763.bcftools_stats.txt
        │   ├── 214765.AF0.75.bcftools_stats.txt
        │   ├── 214765.bcftools_stats.txt
        │   ├── 214766.AF0.75.bcftools_stats.txt
        │   ├── 214766.bcftools_stats.txt
        │   ├── 214780.AF0.75.bcftools_stats.txt
        │   ├── 214780.bcftools_stats.txt
        │   ├── 214807.AF0.75.bcftools_stats.txt
        │   └── 214807.bcftools_stats.txt
        ├── consensus
        │   ├── 214704.AF0.75.consensus.masked.fa
        │   ├── 214721.AF0.75.consensus.masked.fa
        │   ├── 214724.AF0.75.consensus.masked.fa
        │   ├── 214760.AF0.75.consensus.masked.fa
        │   ├── 214763.AF0.75.consensus.masked.fa
        │   ├── 214765.AF0.75.consensus.masked.fa
        │   ├── 214766.AF0.75.consensus.masked.fa
        │   ├── 214780.AF0.75.consensus.masked.fa
        │   ├── 214807.AF0.75.consensus.masked.fa
        │   ├── base_qc
        │   │   ├── 214704.AF0.75.ACTG_density.pdf
        │   │   ├── 214704.AF0.75.base_counts.pdf
        │   │   ├── 214704.AF0.75.base_counts.tsv
        │   │   ├── 214704.AF0.75.N_density.pdf
        │   │   ├── 214704.AF0.75.N_run.tsv
        │   │   ├── 214721.AF0.75.ACTG_density.pdf
        │   │   ├── 214721.AF0.75.base_counts.pdf
        │   │   ├── 214721.AF0.75.base_counts.tsv
        │   │   ├── 214721.AF0.75.N_density.pdf
        │   │   ├── 214721.AF0.75.N_run.tsv
        │   │   ├── 214724.AF0.75.ACTG_density.pdf
        │   │   ├── 214724.AF0.75.base_counts.pdf
        │   │   ├── 214724.AF0.75.base_counts.tsv
        │   │   ├── 214724.AF0.75.N_density.pdf
        │   │   ├── 214724.AF0.75.N_run.tsv
        │   │   ├── 214760.AF0.75.ACTG_density.pdf
        │   │   ├── 214760.AF0.75.base_counts.pdf
        │   │   ├── 214760.AF0.75.base_counts.tsv
        │   │   ├── 214760.AF0.75.N_density.pdf
        │   │   ├── 214760.AF0.75.N_run.tsv
        │   │   ├── 214763.AF0.75.ACTG_density.pdf
        │   │   ├── 214763.AF0.75.base_counts.pdf
        │   │   ├── 214763.AF0.75.base_counts.tsv
        │   │   ├── 214763.AF0.75.N_density.pdf
        │   │   ├── 214763.AF0.75.N_run.tsv
        │   │   ├── 214765.AF0.75.ACTG_density.pdf
        │   │   ├── 214765.AF0.75.base_counts.pdf
        │   │   ├── 214765.AF0.75.base_counts.tsv
        │   │   ├── 214765.AF0.75.N_density.pdf
        │   │   ├── 214765.AF0.75.N_run.tsv
        │   │   ├── 214766.AF0.75.ACTG_density.pdf
        │   │   ├── 214766.AF0.75.base_counts.pdf
        │   │   ├── 214766.AF0.75.base_counts.tsv
        │   │   ├── 214766.AF0.75.N_density.pdf
        │   │   ├── 214766.AF0.75.N_run.tsv
        │   │   ├── 214780.AF0.75.ACTG_density.pdf
        │   │   ├── 214780.AF0.75.base_counts.pdf
        │   │   ├── 214780.AF0.75.base_counts.tsv
        │   │   ├── 214780.AF0.75.N_density.pdf
        │   │   ├── 214780.AF0.75.N_run.tsv
        │   │   ├── 214807.AF0.75.ACTG_density.pdf
        │   │   ├── 214807.AF0.75.base_counts.pdf
        │   │   ├── 214807.AF0.75.base_counts.tsv
        │   │   ├── 214807.AF0.75.N_density.pdf
        │   │   └── 214807.AF0.75.N_run.tsv
        │   └── snpSift_template.txt
        ├── log
        │   ├── 214704.varscan2.log
        │   ├── 214721.varscan2.log
        │   ├── 214724.varscan2.log
        │   ├── 214760.varscan2.log
        │   ├── 214763.varscan2.log
        │   ├── 214765.varscan2.log
        │   ├── 214766.varscan2.log
        │   ├── 214780.varscan2.log
        │   └── 214807.varscan2.log
        ├── quast
        │   └── AF0.75
        │       ├── aligned_stats
        │       │   ├── cumulative_plot.pdf
        │       │   ├── NAx_plot.pdf
        │       │   └── NGAx_plot.pdf
        │       ├── basic_stats
        │       │   ├── 214704.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214721.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214724.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214760.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214763.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214765.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214766.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214780.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── 214807.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── cumulative_plot.pdf
        │       │   ├── GC_content_plot.pdf
        │       │   ├── gc.icarus.txt
        │       │   ├── NGx_plot.pdf
        │       │   └── Nx_plot.pdf
        │       ├── contigs_reports
        │       │   ├── 214704_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214721_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214724_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214760_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214763_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214765_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214766_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214780_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── 214807_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── all_alignments_214704-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214721-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214724-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214760-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214763-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214765-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214766-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214780-AF0-75-consensus-masked.tsv
        │       │   ├── all_alignments_214807-AF0-75-consensus-masked.tsv
        │       │   ├── contigs_report_214704-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214704-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214704-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214704-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214721-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214721-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214721-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214721-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214724-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214724-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214724-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214724-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214760-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214760-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214760-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214760-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214763-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214763-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214763-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214763-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214765-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214765-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214765-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214765-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214766-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214766-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214766-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214766-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214780-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214780-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214780-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214780-AF0-75-consensus-masked.unaligned.info
        │       │   ├── contigs_report_214807-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214807-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214807-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214807-AF0-75-consensus-masked.unaligned.info
        │       │   ├── minimap_output
        │       │   │   ├── 214704-AF0-75-consensus-masked.coords
        │       │   │   ├── 214704-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214704-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214704-AF0-75-consensus-masked.sf
        │       │   │   ├── 214704-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214704-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214721-AF0-75-consensus-masked.coords
        │       │   │   ├── 214721-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214721-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214721-AF0-75-consensus-masked.sf
        │       │   │   ├── 214721-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214721-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214724-AF0-75-consensus-masked.coords
        │       │   │   ├── 214724-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214724-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214724-AF0-75-consensus-masked.sf
        │       │   │   ├── 214724-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214724-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214760-AF0-75-consensus-masked.coords
        │       │   │   ├── 214760-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214760-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214760-AF0-75-consensus-masked.sf
        │       │   │   ├── 214760-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214760-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214763-AF0-75-consensus-masked.coords
        │       │   │   ├── 214763-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214763-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214763-AF0-75-consensus-masked.sf
        │       │   │   ├── 214763-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214763-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214765-AF0-75-consensus-masked.coords
        │       │   │   ├── 214765-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214765-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214765-AF0-75-consensus-masked.sf
        │       │   │   ├── 214765-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214765-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214766-AF0-75-consensus-masked.coords
        │       │   │   ├── 214766-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214766-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214766-AF0-75-consensus-masked.sf
        │       │   │   ├── 214766-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214766-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214780-AF0-75-consensus-masked.coords
        │       │   │   ├── 214780-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214780-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214780-AF0-75-consensus-masked.sf
        │       │   │   ├── 214780-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214780-AF0-75-consensus-masked.used_snps.gz
        │       │   │   ├── 214807-AF0-75-consensus-masked.coords
        │       │   │   ├── 214807-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214807-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214807-AF0-75-consensus-masked.sf
        │       │   │   ├── 214807-AF0-75-consensus-masked.unaligned
        │       │   │   └── 214807-AF0-75-consensus-masked.used_snps.gz
        │       │   ├── misassemblies_frcurve_plot.pdf
        │       │   ├── misassemblies_plot.pdf
        │       │   ├── misassemblies_report.tex
        │       │   ├── misassemblies_report.tsv
        │       │   ├── misassemblies_report.txt
        │       │   ├── transposed_report_misassemblies.tex
        │       │   ├── transposed_report_misassemblies.tsv
        │       │   ├── transposed_report_misassemblies.txt
        │       │   ├── unaligned_report.tex
        │       │   ├── unaligned_report.tsv
        │       │   └── unaligned_report.txt
        │       ├── genome_stats
        │       │   ├── 214704-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214704-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214721-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214721-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214724-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214724-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214760-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214760-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214763-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214763-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214765-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214765-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214766-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214766-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214780-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214780-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── 214807-AF0-75-consensus-masked_gaps.txt
        │       │   ├── 214807-AF0-75-consensus-masked_genomic_features_any.txt
        │       │   ├── complete_features_histogram.pdf
        │       │   ├── features_cumulative_plot.pdf
        │       │   ├── features_frcurve_plot.pdf
        │       │   ├── genome_fraction_histogram.pdf
        │       │   └── genome_info.txt
        │       ├── icarus.html
        │       ├── icarus_viewers
        │       │   ├── alignment_viewer.html
        │       │   └── contig_size_viewer.html
        │       ├── quast.log
        │       ├── report.html
        │       ├── report.pdf
        │       ├── report.tex
        │       ├── report.tsv
        │       ├── report.txt
        │       ├── transposed_report.tex
        │       ├── transposed_report.tsv
        │       └── transposed_report.txt
        └── snpeff
            ├── 214704.AF0.75.snpEff.csv
            ├── 214704.AF0.75.snpEff.genes.txt
            ├── 214704.AF0.75.snpEff.summary.html
            ├── 214704.AF0.75.snpEff.vcf.gz
            ├── 214704.AF0.75.snpEff.vcf.gz.tbi
            ├── 214704.AF0.75.snpSift.table.txt
            ├── 214704.snpEff.csv
            ├── 214704.snpEff.genes.txt
            ├── 214704.snpEff.summary.html
            ├── 214704.snpEff.vcf.gz
            ├── 214704.snpEff.vcf.gz.tbi
            ├── 214704.snpSift.table.txt
            ├── 214721.AF0.75.snpEff.csv
            ├── 214721.AF0.75.snpEff.genes.txt
            ├── 214721.AF0.75.snpEff.summary.html
            ├── 214721.AF0.75.snpEff.vcf.gz
            ├── 214721.AF0.75.snpEff.vcf.gz.tbi
            ├── 214721.AF0.75.snpSift.table.txt
            ├── 214721.snpEff.csv
            ├── 214721.snpEff.genes.txt
            ├── 214721.snpEff.summary.html
            ├── 214721.snpEff.vcf.gz
            ├── 214721.snpEff.vcf.gz.tbi
            ├── 214721.snpSift.table.txt
            ├── 214724.AF0.75.snpEff.csv
            ├── 214724.AF0.75.snpEff.genes.txt
            ├── 214724.AF0.75.snpEff.summary.html
            ├── 214724.AF0.75.snpEff.vcf.gz
            ├── 214724.AF0.75.snpEff.vcf.gz.tbi
            ├── 214724.AF0.75.snpSift.table.txt
            ├── 214724.snpEff.csv
            ├── 214724.snpEff.genes.txt
            ├── 214724.snpEff.summary.html
            ├── 214724.snpEff.vcf.gz
            ├── 214724.snpEff.vcf.gz.tbi
            ├── 214724.snpSift.table.txt
            ├── 214760.AF0.75.snpEff.csv
            ├── 214760.AF0.75.snpEff.genes.txt
            ├── 214760.AF0.75.snpEff.summary.html
            ├── 214760.AF0.75.snpEff.vcf.gz
            ├── 214760.AF0.75.snpEff.vcf.gz.tbi
            ├── 214760.AF0.75.snpSift.table.txt
            ├── 214760.snpEff.csv
            ├── 214760.snpEff.genes.txt
            ├── 214760.snpEff.summary.html
            ├── 214760.snpEff.vcf.gz
            ├── 214760.snpEff.vcf.gz.tbi
            ├── 214760.snpSift.table.txt
            ├── 214763.AF0.75.snpEff.csv
            ├── 214763.AF0.75.snpEff.genes.txt
            ├── 214763.AF0.75.snpEff.summary.html
            ├── 214763.AF0.75.snpEff.vcf.gz
            ├── 214763.AF0.75.snpEff.vcf.gz.tbi
            ├── 214763.AF0.75.snpSift.table.txt
            ├── 214763.snpEff.csv
            ├── 214763.snpEff.genes.txt
            ├── 214763.snpEff.summary.html
            ├── 214763.snpEff.vcf.gz
            ├── 214763.snpEff.vcf.gz.tbi
            ├── 214763.snpSift.table.txt
            ├── 214765.AF0.75.snpEff.csv
            ├── 214765.AF0.75.snpEff.genes.txt
            ├── 214765.AF0.75.snpEff.summary.html
            ├── 214765.AF0.75.snpEff.vcf.gz
            ├── 214765.AF0.75.snpEff.vcf.gz.tbi
            ├── 214765.AF0.75.snpSift.table.txt
            ├── 214765.snpEff.csv
            ├── 214765.snpEff.genes.txt
            ├── 214765.snpEff.summary.html
            ├── 214765.snpEff.vcf.gz
            ├── 214765.snpEff.vcf.gz.tbi
            ├── 214765.snpSift.table.txt
            ├── 214766.AF0.75.snpEff.csv
            ├── 214766.AF0.75.snpEff.genes.txt
            ├── 214766.AF0.75.snpEff.summary.html
            ├── 214766.AF0.75.snpEff.vcf.gz
            ├── 214766.AF0.75.snpEff.vcf.gz.tbi
            ├── 214766.AF0.75.snpSift.table.txt
            ├── 214766.snpEff.csv
            ├── 214766.snpEff.genes.txt
            ├── 214766.snpEff.summary.html
            ├── 214766.snpEff.vcf.gz
            ├── 214766.snpEff.vcf.gz.tbi
            ├── 214766.snpSift.table.txt
            ├── 214780.AF0.75.snpEff.csv
            ├── 214780.AF0.75.snpEff.genes.txt
            ├── 214780.AF0.75.snpEff.summary.html
            ├── 214780.AF0.75.snpEff.vcf.gz
            ├── 214780.AF0.75.snpEff.vcf.gz.tbi
            ├── 214780.AF0.75.snpSift.table.txt
            ├── 214780.snpEff.csv
            ├── 214780.snpEff.genes.txt
            ├── 214780.snpEff.summary.html
            ├── 214780.snpEff.vcf.gz
            ├── 214780.snpEff.vcf.gz.tbi
            ├── 214780.snpSift.table.txt
            ├── 214807.AF0.75.snpEff.csv
            ├── 214807.AF0.75.snpEff.genes.txt
            ├── 214807.AF0.75.snpEff.summary.html
            ├── 214807.AF0.75.snpEff.vcf.gz
            ├── 214807.AF0.75.snpEff.vcf.gz.tbi
            ├── 214807.AF0.75.snpSift.table.txt
            ├── 214807.snpEff.csv
            ├── 214807.snpEff.genes.txt
            ├── 214807.snpEff.summary.html
            ├── 214807.snpEff.vcf.gz
            ├── 214807.snpEff.vcf.gz.tbi
            └── 214807.snpSift.table.txt
----------------------

├── ANALYSIS

├── DOC
├── RAW
│   └── MiSeq_GEN_216_10_samples
├── REFERENCES
├── RESULTS
└── TMP

date_ANALYSIS02_MET.    

    Lablog. Includes the commands to construct symbolic links with 00-reads and samples_id.txt files from the main "Analysis" directory. It also has the nextflow run command neccesary to perform the metagenomic analysis using the main.nf file using the mag repository data taken from "/processing_Data/bioinformatics/pipelines/mag/main.nf", to make date_mag folder and performed the metagnomic classification using Kraken from "/processing_Data/bioinformatics/references/kraken/minikraken_8GB_20200312.tgz" included in the _01_mag.sh. In this _01_mag.sh, the output will be redirected to a .log in order to follow the run. 


## Pipeline summary

Regarding the Analysis of variants (folder Analysis_01), samples are processed and analysed using viralrecon pipeline, which is run as follows: 

1. Download samples via SRA, ENA or GEO ids ([`ENA FTP`](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html), [`parallel-fastq-dump`](https://github.com/rvalieris/parallel-fastq-dump); *if required*)
2. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html); *if required*)
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Adapter trimming ([`fastp`](https://github.com/OpenGene/fastp))
5. Variant calling
    1. Read alignment ([`Bowtie 2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
    2. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    3. Primer sequence removal ([`iVar`](https://github.com/andersen-lab/ivar); *amplicon data only*)
    4. Duplicate read marking ([`picard`](https://broadinstitute.github.io/picard/); *removal optional*)
    5. Alignment-level QC ([`picard`](https://broadinstitute.github.io/picard/), [`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    6. Genome-wide and amplicon coverage QC plots ([`mosdepth`](https://github.com/brentp/mosdepth/))
    7. Choice of multiple variant calling and consensus sequence generation routes ([`VarScan 2`](http://dkoboldt.github.io/varscan/), [`BCFTools`](http://samtools.github.io/bcftools/bcftools.html), [`BEDTools`](https://github.com/arq5x/bedtools2/) *||* [`iVar variants and consensus`](https://github.com/andersen-lab/ivar) *||* [`BCFTools`](http://samtools.github.io/bcftools/bcftools.html), [`BEDTools`](https://github.com/arq5x/bedtools2/))
        * Variant annotation ([`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html), [`SnpSift`](http://snpeff.sourceforge.net/SnpSift.html))
        * Consensus assessment report ([`QUAST`](http://quast.sourceforge.net/quast))
    8. Intersect variants across callers ([`BCFTools`](http://samtools.github.io/bcftools/bcftools.html))
6. _De novo_ assembly
    1. Primer trimming ([`Cutadapt`](https://cutadapt.readthedocs.io/en/stable/guide.html); *amplicon data only*)
    2. Removal of host reads ([`Kraken 2`](http://ccb.jhu.edu/software/kraken2/))
    3. Choice of multiple assembly tools ([`SPAdes`](http://cab.spbu.ru/software/spades/) *||* [`metaSPAdes`](http://cab.spbu.ru/software/meta-spades/) *||* [`Unicycler`](https://github.com/rrwick/Unicycler) *||* [`minia`](https://github.com/GATB/minia))
        * Blast to reference genome ([`blastn`](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch))
        * Contiguate assembly ([`ABACAS`](https://www.sanger.ac.uk/science/tools/pagit))
        * Assembly report ([`PlasmidID`](https://github.com/BU-ISCIII/plasmidID))
        * Assembly assessment report ([`QUAST`](http://quast.sourceforge.net/quast))
        * Call variants relative to reference ([`Minimap2`](https://github.com/lh3/minimap2), [`seqwish`](https://github.com/ekg/seqwish), [`vg`](https://github.com/vgteam/vg), [`Bandage`](https://github.com/rrwick/Bandage))
        * Variant annotation ([`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html), [`SnpSift`](http://snpeff.sourceforge.net/SnpSift.html))
7. Present QC and visualisation for raw read, alignment, assembly and variant calling results ([`MultiQC`](http://multiqc.info/))

Summary of resulting data: 
.
├── bam
│   ├── SAMPLE_NAME.bam
│   ├── SAMPLE_NAME.sorted.bam
│   ├── SAMPLE_NAME.sorted.bam.bai
│   ├── SAMPLE_NAME.trim.sorted.bam
│   ├── SAMPLE_NAME.trim.sorted.bam.bai
│   ├── log
│   │   ├── SAMPLE_NAME.bowtie2.log
│   │   ├── SAMPLE_NAME.trim.ivar.log
│   ├── mosdepth
│   │   ├── amplicon
│   │   │   ├── SAMPLE_NAME.trim.amplicon.mosdepth.global.dist.txt
│   │   │   ├── SAMPLE_NAME.trim.amplicon.mosdepth.region.dist.txt
│   │   │   ├── SAMPLE_NAME.trim.amplicon.mosdepth.summary.txt
│   │   │   ├── SAMPLE_NAME.trim.amplicon.per-base.bed.gz
│   │   │   ├── SAMPLE_NAME.trim.amplicon.per-base.bed.gz.csi
│   │   │   ├── SAMPLE_NAME.trim.amplicon.regions.bed.gz
│   │   │   ├── SAMPLE_NAME.trim.amplicon.regions.bed.gz.csi
│   │   │   ├── SAMPLE_NAME.trim.amplicon.thresholds.bed.gz
│   │   │   ├── SAMPLE_NAME.trim.amplicon.thresholds.bed.gz.csi
│   │   │   └── plots
│   │   │       ├── SAMPLE_NAME.trim.amplicon.regions.coverage.pdf
│   │   │       ├── SAMPLE_NAME.trim.amplicon.regions.coverage.tsv
│   │   │       ├── all_samples.trim.amplicon.regions.coverage.tsv
│   │   │       └── all_samples.trim.amplicon.regions.heatmap.pdf
│   │   └── genome
│   │       ├── SAMPLE_NAME.trim.genome.mosdepth.global.dist.txt
│   │       ├── SAMPLE_NAME.trim.genome.mosdepth.region.dist.txt
│   │       ├── SAMPLE_NAME.trim.genome.mosdepth.summary.txt
│   │       ├── SAMPLE_NAME.trim.genome.per-base.bed.gz
│   │       ├── SAMPLE_NAME.trim.genome.per-base.bed.gz.csi
│   │       ├── SAMPLE_NAME.trim.genome.regions.bed.gz
│   │       ├── SAMPLE_NAME.trim.genome.regions.bed.gz.csi
│   │       └── plots
│   │           ├── SAMPLE_NAME.trim.genome.regions.coverage.pdf
│   │           ├── SAMPLE_NAME.trim.genome.regions.coverage.tsv
│   │           └── all_samples.trim.genome.regions.coverage.tsv
│   ├── mpileup
│   │   ├── SAMPLE_NAME.trim.mpileup
│   ├── picard_metrics
│   │   ├── SAMPLE_NAME.trim.CollectMultipleMetrics.alignment_summary_metrics
│   │   ├── SAMPLE_NAME.trim.CollectMultipleMetrics.base_distribution_by_cycle_metrics
│   │   ├── SAMPLE_NAME.trim.CollectMultipleMetrics.base_distribution_by_cycle.pdf
│   │   ├── SAMPLE_NAME.trim.CollectMultipleMetrics.insert_size_histogram.pdf
│   │   ├── SAMPLE_NAME.trim.CollectMultipleMetrics.insert_size_metrics
│   │   ├── SAMPLE_NAME.trim.CollectMultipleMetrics.quality_by_cycle_metrics
│   │   ├── SAMPLE_NAME.trim.CollectMultipleMetrics.quality_by_cycle.pdf
│   │   ├── SAMPLE_NAME.trim.CollectMultipleMetrics.quality_distribution_metrics
│   │   ├── SAMPLE_NAME.trim.CollectMultipleMetrics.quality_distribution.pdf
│   │   ├── SAMPLE_NAME.trim.CollectWgsMetrics.coverage_metrics
│   └── samtools_stats
│       ├── SAMPLE_NAME.sorted.bam.flagstat
│       ├── SAMPLE_NAME.sorted.bam.idxstats
│       ├── SAMPLE_NAME.sorted.bam.stats
│       ├── SAMPLE_NAME.trim.sorted.bam.flagstat
│       ├── SAMPLE_NAME.trim.sorted.bam.idxstats
│       ├── SAMPLE_NAME.trim.sorted.bam.stats
├── bcftools
│   ├── SAMPLE_NAME.vcf.gz
│   ├── SAMPLE_NAME.vcf.gz.tbi
│   ├── bcftools_stats
│   │   ├── SAMPLE_NAME.bcftools_stats.txt
│   ├── consensus
│   │   ├── SAMPLE_NAME.consensus.masked.fa
│   │   └── base_qc
│   │       ├── SAMPLE_NAME.ACTG_density.pdf
│   │       ├── SAMPLE_NAME.base_counts.pdf
│   │       ├── SAMPLE_NAME.base_counts.tsv
│   │       ├── SAMPLE_NAME.N_density.pdf
│   │       ├── SAMPLE_NAME.N_run.tsv
│   ├── quast
│   │   ├── aligned_stats
│   │   │   ├── cumulative_plot.pdf
│   │   │   ├── NAx_plot.pdf
│   │   │   └── NGAx_plot.pdf
│   │   ├── basic_stats
│   │   │   ├── SAMPLE_NAME.consensus.masked_GC_content_plot.pdf
│   │   │   ├── cumulative_plot.pdf
│   │   │   ├── GC_content_plot.pdf
│   │   │   ├── gc.icarus.txt
│   │   │   ├── NGx_plot.pdf
│   │   │   └── Nx_plot.pdf
│   │   ├── contigs_reports
│   │   │   ├── SAMPLE_NAME_consensus_masked.mis_contigs.fa
│   │   │   ├── all_alignments_SAMPLE_NAME-consensus-masked.tsv
│   │   │   ├── contigs_report_SAMPLE_NAME-consensus-masked.mis_contigs.info
│   │   │   ├── contigs_report_SAMPLE_NAME-consensus-masked.stderr
│   │   │   ├── contigs_report_SAMPLE_NAME-consensus-masked.stdout
│   │   │   ├── contigs_report_SAMPLE_NAME-consensus-masked.unaligned.info
│   │   │   ├── minimap_output
│   │   │   │   ├── SAMPLE_NAME-consensus-masked.coords
│   │   │   │   ├── SAMPLE_NAME-consensus-masked.coords.filtered
│   │   │   │   ├── SAMPLE_NAME-consensus-masked.coords_tmp
│   │   │   │   ├── SAMPLE_NAME-consensus-masked.sf
│   │   │   │   ├── SAMPLE_NAME-consensus-masked.unaligned
│   │   │   │   ├── SAMPLE_NAME-consensus-masked.used_snps.gz
│   │   │   ├── misassemblies_frcurve_plot.pdf
│   │   │   ├── misassemblies_plot.pdf
│   │   │   ├── misassemblies_report.tex
│   │   │   ├── misassemblies_report.tsv
│   │   │   ├── misassemblies_report.txt
│   │   │   ├── transposed_report_misassemblies.tex
│   │   │   ├── transposed_report_misassemblies.tsv
│   │   │   ├── transposed_report_misassemblies.txt
│   │   │   ├── unaligned_report.tex
│   │   │   ├── unaligned_report.tsv
│   │   │   └── unaligned_report.txt
│   │   ├── genome_stats
│   │   │   ├── SAMPLE_NAME-consensus-masked_gaps.txt
│   │   │   ├── SAMPLE_NAME-consensus-masked_genomic_features_any.txt
│   │   │   ├── complete_features_histogram.pdf
│   │   │   ├── features_cumulative_plot.pdf
│   │   │   ├── features_frcurve_plot.pdf
│   │   │   ├── genome_fraction_histogram.pdf
│   │   │   └── genome_info.txt
│   │   ├── icarus.html
│   │   ├── icarus_viewers
│   │   │   ├── alignment_viewer.html
│   │   │   └── contig_size_viewer.html
│   │   ├── quast.log
│   │   ├── report.html
│   │   ├── report.pdf
│   │   ├── report.tex
│   │   ├── report.tsv
│   │   ├── report.txt
│   │   ├── transposed_report.tex
│   │   ├── transposed_report.tsv
│   │   └── transposed_report.txt
│   └── snpeff
│       ├── SAMPLE_NAME.snpEff.csv
│       ├── SAMPLE_NAME.snpEff.genes.txt
│       ├── SAMPLE_NAME.snpEff.summary.html
│       ├── SAMPLE_NAME.snpEff.vcf.gz
│       ├── SAMPLE_NAME.snpEff.vcf.gz.tbi
│       ├── SAMPLE_NAME.snpSift.table.txt
├── intersect
│   ├── SAMPLE_NAME
│   │   ├── 0000.vcf.gz
│   │   ├── 0000.vcf.gz.tbi
│   │   ├── 0001.vcf.gz
│   │   ├── 0001.vcf.gz.tbi
│   │   ├── 0002.vcf.gz
│   │   ├── 0002.vcf.gz.tbi
│   │   ├── README.txt
│   │   └── sites.txt
├── ivar
│   ├── SAMPLE_NAME.AF0.75.vcf.gz
│   ├── SAMPLE_NAME.AF0.75.vcf.gz.tbi
│   ├── SAMPLE_NAME.tsv
│   ├── SAMPLE_NAME.vcf.gz
│   ├── SAMPLE_NAME.vcf.gz.tbi
│   ├── bcftools_stats
│   │   ├── SAMPLE_NAME.AF0.75.bcftools_stats.txt
│   │   ├── SAMPLE_NAME.bcftools_stats.txt
│   ├── consensus
│   │   ├── SAMPLE_NAME.AF0.75.consensus.fa
│   │   ├── SAMPLE_NAME.AF0.75.consensus.qual.txt
│   │   └── base_qc
│   │       ├── SAMPLE_NAME.AF0.75.ACTG_density.pdf
│   │       ├── SAMPLE_NAME.AF0.75.base_counts.pdf
│   │       ├── SAMPLE_NAME.AF0.75.base_counts.tsv
│   │       ├── SAMPLE_NAME.AF0.75.N_density.pdf
│   │       ├── SAMPLE_NAME.AF0.75.N_run.tsv
│   │       ├── SAMPLE_NAME.AF0.75.R_density.pdf
│   │       ├── SAMPLE_NAME.AF0.75.Y_density.pdf
│   ├── log
│   │   ├── SAMPLE_NAME.AF0.75.variant.counts.log
│   │   ├── SAMPLE_NAME.variant.counts.log
│   ├── quast
│   │   └── AF0.75
│   │       ├── aligned_stats
│   │       │   ├── cumulative_plot.pdf
│   │       │   ├── NAx_plot.pdf
│   │       │   └── NGAx_plot.pdf
│   │       ├── basic_stats
│   │       │   ├── SAMPLE_NAME.AF0.75.consensus_GC_content_plot.pdf
│   │       │   ├── cumulative_plot.pdf
│   │       │   ├── GC_content_plot.pdf
│   │       │   ├── gc.icarus.txt
│   │       │   ├── NGx_plot.pdf
│   │       │   └── Nx_plot.pdf
│   │       ├── contigs_reports
│   │       │   ├── SAMPLE_NAME_AF0_75_consensus.mis_contigs.fa
│   │       │   ├── all_alignments_SAMPLE_NAME-AF0-75-consensus.tsv
│   │       │   ├── contigs_report_SAMPLE_NAME-AF0-75-consensus.mis_contigs.info
│   │       │   ├── contigs_report_SAMPLE_NAME-AF0-75-consensus.stderr
│   │       │   ├── contigs_report_SAMPLE_NAME-AF0-75-consensus.stdout
│   │       │   ├── contigs_report_SAMPLE_NAME-AF0-75-consensus.unaligned.info
│   │       │   ├── minimap_output
│   │       │   │   ├── SAMPLE_NAME-AF0-75-consensus.coords
│   │       │   │   ├── SAMPLE_NAME-AF0-75-consensus.coords.filtered
│   │       │   │   ├── SAMPLE_NAME-AF0-75-consensus.coords_tmp
│   │       │   │   ├── SAMPLE_NAME-AF0-75-consensus.sf
│   │       │   │   ├── SAMPLE_NAME-AF0-75-consensus.unaligned
│   │       │   │   ├── SAMPLE_NAME-AF0-75-consensus.used_snps.gz
│   │       │   ├── misassemblies_frcurve_plot.pdf
│   │       │   ├── misassemblies_plot.pdf
│   │       │   ├── misassemblies_report.tex
│   │       │   ├── misassemblies_report.tsv
│   │       │   ├── misassemblies_report.txt
│   │       │   ├── transposed_report_misassemblies.tex
│   │       │   ├── transposed_report_misassemblies.tsv
│   │       │   ├── transposed_report_misassemblies.txt
│   │       │   ├── unaligned_report.tex
│   │       │   ├── unaligned_report.tsv
│   │       │   └── unaligned_report.txt
│   │       ├── genome_stats
│   │       │   ├── SAMPLE_NAME-AF0-75-consensus_gaps.txt
│   │       │   ├── SAMPLE_NAME-AF0-75-consensus_genomic_features_any.txt
│   │       │   ├── complete_features_histogram.pdf
│   │       │   ├── features_cumulative_plot.pdf
│   │       │   ├── features_frcurve_plot.pdf
│   │       │   ├── genome_fraction_histogram.pdf
│   │       │   └── genome_info.txt
│   │       ├── icarus.html
│   │       ├── icarus_viewers
│   │       │   ├── alignment_viewer.html
│   │       │   └── contig_size_viewer.html
│   │       ├── quast.log
│   │       ├── report.html
│   │       ├── report.pdf
│   │       ├── report.tex
│   │       ├── report.tsv
│   │       ├── report.txt
│   │       ├── transposed_report.tex
│   │       ├── transposed_report.tsv
│   │       └── transposed_report.txt
│   └── snpeff
│       ├── SAMPLE_NAME.AF0.75.snpEff.csv
│       ├── SAMPLE_NAME.AF0.75.snpEff.genes.txt
│       ├── SAMPLE_NAME.AF0.75.snpEff.summary.html
│       ├── SAMPLE_NAME.AF0.75.snpEff.vcf.gz
│       ├── SAMPLE_NAME.AF0.75.snpEff.vcf.gz.tbi
│       ├── SAMPLE_NAME.AF0.75.snpSift.table.txt
│       ├── SAMPLE_NAME.snpEff.csv
│       ├── SAMPLE_NAME.snpEff.genes.txt
│       ├── SAMPLE_NAME.snpEff.summary.html
│       ├── SAMPLE_NAME.snpEff.vcf.gz
│       ├── SAMPLE_NAME.snpEff.vcf.gz.tbi
│       ├── SAMPLE_NAME.snpSift.table.txt
│       ├── snpSift_template2.txt
│       ├── snpSift_template_filtered.txt
│       └── snpSift_template.txt
├── summary_variants_metrics_mqc.tsv
└── varscan2
    ├── SAMPLE_NAME.AF0.75.vcf.gz
    ├── SAMPLE_NAME.AF0.75.vcf.gz.tbi
    ├── SAMPLE_NAME.vcf.gz
    ├── SAMPLE_NAME.vcf.gz.tbi
    ├── bcftools_stats
    │   ├── SAMPLE_NAME.AF0.75.bcftools_stats.txt
    │   ├── SAMPLE_NAME.bcftools_stats.txt
    ├── consensus
    │   ├── SAMPLE_NAME.AF0.75.consensus.masked.fa
    │   ├── base_qc
    │   │   ├── SAMPLE_NAME.AF0.75.ACTG_density.pdf
    │   │   ├── SAMPLE_NAME.AF0.75.base_counts.pdf
    │   │   ├── SAMPLE_NAME.AF0.75.base_counts.tsv
    │   │   ├── SAMPLE_NAME.AF0.75.N_density.pdf
    │   │   ├── SAMPLE_NAME.AF0.75.N_run.tsv
    │   └── snpSift_template.txt
    ├── log
    │   ├── SAMPLE_NAME.varscan2.log
    ├── quast
    │   └── AF0.75
    │       ├── aligned_stats
    │       │   ├── cumulative_plot.pdf
    │       │   ├── NAx_plot.pdf
    │       │   └── NGAx_plot.pdf
    │       ├── basic_stats
    │       │   ├── SAMPLE_NAME.AF0.75.consensus.masked_GC_content_plot.pdf
    │       │   ├── cumulative_plot.pdf
    │       │   ├── GC_content_plot.pdf
    │       │   ├── gc.icarus.txt
    │       │   ├── NGx_plot.pdf
    │       │   └── Nx_plot.pdf
    │       ├── contigs_reports
    │       │   ├── SAMPLE_NAME_AF0_75_consensus_masked.mis_contigs.fa
    │       │   ├── all_alignments_SAMPLE_NAME-AF0-75-consensus-masked.tsv
    │       │   ├── contigs_report_SAMPLE_NAME-AF0-75-consensus-masked.mis_contigs.info
    │       │   ├── contigs_report_SAMPLE_NAME-AF0-75-consensus-masked.stderr
    │       │   ├── contigs_report_SAMPLE_NAME-AF0-75-consensus-masked.stdout
    │       │   ├── contigs_report_SAMPLE_NAME-AF0-75-consensus-masked.unaligned.info
    │       │   ├── minimap_output
    │       │   │   ├── SAMPLE_NAME-AF0-75-consensus-masked.coords
    │       │   │   ├── SAMPLE_NAME-AF0-75-consensus-masked.coords.filtered
    │       │   │   ├── SAMPLE_NAME-AF0-75-consensus-masked.coords_tmp
    │       │   │   ├── SAMPLE_NAME-AF0-75-consensus-masked.sf
    │       │   │   ├── SAMPLE_NAME-AF0-75-consensus-masked.unaligned
    │       │   │   ├── SAMPLE_NAME-AF0-75-consensus-masked.used_snps.gz
    │       │   ├── misassemblies_frcurve_plot.pdf
    │       │   ├── misassemblies_plot.pdf
    │       │   ├── misassemblies_report.tex
    │       │   ├── misassemblies_report.tsv
    │       │   ├── misassemblies_report.txt
    │       │   ├── transposed_report_misassemblies.tex
    │       │   ├── transposed_report_misassemblies.tsv
    │       │   ├── transposed_report_misassemblies.txt
    │       │   ├── unaligned_report.tex
    │       │   ├── unaligned_report.tsv
    │       │   └── unaligned_report.txt
    │       ├── genome_stats
    │       │   ├── SAMPLE_NAME-AF0-75-consensus-masked_gaps.txt
    │       │   ├── SAMPLE_NAME-AF0-75-consensus-masked_genomic_features_any.txt
    │       │   ├── complete_features_histogram.pdf
    │       │   ├── features_cumulative_plot.pdf
    │       │   ├── features_frcurve_plot.pdf
    │       │   ├── genome_fraction_histogram.pdf
    │       │   └── genome_info.txt
    │       ├── icarus.html
    │       ├── icarus_viewers
    │       │   ├── alignment_viewer.html
    │       │   └── contig_size_viewer.html
    │       ├── quast.log
    │       ├── report.html
    │       ├── report.pdf
    │       ├── report.tex
    │       ├── report.tsv
    │       ├── report.txt
    │       ├── transposed_report.tex
    │       ├── transposed_report.tsv
    │       └── transposed_report.txt
    └── snpeff
        ├── SAMPLE_NAME.AF0.75.snpEff.csv
        ├── SAMPLE_NAME.AF0.75.snpEff.genes.txt
        ├── SAMPLE_NAME.AF0.75.snpEff.summary.html
        ├── SAMPLE_NAME.AF0.75.snpEff.vcf.gz
        ├── SAMPLE_NAME.AF0.75.snpEff.vcf.gz.tbi
        ├── SAMPLE_NAME.AF0.75.snpSift.table.txt

> **NB:** The pipeline has a number of options to allow you to run only specific aspects of the workflow if you so wish.
For example, you can skip all of the assembly steps with the `--skip_assembly` parameter.
See the [usage docs](docs/usage.md) for all of the available options when running the pipeline.

## Pipeline reporting

Numerous QC and reporting steps are included in the pipeline in order to collate a full summary of the analysis within a single [MultiQC](https://multiqc.info/) report. You can see [an example MultiQC report here](https://raw.githack.com/nf-core/viralrecon/master/docs/html/multiqc_report.html), generated using the parameters defined in [this configuration file](https://github.com/nf-core/viralrecon/blob/master/conf/test_full.config). The pipeline was run with [these samples](https://zenodo.org/record/3735111), prepared from the [ncov-2019 ARTIC Network V1 amplicon set](https://artic.network/ncov-2019) and sequenced on the Illumina MiSeq platform in 301bp paired-end format.

## Results and statistics. 

From the files generated after the variant obtention, the next data are extracted and mounted in a tab separated file that will be used to analysed the statistics from every sample sequenced. 

Host (Human/bison)
VirusSequence (reference genome, the first Wuhan SARS-CoV2 sequenced: NC_045512.2
Sample
Total reads (*.sorted.bam.flagstat)
Reads_host R1 (*.kraken2.report.txt)
Reads_host_total 
%_reads_host (*.kraken.report.txt)
Reads_virus_total (.sorted.bam.flagstat)
%_reads_virus 
Unmapped_reads 
%Unmapped_reads
Mean DP Coverage (*.trim.CollectWgsMetrics.coverage_metrics)
PCT_10X (*.trim.CollectWgsMetrics.coverage_metrics)
Variants_consensusx10 (*.AF0.75.snpSift.table.txt)
Missense_variants (*.AF0.75.snpSift.table.txt)
%Ns_10x (_Python_script)
Lineage (Pangolin software)

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/viralrecon -profile test,<docker/singularity/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

See the [usage documentation](docs/usage.md) for all of the available options when running the pipeline.

## Documentation

The nf-core/viralrecon pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](docs/usage.md#reference-genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits
