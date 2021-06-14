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
