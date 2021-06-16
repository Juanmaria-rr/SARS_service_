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

In order to make reproducible the analysis run by BU_ISCIII, here we describe step by step how SARS_service is performed in BU_ISCIII using the combination of viralrecon pipeline with Pangolin and Kraken, taking adventage of the Nextflow workflow framework. 

## Viene de DSL2??# 
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with Docker containers making installation trivial and results highly reproducible. Furthermore, automated continuous integration tests that run the pipeline on a full-sized dataset using AWS cloud ensure that the code is stable.

## Local data organization summary

Once every service is ordered, a template of data structure is created, starting with a new folder with the name of the service, the date, the type of project (i.e. SARS-Cov) and the name of the researcher, typically "SRVxxxx_date_SARS_name-researcher". 
In this folder there will be other directories: 

```ruby
├── ANALYSIS
│   ├── 00-reads
│   ├── $(date '+%Y%m%d')_ANALYSIS01_AMPLICONS_HUMAN
│   └── $(date '+%Y%m%d')_ANALYSIS02_MET
├── DOC
├── RAW
│   └── <samples_folder>
├── REFERENCES
├── RESULTS
└── TMP

```

In ANALYSIS folder is where all the analysis process will run. In RAW, inside a folder named with the RUN  there are the R1/2 fastq.gz files, where they will be taken from in order to perform the analysis. DOC, REFERENCES, RESULTS and TMP will be empty. 

**ANALYSIS_folder**

```ruby
ANALYSIS
├── 00-reads
├── 20210531_ANALYSIS01_AMPLICONS_HUMAN
├── 20210531_ANALYSIS02_MET
├── lablog
└── samples_id.txt
```
Lablog - 

It contains the commands necesary to make the sample_id.txt file, which is derived from all the samples located in the RAW folder. It also has the commands to create the other two folders in which the analysis will be splited:

```ruby
ls ../RAW/*/* | tr '\/' '\t' | cut -f4 | cut -d "_" -f 1 | sort -u | grep -v "md5" > samples_id.txt
mkdir -p 00-reads
mkdir -p $(date '+%Y%m%d')_ANALYSIS01_AMPLICONS_HUMAN
mkdir -p $(date '+%Y%m%d')_ANALYSIS02_MET
```
samples_id.txt - 

Is a list with the <sample_name> derived from the samples files .fastq.gz stored in the folder RAW.

**/00-reads. 

```ruby
ANALYSIS/00-reads/
└── lablog
```
Its Lablog file contains the commands that will make symbolic links to every RAW/<RUN_name>/R1/2.fastq.gz files to be analysed.

```ruby
cat ../samples_id.txt | xargs -I % echo "ln -s ../../RAW/*/%_*R1*.fastq.gz %_R1.fastq.gz" | bash
cat ../samples_id.txt | xargs -I % echo "ln -s ../../RAW/*/%_*R2*.fastq.gz %_R2.fastq.gz" | bash
```
Once these symbolic links are constructed, the folder schema will look like this: 

```ruby
ANALYSIS/00-reads/
├── <sample_name>_R1.fastq.gz -> ../../RAW/<RUN_name>_samples/<sample_name>_R1_001.fastq.gz
├── <sample_name>_R2.fastq.gz -> ../../RAW/<RUN_name>_samples/<sample_name>_001.fastq.gz
└── lablog
```

/date_ANALYSIS01_AMPLICONS_HUMAN

```ruby
20210531_ANALYSIS01_AMPLICONS_HUMAN/
└── lablog
```
**Lablog

To note: this lablog has to be copied from the previous service performed. This lablog is the essential file which all the analysis and directories organization will be contructed from. IThis lablog perform the next tasks: 


date_ANALYSIS01_AMPLICONS_HUMAN.  


1) It makes symbolic links with 00-reads folder and samples_id.txt file, which contain the raw data (Fastqc format) and sample names. Moreover, it creates the samplesheet.csv needed for the analyses. 

```ruby
ln -s ../00-reads .
ln -s ../samples_id.txt .
echo "sample,fastq_1,fastq_2" > samplesheet.csv
cat samples_id.txt | while read in; do echo "${in},00-reads/${in}_R1.fastq.gz,00-reads/${in}_R2.fastq.gz"; done >> samplesheet.csv
```
After running those commands, the directory will look as shown below:

```ruby
date_ANALYSIS01_AMPLICONS_HUMAN/
├── 00-reads -> ../00-reads
├── lablog
├── samplesheet.csv
└── samples_id.txt -> ../samples_id.txt

```

2) In addition, the lablog contains the commands necessary to run viralrecon_prod using main.nf from the repository located locally in processing_Data/bioinformatics/pipelines/viralrecon_prod/main.nf. The commands will be written in a bash script file in order to run the analysis by bash. The output from the run is allways stored in a .log file, that allows us to follow the process apart form the kernell. 

```ruby
echo "nextflow run /processing_Data/bioinformatics/pipelines/viralrecon_prod/main.nf -bg --input samplesheet.csv -profile conda --outdir $(date '+%Y%m%d')_viralrecon_mapping --assemblers none --amplicon_fasta 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/genome/NC_045512.2/amplicon/nCoV-2019.artic.V3.primer.fasta' --amplicon_bed 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/genome/NC_045512.2/amplicon/nCoV-2019.artic.V3.bed' --kraken2_db /processing_Data/bioinformatics/references/eukaria/homo_sapiens/hg38/UCSC/kraken2/kraken2_human.tar.gz --protocol amplicon --genome NC_045512.2 --save_align_intermeds --skip_vg --skip_markduplicates --save_mpileup --min_allele_freq 0 -resume" > _01_viralrecon_mapping.sh
```

3) Moreover, a bash script to create a summary_report and a python script to analyse the percentages of N in every analysed samples are copied in this folder using two commands in the lablog:

```ruby
cp /processing_Data/bioinformatics/services_and_colaborations/CNM/virologia/SRVCNM327_20210201_SARSCOV228_icasas_S/ANALYSIS/20210201_ANALYSIS01_AMPLICONS_HUMAN/create_summary_report.sh .
cp /processing_Data/bioinformatics/services_and_colaborations/CNM/virologia/SRVCNM327_20210201_SARSCOV228_icasas_S/ANALYSIS/20210201_ANALYSIS01_AMPLICONS_HUMAN/percentajeNs.py .
```
`

4) Finally, the directories of $(date '+%Y%m%d')_Pangolin analyses and $(date '+%Y%m%d')variants_table are created. This folderds will include the following results from Pangoling lineages classification and the results of SARS-CoV2 variants. 

```ruby
mkdir -p $(date '+%Y%m%d')_pangolin
mkdir -p $(date '+%Y%m%d')_variants_table
```
At the end of using lablog file commands, the full directory will have the next structure:

```ruby
date_ANALYSIS01_AMPLICONS_HUMAN/
├── 00-reads -> ../00-reads
├── _01_viralrecon_mapping.sh
├── create_summary_report.sh**
├── lablog
├── pangolin
├── percentajeNs.py
├── samplesheet.csv
├── samples_id.txt -> ../samples_id.txt
├── variants_table
└── work
``` 
***Launching Viralrecon analysis

Once this, we are ready to run viralrecon analysis, which will be launched using nohup command in order to allow the computer to run other tasks independetly of our ssh conection: 

```ruby
nohup bash _01_viralrecon_mapping.sh &> $(date '+%Y%m%d')_mapping01.log &
```
As you may note, the output from the viralrecon pipeline is directed to a .log file, that will allow us to follow in real time the process while our console will be free to run other tasks. Do ```ruby tail -f $(date '+%Y%m%d')_mapping01.log ``` to supervise the pipeline. 

After this, a new folder will be created, named /$(date '+%Y%m%d')_viralrecon_mapping, where every output from the pipeline will be located, with the following structure. In addition, nextflow will create the folder /work, in which every task is isolated in folders. 


```ruby
date_ANALYSIS01_AMPLICONS_HUMAN/
├── 00-reads -> ../00-reads
├── _01_viralrecon_mapping.sh
├── 20210603_mapping01.log
├── 20210603_viralrecon_mapping
├── create_summary_report.sh**
├── lablog
├── pangolin
├── percentajeNs.py
├── samplesheet.csv
├── samples_id.txt -> ../samples_id.txt
├── variants_table
└── work
``` 

Now, all the data coming from the viralrecon analysis is located in /20210603_viralrecon_mapping folder, which will take the following structure if the analysis is performed correctly (only directories shown).

```ruby
date_ANALYSIS01_AMPLICONS_HUMAN/20210603_viralrecon_mapping
├── assembly
│   ├── cutadapt
│   └── kraken2
├── multiqc
│   └── multiqc_data
├── pipeline_info
├── preprocess
│   ├── fastp
│   └── fastqc
└── variants
    ├── bam
    ├── bcftools
    ├── intersect
    ├── ivar
    └── varscan2
```
    date_viralrecon_mapping

** Only directories.
```ruby
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
    │   ├── Sample_1
    │   ├── Sample_2
    │   ├── Sample_3
    │   ├── Sample_4
    │   ├── Sample_5
    │   ├── Sample_6
    │   ├── Sample_7
    │   ├── Sample_8
    │   └── Sample_9
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
        
**Directories + sample_name files
```ruby

├── assembly
│   ├── cutadapt
│   │   ├── fastqc
│   │   │   ├── 214704_1.ptrim_fastqc.html
│   │   │   ├── 214704_2.ptrim_fastqc.html
│   │   │   └── zips
│   │   │       ├── 214704_1.ptrim_fastqc.zip
│   │   │       ├── 214704_2.ptrim_fastqc.zip
│   │   └── log
│   │       ├── 214704.cutadapt.log
│   ├── kraken2
│   │   ├── 214704.kraken2.report.txt
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
│   │   ├── fastqc
│   │   │   ├── 214704_1.trim_fastqc.html
│   │   │   ├── 214704_2.trim_fastqc.html
│   │   │   └── zips
│   │   │       ├── 214704_1.trim_fastqc.zip
│   │   │       ├── 214704_2.trim_fastqc.zip
│   │   └── log
│   │       ├── 214704.fastp.log
│   └── fastqc
│       ├── 214704_1.merged_fastqc.html
│       ├── 214704_2.merged_fastqc.html
│       └── zips
│           ├── 214704_1.merged_fastqc.zip
│           ├── 214704_2.merged_fastqc.zip
└── variants
    ├── bam
    │   ├── 214704.bam
    │   ├── 214704.sorted.bam
    │   ├── 214704.sorted.bam.bai
    │   ├── 214704.trim.sorted.bam
    │   ├── 214704.trim.sorted.bam.bai
    │   ├── log
    │   │   ├── 214704.bowtie2.log
    │   │   ├── 214704.trim.ivar.log
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
    │   │   │   └── plots
    │   │   │       ├── 214704.trim.amplicon.regions.coverage.pdf
    │   │   │       ├── 214704.trim.amplicon.regions.coverage.tsv
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
    │   │       └── plots
    │   │           ├── 214704.trim.genome.regions.coverage.pdf
    │   │           ├── 214704.trim.genome.regions.coverage.tsv
    │   │           └── all_samples.trim.genome.regions.coverage.tsv
    │   ├── mpileup
    │   │   ├── 214704.trim.mpileup
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
    │   └── samtools_stats
    │       ├── 214704.sorted.bam.flagstat
    │       ├── 214704.sorted.bam.idxstats
    │       ├── 214704.sorted.bam.stats
    │       ├── 214704.trim.sorted.bam.flagstat
    │       ├── 214704.trim.sorted.bam.idxstats
    │       ├── 214704.trim.sorted.bam.stats
    ├── bcftools
    │   ├── 214704.vcf.gz
    │   ├── 214704.vcf.gz.tbi
    │   ├── bcftools_stats
    │   │   ├── 214704.bcftools_stats.txt
    │   ├── consensus
    │   │   ├── 214704.consensus.masked.fa
    │   │   └── base_qc
    │   │       ├── 214704.ACTG_density.pdf
    │   │       ├── 214704.base_counts.pdf
    │   │       ├── 214704.base_counts.tsv
    │   │       ├── 214704.N_density.pdf
    │   │       ├── 214704.N_run.tsv
    │   ├── quast
    │   │   ├── aligned_stats
    │   │   │   ├── cumulative_plot.pdf
    │   │   │   ├── NAx_plot.pdf
    │   │   │   └── NGAx_plot.pdf
    │   │   ├── basic_stats
    │   │   │   ├── 214704.consensus.masked_GC_content_plot.pdf
    │   │   │   ├── cumulative_plot.pdf
    │   │   │   ├── GC_content_plot.pdf
    │   │   │   ├── gc.icarus.txt
    │   │   │   ├── NGx_plot.pdf
    │   │   │   └── Nx_plot.pdf
    │   │   ├── contigs_reports
    │   │   │   ├── 214704_consensus_masked.mis_contigs.fa
    │   │   │   ├── all_alignments_214704-consensus-masked.tsv
    │   │   │   ├── contigs_report_214704-consensus-masked.mis_contigs.info
    │   │   │   ├── contigs_report_214704-consensus-masked.stderr
    │   │   │   ├── contigs_report_214704-consensus-masked.stdout
    │   │   │   ├── contigs_report_214704-consensus-masked.unaligned.info
    │   │   │   ├── minimap_output
    │   │   │   │   ├── 214704-consensus-masked.coords
    │   │   │   │   ├── 214704-consensus-masked.coords.filtered
    │   │   │   │   ├── 214704-consensus-masked.coords_tmp
    │   │   │   │   ├── 214704-consensus-masked.sf
    │   │   │   │   ├── 214704-consensus-masked.unaligned
    │   │   │   │   ├── 214704-consensus-masked.used_snps.gz
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
    ├── intersect
    │   ├── 214704
    │       ├── 0000.vcf.gz
    │       ├── 0000.vcf.gz.tbi
    │       ├── 0001.vcf.gz
    │       ├── 0001.vcf.gz.tbi
    │       ├── 0002.vcf.gz
    │       ├── 0002.vcf.gz.tbi
    │       ├── README.txt
    │       └── sites.txt
    ├── ivar
    │   ├── 214704.AF0.75.vcf.gz
    │   ├── 214704.AF0.75.vcf.gz.tbi
    │   ├── 214704.tsv
    │   ├── 214704.vcf.gz
    │   ├── 214704.vcf.gz.tbi
    │   ├── bcftools_stats
    │   │   ├── 214704.AF0.75.bcftools_stats.txt
    │   │   ├── 214704.bcftools_stats.txt
    │   ├── consensus
    │   │   ├── 214704.AF0.75.consensus.fa
    │   │   ├── 214704.AF0.75.consensus.qual.txt
    │   │   └── base_qc
    │   │       ├── 214704.AF0.75.ACTG_density.pdf
    │   │       ├── 214704.AF0.75.base_counts.pdf
    │   │       ├── 214704.AF0.75.base_counts.tsv
    │   │       ├── 214704.AF0.75.N_density.pdf
    │   │       ├── 214704.AF0.75.N_run.tsv
    │   │       ├── 214704.AF0.75.R_density.pdf
    │   │       ├── 214704.AF0.75.Y_density.pdf
    │   ├── log
    │   │   ├── 214704.AF0.75.variant.counts.log
    │   │   ├── 214704.variant.counts.log
    │   ├── quast
    │   │   └── AF0.75
    │   │       ├── aligned_stats
    │   │       │   ├── cumulative_plot.pdf
    │   │       │   ├── NAx_plot.pdf
    │   │       │   └── NGAx_plot.pdf
    │   │       ├── basic_stats
    │   │       │   ├── 214704.AF0.75.consensus_GC_content_plot.pdf
    │   │       │   ├── cumulative_plot.pdf
    │   │       │   ├── GC_content_plot.pdf
    │   │       │   ├── gc.icarus.txt
    │   │       │   ├── NGx_plot.pdf
    │   │       │   └── Nx_plot.pdf
    │   │       ├── contigs_reports
    │   │       │   ├── 214704_AF0_75_consensus.mis_contigs.fa
    │   │       │   ├── all_alignments_214704-AF0-75-consensus.tsv
    │   │       │   ├── contigs_report_214704-AF0-75-consensus.mis_contigs.info
    │   │       │   ├── contigs_report_214704-AF0-75-consensus.stderr
    │   │       │   ├── contigs_report_214704-AF0-75-consensus.stdout
    │   │       │   ├── contigs_report_214704-AF0-75-consensus.unaligned.info
    │   │       │   ├── minimap_output
    │   │       │   │   ├── 214704-AF0-75-consensus.coords
    │   │       │   │   ├── 214704-AF0-75-consensus.coords.filtered
    │   │       │   │   ├── 214704-AF0-75-consensus.coords_tmp
    │   │       │   │   ├── 214704-AF0-75-consensus.sf
    │   │       │   │   ├── 214704-AF0-75-consensus.unaligned
    │   │       │   │   ├── 214704-AF0-75-consensus.used_snps.gz
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
    │       ├── snpSift_template2.txt
    │       ├── snpSift_template_filtered.txt
    │       └── snpSift_template.txt
    ├── summary_variants_metrics_mqc.tsv
    └── varscan2
        ├── 214704.AF0.75.vcf.gz
        ├── 214704.AF0.75.vcf.gz.tbi
        ├── 214704.vcf.gz
        ├── 214704.vcf.gz.tbi
        ├── bcftools_stats
        │   ├── 214704.AF0.75.bcftools_stats.txt
        │   ├── 214704.bcftools_stats.txt
        ├── consensus
        │   ├── 214704.AF0.75.consensus.masked.fa
        │   ├── base_qc
        │   │   ├── 214704.AF0.75.ACTG_density.pdf
        │   │   ├── 214704.AF0.75.base_counts.pdf
        │   │   ├── 214704.AF0.75.base_counts.tsv
        │   │   ├── 214704.AF0.75.N_density.pdf
        │   │   ├── 214704.AF0.75.N_run.tsv
        │   └── snpSift_template.txt
        ├── log
        │   ├── 214704.varscan2.log
        ├── quast
        │   └── AF0.75
        │       ├── aligned_stats
        │       │   ├── cumulative_plot.pdf
        │       │   ├── NAx_plot.pdf
        │       │   └── NGAx_plot.pdf
        │       ├── basic_stats
        │       │   ├── 214704.AF0.75.consensus.masked_GC_content_plot.pdf
        │       │   ├── cumulative_plot.pdf
        │       │   ├── GC_content_plot.pdf
        │       │   ├── gc.icarus.txt
        │       │   ├── NGx_plot.pdf
        │       │   └── Nx_plot.pdf
        │       ├── contigs_reports
        │       │   ├── 214704_AF0_75_consensus_masked.mis_contigs.fa
        │       │   ├── all_alignments_214704-AF0-75-consensus-masked.tsv
        │       │   ├── contigs_report_214704-AF0-75-consensus-masked.mis_contigs.info
        │       │   ├── contigs_report_214704-AF0-75-consensus-masked.stderr
        │       │   ├── contigs_report_214704-AF0-75-consensus-masked.stdout
        │       │   ├── contigs_report_214704-AF0-75-consensus-masked.unaligned.info
        │       │   ├── minimap_output
        │       │   │   ├── 214704-AF0-75-consensus-masked.coords
        │       │   │   ├── 214704-AF0-75-consensus-masked.coords.filtered
        │       │   │   ├── 214704-AF0-75-consensus-masked.coords_tmp
        │       │   │   ├── 214704-AF0-75-consensus-masked.sf
        │       │   │   ├── 214704-AF0-75-consensus-masked.unaligned
        │       │   │   ├── 214704-AF0-75-consensus-masked.used_snps.gz
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
```

----------------------


date_ANALYSIS02_MET.    
└──  lablog


**Lablog

Includes the commands to construct symbolic links with 00-reads and samples_id.txt files from the main "Analysis" directory. It also has the nextflow run command neccesary to perform the metagenomic analysis using the main.nf file using the mag repository data taken from "/processing_Data/bioinformatics/pipelines/mag/main.nf", to make date_mag folder and performed the metagnomic classification using Kraken from "/processing_Data/bioinformatics/references/kraken/minikraken_8GB_20200312.tgz" included in the _01_mag.sh. As in viralrecon analysis, the output of _01_mag.sh will be redirected to a .log in order to follow the run. The content of the lablog is shown in the following lines: 

```ruby
ln -s ../00-reads .
ln -s ../samples_id.txt .
echo "nextflow run /processing_Data/bioinformatics/pipelines/mag/main.nf -bg --reads '00-reads/*_R{1,2}.fastq.gz' -profile hpc_isciii --outdir $(date '+%Y%m%d')_mag --kraken2_db /processing_Data/bioinformatics/references/kraken/minikraken_8GB_20200312.tgz --skip_busco --skip_spades --skip_spadeshybrid --skip_megahit -resume" > _01_mag.sh
```

After this, the directory should be like this: 

```ruby
.
├── 00-reads -> ../00-reads
├── _01_mag.sh
├── lablog
└── samples_id.txt -> ../samples_id.txt
```

***Launching metagenomic analysis analyses

As seen for viralrecon analyses, we run this pipeline using nextflow and the nohup command: 

```ruby
#nohup bash _01_mag.sh &> $(date '+%Y%m%d')_mag01.log &
```
Similarly to what viralrecon, the directory is modified with new folders that store the data and results derived from the analysis: 

```ruby
.
├── 00-reads -> ../00-reads
├── _01_mag.sh
├── _02_mag.sh
├── 20210603_mag
├── 20210607_mag01.log
├── lablog
├── samples_id.txt -> ../samples_id.txt
└── work
```
Now, /$(date '+%Y%m%d')_mag folder include the next directories regarding metagnomic analyses: 

```ruby
.
├── pipeline_info
├── QC_shortreads
│   ├── fastp
│   │   ├── <sample_name>
│   ├── fastqc
│   └── remove_phix
└── Taxonomy
    └── kraken2
        ├── <sample_name>
```


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

<details>
  <summary>Click here to see the Summary of the resulting data structure from Assembly folder</summary>

```ruby
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
```
</details>

## Results and statistics. 

From the files generated after the two analyses (viralrecon and metagenomics), the next aim is to extract and mount in a tab separated file different data from both analyses that will give an overview of the metrics and summary results. 


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
