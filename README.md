<p align="left">
    <img src="docs/images/loma_logo.png" alt="loma logo" width="30%">
</p>

***

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.4-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)

# Table of contents 
* [Introduction](#Introduction)
* [Summary](#summary)
* [Installation](#install)
* [Running](#run)
* [Output](#output)
* [Tips for improving speed and efficiency](#tips)
* [Troubleshooting and errors](#troubleshoot)

# Introduction <a name="Introduction"></a>

LOMA is a Nextflow pipeline designed to comprehensively assess metagenomic samples sequenced using the Oxford Nanopore (long-read) platform.

The pipeline has two primary approaches to analysis:\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1). Read-based taxonomic classification - assign a taxonomic designation to individual sequencing reads to assess the sample composition.\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2). Assembly of metagenome assembled genomes (MAGs) followed by general/taxonomy-specific *in silico* typing.

The results of these analyses are then available of per-metric files and further compilied into a series of summary reports (HTML and TSV files).

A general overview is provided below. Detailed guidance on the installation, usage and function of LOMA can be found in the **[`wiki`]()**, example outputs can be found here.
# Pipeline summary <a name="summary"></a>

## Simplified schematic overview 
<p align="center">
    <img src="docs/images/loma_schematic_simplified.png" alt="loma overview" width="90%">
</p>

## Description

The pipeline will perform the following steps: 

**1). Read quality control** - Assess read quality, remove adapters and filter long-reads by quality.\
**2). Read-based taxonomic annotation** - Read-based taxonomic classification, standardization and summary reporting.\
**3). Host read removal** - Identify and remove host-contaminant reads followed by merging of quality control results into summary reports.\
**4). Assembly** - Assembly of reads into metagenome and polishing (contig error correction).\
**5). Contig analysis** - Per-contig identification of closest taxonomic hits, identification of mobile genetic elements and calculation of contig summary statistics.\
**6). Assembly binning** - Classify and bin contigs into individual metagenome assembled genomes (MAGs).\
**7). Bin quality control** - Assess the quality of MAGs and merge bin QC and contig QC results into summary reports.\
**8). Typing** - Subset MAGs of interest (target species) and dependent on taxonomic classification pass them on to individual subworkflows (run per-MAG).\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**8a). Bacteria** - Identification of genes of interest, multi-locus sequence type and screen for plasmids.\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**8b). *Listeria monocytogenes*** - Perform *in silico* serogroup prediction.\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**8c). *Salmonella*** - Perform *in silico* *Salmonella* serotyping, identify cgMLST alleles, AMR genes and  lineages (*S. Typhi* and *S. Paratyphi B* only).\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**8d). *Escherichia coli* / *Shigella spp.*** - Identify pathotype, serotype, AMR genes and lineage (*S. sonnei*) of *Escherichia coli* / *Shigella spp.*\
**9). Antimicrobial resistance** - Identify AMR genes (incl. point mutations), virulence/stress resistance genes and merge results into a summary report.

# Installation <a name="install"></a>
LOMA has two dependencies:

- A container runtime, both Apptainer and Singularity are supported.
- Nextflow.

Only one database is mandatory to run LOMA - the host reference genome assembly or the host reference Kraken2 database (can be either or both). There are 14 optional databases required to run certain stages if these are not downloaded/installed then parts of the pipeline will be skipped. To simplify database installation, a [`script`](https://github.com/ukhsa-collaboration/LOMA/blob/main/bin/get_dbs.py) is provided which will download any requested databases and update the relevant config files. 

Detailed installation instructions for LOMA and associated databases can be found on the **[`wiki`]()**.

# Running  <a name="run"></a>

### Usage

There is only one mandatory parameter for running LOMA, an input file (format detailed below).

```
./run_loma --input input.tsv
```

### Input file structure

The input file (e.g. 'input.tsv') is a five column tab-separated file with the following structure:
```
RUN_ID  BARCODE_ID  SAMPLE_ID   SAMPLE_TYPE /FULL/PATH/TO/FASTQ_FILE
```

- RUN_ID:                       Run identifier, will determine the highest level directory name in the results directory
- BARCODE_ID:                   Sample barcode 
- SAMPLE_ID:                    Sample identifier, will determine the subdirectory where results are stored per-sample
- SAMPLE_TYPE:                  Sample description, will be added to the reports, but doesn't change how the sample is processed. 
- /FULL/PATH/TO/FASTQ_FILE:     Location of input FASTQ files.

Any number of samples can be included provided they do not have both identical RUN_ID and SAMPLE_ID's (either is fine though). If any of the columns contain a period ('.'), they'll be automatically replaced with an underscore ('_') in the output. 

Example input file:
```
RUN01	RB01	SAMPLE_1	BLOOD	/data/projects/metagenome_ont/SAMPLE_1.BLOOD.fq.gz
RUN01	RB02	SAMPLE_2	BLOOD	/data/projects/metagenome_ont/SAMPLE_2.BLOOD.fq.gz
RUN02	RB01	SAMPLE_3	SALIVA	/data/projects/metagenome_ont/SAMPLE_3.NASOPHARYNGEAL.fq.gz
RUN03	XBD     SAMPLE_1	SKIN	/data/projects/metagenome_ont/SAMPLE_3.SKIN.fq.gz

```
Further examples can be found [`here`](). 

Detailed running instructions for LOMA can be found on the **[`wiki`]()**.

### Optional parameters





Detailed running instructions for LOMA can be found on the **[`wiki`]()**.



# Output  <a name="output"></a>

## Output folder structure

LOMA will output files per-metric/tool as well as a series of summary reports. Outputs can be found in the directory specified by the '--outdir' parameter (default: 'results'). Results are separated by RUN_ID, then SAMPLE_ID and finally by approach and metric. A summary of the folder structure can be found below. 

Detailed descriptions of the output folder structure and summary reportscan be found on the **[`wiki`]()**.


## HTML reports

A series of HTML summary reports can be found in the 'summary' directory, these cover the major areas of interest when sequencing metagenomic samples and should assist in assessing saequencing data quality, sample composition and perform various types of *in silico* phenotyping.


#### **<SAMPLE_ID>.<RUN_ID>.summary_report.html** ([`example report`]())
- Simiplified summary detailing all major metrics/results of interest, for an overall sample summary.

#### **<SAMPLE_ID>.<RUN_ID>.readqc_report.html** ([`example report`]())
- Read quality metrics, pre- and post- quality control.

#### **<SAMPLE_ID>.<RUN_ID>.taxonomy_report.html** ([`example report`]())
- Read-based taxonomic abundance for Kraken2, Centrifuger and/or Sylph (depending on which tools were used).

####  **<SAMPLE_ID>.<RUN_ID>.amr_report.html** ([`example report`]())
- Results of AMR typing tools (ABRicate, AMRFinderPlus, ResFinder and RGI).

#### **<SAMPLE_ID>.<RUN_ID>.summary_binning_report.html** ([`example report`]())
- Summary of binning results.

# Tips for improving speed and efficiency  <a name="tips"></a>

### Skipping analysis steps

When specified, the following parameters will skip substantial sections of the pipeline, saving resources if the results are not of interest:

```
  --skip_assembly                                 Skip read assembly.
  --skip_taxonomic_profiling                      Skip read-based taxonomic profiling.
  --skip_prokarya_typing                          Skip metagenome assembled genome analyses.
```  

### Skipping Read-based taxonomic annotation


```
  TAXONOMIC_PROFILING.krakendb = "/data/databases/kraken2_databases/kraken2_standardPlusPF_57_100923/"
  TAXONOMIC_PROFILING.centrifugerdb = "/data/databases/centrifuger/"
  TAXONOMIC_PROFILING.sylphdb = "/data/databases/gtdb_tk_databases/sylph_db/gtdb-r220-c200-dbv1.syldb"
```



### Skipping polishing

Assembly error correction is a time consuming step. To save time you can reduce the number of rounds of Racon polishing (default: 4, range: 0-4), e.g.:
```
--ASSEMBLY.racon_rounds 1
```
Medaka polishing is very slow and is disabled by default it can be enabled by specifying: 
```
--ASSEMBLY.medaka
```
If you find the per-base accuracy of your MAGs are low, even after Racon polishing but you want to skip running Medaka on the entire metagenomic assembly, a quicker approach is to polish only the MAGs of interest. This can be done by specifying:
```
--BIN_TAXONOMY.medaka_mag
```

Further tips for optimization can be found on the **[`wiki`]()**.


# Troubleshooting and errors  <a name="troubleshoot"></a>
