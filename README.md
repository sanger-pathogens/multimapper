# Multimapper

[![License: GPL v3](https://img.shields.io/badge/license-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/bact-gen-scripts/blob/master/LICENSE)   


This pipeline is a re-write of the multiple_mappings_to_bam.py by Simon Harris, formerly of the Bacterial Genomics team at the Wellcome Sanger Institute.

This pipeline focuses on being a quick and easy mapping software. It can run in a directory of paired lanes. So long as the pairs are identicle in prefix they can be paired.

## Dependencies

- bwa-mem2
- python 3.7+
- raxml-ng 
- smalt  
- samtools  1.6+
- bowtie2
- 


## Install

To install by Dockerfile please pull the docker image from the Pathogen Informatics repository[https://hub.docker.com/u/sangerpathogens]

Otherwise clone by running the following and installing all listed dependencies.

```bash
git clone https://github.com/sanger-pathogens/multimapper.git
cd multimapper
pip3 install requirements.txt
```

## Basic Usage

To run the pipeline in its simplest form you need a reference to map to and a folder containing paired reads.

## Surface Usage Options

The pipeline has a small group of options shared across it's mappers and vcf tools. These include 
```commandline
--qc      if set to 'on' will produce qc readings for mapped reads
--mapper  sets the mapper from between 'smalt', 'bwa-mem2' and 'bowtie2'
--snp     either sh16 for the original conversions to bcf, or freebayes for a simpler variant calling
--k       scanning kmer length for mapping and snp calling processes
WIP: not fully functional
--pseudo  creation of pseudosequence from 
```

## Pre-set configs

The pipeline can use alternate configs, accessible using -c <my_config> on the command line. The available configs included are fast.config and a WIP full.config. 
Fast config removes qc and all Simon Harris processes to quickly process the data through BWA-mem2 and Freebayes. WIP Full config applies all available processing; this includes QC filtering, mapping, variant calling using Simon Harris methods and pseudosequence generation, indel trimming and phylogeny tree generation.

It can be run as shown below:

```commandline
nextflow multimapper.nf -c <fast/full>.config -r <reference.fa> 
```

## Setting config guide

The nextflow.config file contained in this repo has the basic set-up showing the available tools. To build a custom run, download the nextflow.config file and edit it so suit specific pipeline requirement.

### Links to external software

SMALT: documentation available near the bottom of this[https://www.sanger.ac.uk/tool/smalt-0/] page
