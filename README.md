# Multimapper

[![License: GPL v3](https://img.shields.io/badge/license-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/bact-gen-scripts/blob/master/LICENSE)   


This pipeline is a re-write of the multiple_mappings_to_bam.py by Simon Harris, formerly of the Bacterial Genomics team at the Wellcome Sanger Institute.

This pipeline foccuses on

## Dependencies

- bwa  0.7.17
- gatk  3.7.0
- picard  1.126
- python 2.7
- raxml-ng  8.2.8
- smalt  0.7.6
- samtools  1.8


## Install

To install by Dockerfile please pull the docker image from the Pathogen Informatics repository[https://hub.docker.com/u/sangerpathogens]

Otherwise clone by running

```bash
git clone [githublinks_101]
```

## Basic Usage

To run the pipeline in its simplest form you need

## Surface Usage Options

The pipeline has a small group of options shared across it's 

## Pre-set configs

The pipeline can use alternate configs, accessible using -c <my_config> on the command line.

## Setting config guide

The nextflow.config file contained in this repo has the basic 

