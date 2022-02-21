[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.2-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


Snakemake Workflow for somatic mutation calling 
==========================================================================

The is GATK4/Mutect2 pipeline for Somatic Mutation. 

#### Edit config file 

The config file has several sections. You will need to edit the general section and add your tumor/normal paired samples in samples.tsv.
You can change **samples.tsv** in the *config file*. AS Well as edit the config files for many other parameters.

The pipeline expects samples with suffix ".r_1.fq.gz" and ".r_2.fq.gz" if samples are paired-end.
Any prefix before this suffix is the sample name and to be written in the "samples.tsv". For single-end reads, the samples suffix is ".fq.gz" and any prefix before this suffix is written in the **"samples.tsv"**.
For example, if your sample name is sample1.s_1.r_1.fq.gz, then your sample name in the samples file should be sample1.s_1.

You need to update the *config file* with whether your samples are paired-end or single reads. If your samples are paired-end, then the **PAIRD** entry in the config file should be set to TRUE, otherwise, set the **PAIRED** entry in the config file to FALSE. You can change the **samples.tsv** name in the *config file*.


The **samples.tsv** has the following format:

.. list-table:: samples.tsv
   :widths: 35 35 
   :header-rows: 1

   * - Heading TUMORS
    -  Heading NORMALS
   * - Row 1, SLX-18967.UDP0126.HT3G5DMXX.s_1
     - Row 1, SLX-18967.UDP0129.HT3G5DMXX.s_1
   * - Row 2, SLX-18967.UDP0146.HT3G5DMXX.s_1
     - Row 2, SLX-18967.UDP0149.HT3G5DMXX.s_1

You will need to edit the names and directory of your genome, your genome index, GTF, adapters, read groups in the GENOME, INDEX, GTF, ADAPTERS,  and RG entries in the *config file* respectively. 
You will also need to have your DBSNP vcf, indels vc, gold standard vcfs, and AF only gnomAD in the DBSNP, INDELS, GOLD_STANDARD, and AFONLYGNOMAD entries respectively  in the *config file*. If you are using human genome, this can be pulled from broadinstitute Get the latest news at `Resource Bundle`_.

.. _Resource Bundle: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle


You need to update your interval list, by editing the **intervals.list** file to list only the chromosomes of interest. You can change the name of this file by editing the *config file* entry **INTERVALS**.
	
The pipeline will automatically pull the biallelic gnomAD. You can change its name/location by editing the GNOMAD_BIALLELIC entry in the *config file*. 


#### Run Snakemake pipeline 

Once you edit the config file to match your needs, then:  


    snakemake -jn 

where n is the number of cores for example for 10 cores use:


    snakemake -j10 


For a dry run use: 
  
  
    snakemake -j1 -n 


and to print command in dry run use: 

  
    snakemake -j1 -n -p 

  
#### Use Conda 

For less frooodiness, to pull automatically the same versions of dependencies use:

    snakemake -jn --use-conda

This will pull the same versions of tools we used. Conda has to be installed in your system.

For example, for 10 cores:

    snakemake -j10 --use-conda


#### Dry run 

for a dry run use:

    snakemake -j1 -n

and you can see the command printed on a dry run using:

    snakemake -j1 -n -p


#### Keep going option 


You can try the following to keep going if any issues happen, like no variants is found by one tool:

    snakemake -j1 --keep-going


#### Collect some stats 

    snakemake -j 10 --keep-going --stats run.stats


 
