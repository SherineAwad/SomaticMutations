[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.2-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


Snakemake Workflow for somatic mutation calling 
==========================================================================

The is GATK4/Mutect2 pipeline for Somatic Mutation. 

### Requirments 

- trim-galore=0.6.6
- star=2.7.10a
- picard=2.25.6
- gatk4=4.2.0.0

You can run the pipeline in -use-conda mode to pull these tools automatically. See use conde section below.

### Edit config file 

You will need to edit the *config file* to match your samples and paramteres. 

The pipeline expects samples with suffix ".r_1.fq.gz" and ".r_2.fq.gz" if samples are paired-end.
Any prefix before this suffix is the sample name and to be written in the "samples.tsv". For single-end reads, the samples suffix is ".fq.gz" and any prefix before this suffix is written in the **"samples.tsv"**.
For example, if your sample name is sample1.s_1.r_1.fq.gz, then your sample name in the samples file should be sample1.s_1.

You need to update the *config file* with whether your samples are paired-end or single reads. If your samples are paired-end, then the **PAIRD** entry in the config file should be set to TRUE, otherwise, set the **PAIRED** entry in the config file to FALSE. You can change the **samples.tsv** name in the *config file*.


The **samples.tsv** has the following format:

Tumors | Normals |
-------|---------|
SLX-18967.UDP0126.HT3G5DMXX.s_1| SLX-18967.UDP0129.HT3G5DMXX.s_1 |
SLX-18967.UDP0146.HT3G5DMXX.s_1| SLX-18967.UDP0149.HT3G5DMXX.s_1 |


You will need to edit the names and directory of your genome, your genome index, GTF, adapters, read groups in the **GENOME**, **INDEX**, **GTF**, **ADAPTERS**,  and **RG** entries in the *config file* respectively. 
You will also need to have your DBSNP vcf, indels vcf, gold standard vcf, and AF only gnomAD in the **DBSNP**, **INDELS**, **GOLD_STANDARD**, and **AFONLYGNOMAD** entries respectively  in the *config file*. If you are using human genome, these resources can be pulled from [Broad institute Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

You need to update your interval list, by editing the **intervals.list** file to list only the chromosomes of your interest. You can change the name of this file by editing the *config file* entry **INTERVALS**.
	
The pipeline will automatically pull the biallelic gnomAD. You can change its name/location by editing the **GNOMAD_BIALLELIC** entry in the *config file*. 


### Run Snakemake pipeline 

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


 
