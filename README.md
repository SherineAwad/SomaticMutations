Snakemake Workflow for somatic mutation calling from RNASeq Analysis 
==========================================================================


This is an Somatic mutation calling from RNASeq snakemake pipeline written by Sherine Awad. 
The pipeline so far uses GATK/Mutect2. More tools in the way. 

#### Edit config file 

The config file has several sections. You will need to edit the general section and add your tumor/normal paired samples in inputs/samples.tsv. 

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


 
