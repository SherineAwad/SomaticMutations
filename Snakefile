import pandas as pd

df= pd.read_csv('inputs/samples.tsv', delimiter='\t')
dict = dict(df.values)
TUMORS =[] 
NORMALS = []
for key in dict.keys(): 
     TUMORS.append(key)
     NORMALS.append(dict[key])

print("dict is" , dict) 
print(TUMORS, NORMALS)
def return_pair(wildcards):
       tumor = ''.join({wildcards.sample})
       normal= dict[tumor]
       tumor += ".recalibrated.bam"
       normal += ".recalibrated.bam" 
       return {"tumor": tumor , "normal": normal} 

rule all: 
   input:
       expand("galore/{sample}_trimmed.fq.gz", sample = TUMORS), 
       expand("{sample}.sam", sample = TUMORS),
       expand("{sample}.RG.sam", sample = TUMORS),
       expand("{sample}.dedupped.bam", sample = TUMORS),
       expand("{sample}.split.bam", sample = TUMORS),
       expand("{sample}.recal_data.table", sample = TUMORS),
       expand("{sample}.recalibrated.bam", sample = TUMORS),
       #Same for Normals 
       expand("galore/{sample}_trimmed.fq.gz", sample = NORMALS),
       expand("{sample}.sam", sample = NORMALS),
       expand("{sample}.RG.sam", sample = NORMALS),
       expand("{sample}.dedupped.bam", sample = NORMALS),
       expand("{sample}.split.bam", sample = NORMALS),
       expand("{sample}.recal_data.table", sample = NORMALS),
       expand("{sample}.recalibrated.bam", sample = NORMALS),
       expand("{sample}_pon.vcf.gz", sample =NORMALS),
       #Create Panel of Normals 
       expand("{sample}_somatics.vcf.gz", sample =TUMORS)

rule trim:
    input:
       "{sample}.fq",
    output:
      "galore/{sample}_trimmed.fq.gz",
    conda: "env/env-trim.yaml"
    shell:
         """
          trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore {input} 
         """

rule tosam:
    input:
       "galore/{sample}_trimmed.fq.gz"
    output:
        "{sample}.sam"
    params: 
       genome = config['GENOME'], 
       mem = config['MEMORY']
    shell:
       """
       bbmap.sh {params.mem} in={input} out={output} ref={params.genome}
       """ 

rule AddRG:
    input:
       "{sample}.sam"
    output:
       "{sample}.RG.sam"
    params:
       RG = config['RG']
    conda: "env/env-picard.yaml"
    shell:
        """
        picard AddOrReplaceReadGroups I={input} O={output} SO=coordinate RGID=@{params} RGSM={wildcards.sample} RGPL=Illumina RGLB={wildcards.sample} RGPU={params}_{wildcards.sample} VALIDATION_STRINGENCY=SILENT
        """

rule deduplicate: 
    input: 
        "{sample}.RG.sam" 
    output: 
        "{sample}.dedupped.bam", 
        "{sample}.output.metrics"
    shell: 
        """
        picard MarkDuplicates I={input} O={output[0]}  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={output[1]} 
        """

rule split: 
    input: 
       "{sample}.dedupped.bam" 
    params: 
       genome = expand("{genome}.fa", genome= config['GENOME'])
    output: 
       "{sample}.split.bam" 
    shell: 
       """
       gatk SplitNCigarReads -R {params.genome} -I {input} -O {output} 
       """

rule recalibrate_a: 
    input:
       "{sample}.split.bam"
    params:
        genome = expand("{genome}.fa", genome= config['GENOME']),
        DBSNP = config['DBSNP'],
        INDELS = config['INDELS'],
        GOLD_STANDARD = config['GOLD_STANDARD'] 
    output:
       "{sample}.recal_data.table"
    shell: 
       """
       gatk BaseRecalibrator -I {input} -R {params.genome}  --known-sites {params.DBSNP} --known-sites {params.INDELS} --known-sites {params.GOLD_STANDARD} -O {output} 
       """

rule recalibrate_b: 
    input: 
       "{sample}.recal_data.table"
    params: 
       genome = expand("{genome}.fa", genome= config['GENOME']), 
    output: 
       "{sample}.recalibrated.bam"
    shell: 
      """ 
        gatk ApplyBQSR -I {input} -R {params.genome} --bqsr-recal-file {input} -O {output}
      """

rule Mutect2:
     input: 
         "{sample}.recalibrated.bam"
     params: 
        genome= expand("{genome}.fa", genome= config['GENOME']),
        AFONLYGNOMAD = config['AFONLYGNOMAD']
     output: 
         "{sample}_pon.vcf.gz"
     shell:
        """
        gatk Mutect2 -R {params.genome}  -I {input} --max-mnp-distance 0 --germline-resource {params.AFONLYGNOMAD} -O {output}
        """

rule PON_DB:
     input:
         sample = expand("{sample}_pon.vcf.gz", sample = NORMALS),
     params:
         genome= expand("{genome}.fa", genome= config['GENOME']),
         mem = {"-Xmx100g"},
         intervals = config['INTERVALS'],
         I =  lambda w: "-V " + " -V ".join(expand("{sample}_pon.vcf.gz", sample =NORMALS))
     output: 
        pon_db = config['PON_DB']
     shell:
         """
         gatk --java-options {params.mem}  GenomicsDBImport -R {params.genome} --genomicsdb-workspace-path {output.pon_db} -L {params.intervals}
         """

rule panel_normals:
    input: 
       pon_db = config['PON_DB']
    output:
       normal_panel= config['NORMALS_PANEL']
    params:
       AFONLYGNOMAD = config['AFONLYGNOMAD'],
       genome= expand("{genome}.fa", genome= config['GENOME']),
    shell:
       """
       gatk CreateSomaticPanelOfNormals -R {params.genome} --germline-resource {params.AFONLYGNOMAD} -V gendb://{input.pon_db} -O {output} 
       """
rule somatic_call: 
     input: 
         unpack(return_pair)
     params:
        genome = expand("{genome}.fa", genome= config['GENOME']),
        AFONLYGNOMAD = config['AFONLYGNOMAD'],
        pon=config['NORMALS_PANEL'],
        prefix = lambda wildcards: dict[wildcards.sample] 
     output: 
         "{sample}_somatics.vcf.gz"
     shell: 
        """
         gatk Mutect2 -R {params.genome} -I {input.tumor} -I {input.normal} -normal {params.prefix} --germline-resource {params.AFONLYGNOMAD}  --panel-of-normals {params.pon} -O {output}
        """ 

