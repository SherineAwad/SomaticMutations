import pandas as pd

input = config['SAMPLES']
df= pd.read_csv(input, delimiter='\t')
dict = dict(df.values)
print("df is", df)
print("dict is", dict)
TUMORS =[] 
NORMALS = []
for key in dict.keys(): 
     TUMORS.append(key)
     NORMALS.append(dict[key])

print("Tumors are:", TUMORS)
print("Normals are:",  NORMALS)
def return_pair(wildcards):
       tumor = ''.join({wildcards.sample})
       normal= dict[tumor]
       tumor += ".recalibrated.bam"
       normal += ".recalibrated.bam" 
       return {"tumor": tumor , "normal": normal} 

rule all: 
   input:
       expand("{sample}_Aligned.out.sam", sample = TUMORS),
       expand("{sample}_Aligned.out.sam", sample = NORMALS),
       expand("{sample}.RG.sam", sample = TUMORS),
       expand("{sample}.RG.sam", sample = NORMALS),
       expand("{sample}.dedupped.bam", sample = TUMORS),
       expand("{sample}.dedupped.bam", sample = NORMALS),
       expand("{sample}.recal_data.table", sample = TUMORS),
       expand("{sample}.recal_data.table", sample = NORMALS),
       expand("{sample}.recalibrated.bam", sample = TUMORS),
       expand("{sample}.recalibrated.bam", sample = NORMALS),
       #Create Panel Of Normals
       expand("{panel}", panel=config['NORMALS_PANEL']),
       expand("{sample}_pon.vcf.gz", sample =NORMALS),
       #Make Somatic Calls 
       expand("{sample}_somatics.vcf.gz", sample =TUMORS),
       #Filter Mutect Call
       expand("{sample}_somatic_oncefiltered.vcf.gz", sample =TUMORS) 

if config['PAIRED']:
    rule trim:
       input:
           r1 = "{sample}.r_1.fq.gz",
           r2 = "{sample}.r_2.fq.gz"
       output:
           "galore/{sample}.r_1_val_1.fq.gz",
           "galore/{sample}.r_2_val_2.fq.gz"
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
           """
    rule align:
       input:
          r1 = "galore/{sample}.r_1_val_1.fq.gz",
          r2 = "galore/{sample}.r_2_val_2.fq.gz"
       output:
             "{sample}_Aligned.out.sam"
       params:
            threads = config['THREADS'],
            gtf = config['GTF'],
            prefix = "{sample}_",
            index = config['INDEX']
       shell:
            """
            STAR --genomeDir {params.index} --runThreadN {params.threads} --readFilesCommand zcat --readFilesIn {input.r1} {input.r2}  --outFileNamePrefix {params.prefix} --sjdbGTFfile {params.gtf}  --twopassMode Basic
            """
else:
     rule trim:
       input:
           "{sample}.fq.gz",
       output:
           "galore/{sample}_trimmed.fq.gz",
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore {input}
           """
     rule align:
         input:
             "galore/{sample}_trimmed.fq.gz"
         output:
             "{sample}_Aligned.out.sam"
         params:
            threads = config['THREADS'],
            gtf = config['GTF'],
            prefix = "{sample}_",
            index = config['INDEX']
         shell:
            """
            STAR --genomeDir {params.index} --runThreadN {params.threads} --readFilesCommand zcat --readFilesIn {input}  --outFileNamePrefix {params.prefix} --sjdbGTFfile {params.gtf}  --twopassMode Basic
            """

rule AddRG:
    input:
       "{sample}_Aligned.out.sam"
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


rule recalibrate_a: 
    input:
       "{sample}.dedupped.bam"
    params:
        genome= config['GENOME'],
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
       "{sample}.dedupped.bam", 
       "{sample}.recal_data.table"
    params: 
       genome= config['GENOME'], 
    output: 
       "{sample}.recalibrated.bam"
    shell: 
      """ 
        gatk ApplyBQSR -I {input[0]} -R {params.genome} --bqsr-recal-file {input[1]} --disable-sequence-dictionary-validation -O {output}
      """

rule Mutect2:
     input: 
         "{sample}.recalibrated.bam"
     params: 
        genome= config['GENOME'],
        AFONLYGNOMAD = config['AFONLYGNOMAD'],
        prefix = "{sample}"
     output: 
         "{sample}_pon.vcf.gz"
     shell:
        """
        gatk Mutect2 -R {params.genome}  -I {input} -tumor {params.prefix} --max-mnp-distance 0 --germline-resource {params.AFONLYGNOMAD} -O {output}
        """

rule PON_DB:
     input:
         sample = expand("{sample}_pon.vcf.gz", sample = NORMALS),
     params:
         genome= config['GENOME'],
         mem = {"-Xmx100g"},
         intervals = config['INTERVALS'],
         I =  lambda w: "-V " + " -V ".join(expand("{sample}_pon.vcf.gz", sample =NORMALS))
     output: 
        pon_db = config['PON_DB']
     shell:
         """
         gatk --java-options {params.mem}  GenomicsDBImport -R {params.genome} --genomicsdb-workspace-path {output.pon_db} -L {params.intervals} {params.I} 
         """

rule panel_normals:
    input: 
       pon_db = config['PON_DB']
    output:
       normal_panel= config['NORMALS_PANEL']
    params:
       AFONLYGNOMAD = config['AFONLYGNOMAD'],
       genome= config['GENOME'],
    shell:
       """
       gatk CreateSomaticPanelOfNormals -R {params.genome} --germline-resource {params.AFONLYGNOMAD} -V gendb://{input.pon_db} -O {output} 
       """
rule somatic_call: 
     input: 
         unpack(return_pair), 
         pon = config['NORMALS_PANEL'],
     params:
        genome= config['GENOME'],
        AFONLYGNOMAD = config['AFONLYGNOMAD'],
        #pon=config['NORMALS_PANEL'],
        prefix = lambda wildcards: dict[wildcards.sample] 
     output: 
         "{sample}_somatics.vcf.gz"
     shell: 
        """
         gatk Mutect2 -R {params.genome} -I {input.tumor} -I {input.normal} -normal {params.prefix} --germline-resource {params.AFONLYGNOMAD}  --panel-of-normals {input.pon} -O {output}
        """ 

rule SelectVariants:
   params: 
      gnomad = config['GNOMAD'], 
      genome = config['GENOME']
   output: 
       config['GNOMAD_BIALLELIC']
   shell: 
      """
       gatk SelectVariants -R {params.genome} -V {params.gnomad} --restrict-alleles-to BIALLELIC -O {output}
      """   

rule GetPileupSummaries: 
   input:  
       "{sample}.recalibrated.bam",
       config['GNOMAD_BIALLELIC']
   params: 
        config['INTERVALS'] 
   output: 
       "{sample}_getpileupsummaries.table"
   shell: 
       """
        gatk GetPileupSummaries -I {input[0]} -V {input[1]} -L {params} -O {output}
       """ 

rule CalculateContamination:
    input: 
       "{sample}_getpileupsummaries.table"
    output:
       "{sample}_tumor_calculatecontamination.table"
    shell: 
        """
        gatk CalculateContamination -I {input} -O {output}
        """

rule FilterMutectCalls: 
    input: 
      "{sample}_tumor_calculatecontamination.table",
      "{sample}_somatics.vcf.gz"
    params:
      genome= config['GENOME'] 
    output: 
       "{sample}_somatic_oncefiltered.vcf.gz"
    shell: 
      """
       gatk FilterMutectCalls -R {params} -V {input[1]} --contamination-table {input[0]} -O {output} 
      """

