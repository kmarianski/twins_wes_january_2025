import os
from pathlib import Path

def get_final_samples():
    # Get unique sample IDs (everything before the first underscore)
    return list(set(s.split('_')[0] for s in get_fastq_stems()))

# Configuration
configfile: "config.yaml"

# Helper function to get complete file stems
def get_fastq_stems():
    stems = set()
    # Get stems from raw files
    raw_files = Path("raw").glob("*_1.fq.gz")  # Look for R1 files
    for f in raw_files:
        stem = f.name.replace('_1.fq.gz', '')
        stems.add(stem)
    return list(stems)

def get_sample_lanes():
    sample_dict = {}
    # Use fastq files instead of BAM files to determine samples
    for fq in Path("raw").glob("*_1.fq.gz"):
        # Assuming format: sampleID_laneInfo_1.fq.gz
        sample_id = fq.name.split('_')[0]
        stem = fq.name.replace('_1.fq.gz', '')
        if sample_id not in sample_dict:
            sample_dict[sample_id] = []
        sample_dict[sample_id].append(stem)
    return sample_dict
    
    # Debug print
    print(f"Found samples: {list(sample_dict.keys())}")
    return sample_dict

# Rule all specifies the final targets of the pipeline
rule all:
    input:
        # FastQC outputs for raw reads
        expand("results/fastqc/raw/{stem}_{read}_fastqc.html", 
               stem=get_fastq_stems(),
               read=["1", "2"]),
        # MultiQC reports for raw reads
        "results/multiqc/multiqc_raw_report.html",
        # Trimmomatic outputs
        expand("results/trimmomatic/{stem}_{read}P.fq.gz", 
               stem=get_fastq_stems(),
               read=["1", "2"]),
        # FastQC outputs for trimmed reads
        expand("results/fastqc/trimmed/{stem}_{read}P_fastqc.html",
               stem=get_fastq_stems(),
               read=["1", "2"]),
        # MultiQC reports for trimmed reads
        "results/multiqc/multiqc_trimmed_report.html",
        # BWA mapping outputs
        expand("results/bwa/{stem}.bam", stem=get_fastq_stems()),
        # Final merged BAMs
        expand("results/merged_bams/{sample}.bam", 
               sample=get_sample_lanes().keys()),
        # Index recalibrated BAMs
        expand("results/recal/{sample}.bam.bai", sample=get_sample_lanes().keys()), 
        # Deepvariant
        expand("results/deepvariant/{sample}.vcf.gz", sample=get_sample_lanes().keys()),
        # ANNOVAR outputs
        expand('results/annovar/{sample}.Info.recode.vcf', sample=get_sample_lanes().keys()), #annovar_recode
        expand('results/annovar/{sample}.avinput', sample=get_sample_lanes().keys()), #convert2annovar
        expand('results/annovar/multianno_csv/{sample}.hg38_multianno.csv', sample=get_sample_lanes().keys()), #annotate_variation

# Directory creation rule
rule setup_directories:
    output:
        touch("logs/setup_directories.done")
    run:
        for subdir in config["directories"]:
            Path("results/" + subdir).mkdir(parents=True, exist_ok=True)
        for subdir in config["log_dirs"]:
            Path("logs/" + subdir).mkdir(parents=True, exist_ok=True)
        shell("touch {output}")

# FastQC for raw reads
rule fastqc_raw:
    input:
        "raw/{stem}_{read}.fq.gz"
    output:
        html="results/fastqc/raw/{stem}_{read}_fastqc.html",
        zip="results/fastqc/raw/{stem}_{read}_fastqc.zip"
    log:
        "logs/fastqc/raw/{stem}_{read}.log"
    threads: 1
    resources:
        mem_mb=4000,
        time_min=30,
        cpus=1,
        job_name="fastqc_raw",
        partition="short"
    conda:
        "envs/biotools.yml"
    shell:
        "fastqc --threads {threads} {input} -o results/fastqc/raw &> {log}"

# MultiQC for raw reads
rule multiqc_raw:
    input:
        expand("results/fastqc/raw/{stem}_{read}_fastqc.html",
               stem=get_fastq_stems(),
               read=["1", "2"])
    output:
        "results/multiqc/multiqc_raw_report.html"
    log:
        "logs/multiqc/multiqc_raw.log"
    threads: 1
    resources:
        mem_mb=4000,
        time_min=30,
        cpus=1,
        job_name="multiqc_raw",
        partition="short"
    conda:
        "envs/biotools.yml"
    shell:
        "multiqc results/fastqc/raw -o results/multiqc --filename multiqc_raw_report.html &> {log}"

# Trimmomatic
rule trimmomatic:
    input:
        r1="raw/{stem}_1.fq.gz",
        r2="raw/{stem}_2.fq.gz"
    output:
        r1_paired="results/trimmomatic/{stem}_1P.fq.gz",
        r2_paired="results/trimmomatic/{stem}_2P.fq.gz",
        r1_unpaired="results/trimmomatic/{stem}_1U.fq.gz",
        r2_unpaired="results/trimmomatic/{stem}_2U.fq.gz"
    log:
        "logs/trimmomatic/{stem}.log"
    threads: 8
    resources:
        mem_mb=16000,
        time_min=120,
        cpus=8,
        job_name="trimmomatic",
        partition="short"
    params:
        adapter=config["adapter"]
    conda:
        "envs/biotools.yml"
    shell:
        "trimmomatic PE -threads {threads} {input.r1} {input.r2} "
        "{output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} "
        "ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 "
        "&> {log}"

# FastQC for trimmed reads
rule fastqc_trimmed:
    input:
        "results/trimmomatic/{stem}_{read}P.fq.gz"
    output:
        html="results/fastqc/trimmed/{stem}_{read}P_fastqc.html",
        zip="results/fastqc/trimmed/{stem}_{read}P_fastqc.zip"
    log:
        "logs/fastqc/trimmed/{stem}_{read}P.log"
    threads: 1
    resources:
        mem_mb=4000,
        time_min=30,
        cpus=1,
        job_name="fastqc_trimmed",
        partition="short"
    conda:
        "envs/biotools.yml"
    shell:
        "fastqc --threads {threads} {input} -o results/fastqc/trimmed &> {log}"

# MultiQC for trimmed reads
rule multiqc_trimmed:
    input:
        expand("results/fastqc/trimmed/{stem}_{read}P_fastqc.html",
               stem=get_fastq_stems(),
               read=["1", "2"])
    output:
        "results/multiqc/multiqc_trimmed_report.html"
    log:
        "logs/multiqc/multiqc_trimmed.log"
    threads: 1
    resources:
        mem_mb=4000,
        time_min=30,
        cpus=1,
        job_name="multiqc_trimmed",
        partition="short"
    conda:
        "envs/biotools.yml"
    shell:
        "multiqc results/fastqc/trimmed -o results/multiqc --filename multiqc_trimmed_report.html &> {log}"

# BWA Mapping
rule bwa_map:
    input:
        r1=rules.trimmomatic.output.r1_paired,
        r2=rules.trimmomatic.output.r2_paired
    output:
        bam="results/bwa/{stem}.bam"
    params:
        index=config["reference"]
    threads:
        config["threads"]
    log:
        "logs/bwa/{stem}.log"
    # conda:
    #     "envs/biotools.yml"
    resources:
        mem_mb=32000,
        time_min=240,
        cpus=16,
        job_name='bwa',
        partition="medium"
    shell:
        "bwa mem -t {threads} {params.index} {input.r1} {input.r2} | "
        "samtools sort -@{threads} -o {output.bam}"

# Samtools Merge
rule samtools_merge:
    input:
        # Get all BAM files for a given sample
        lambda wildcards: expand("results/bwa/{stem}.bam", 
                               stem=get_sample_lanes()[wildcards.sample])
    output:
        "results/merged_bams/{sample}.bam"
    # conda:
    #     "envs/biotools.yml"
    threads: 8
    resources:
        mem_mb=16000,
        time_min=120,
        cpus=8,
        job_name="samtools_merge",
        partition="short"
    log:
        "logs/merged_bams/{sample}.log"
    shell:
        """
        if [ $(echo {input} | wc -w) -gt 1 ]; then
            samtools merge -@ {threads} {output} {input} 2> {log}
        else
            ln -sr {input} {output} 2> {log}
        fi
        """

# Replace Read Groups
rule replace_rg:
    input:
        "results/merged_bams/{sample}.bam"
    output:
        "results/replace_rg/{sample}.bam"
    params:
        extra="--RGLB lib1 --RGPL illumina --RGPU {s} --RGSM {s}"
    log:
        "logs/picard/replace_rg/{s}.log"
    wrapper:
        "v1.14.1/bio/picard/addorreplacereadgroups"

# Mark Duplicates
rule markduplicates_bam:
    input:
        bams="results/replace_rg/{sample}.bam"
    output:
        bam="results/dedup/{sample}.bam",
        metrics="results/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        extra="--REMOVE_DUPLICATES true"
    wrapper:
        "v3.3.3/bio/picard/markduplicates"

# Base Recalibrator
rule base_recalibrator:
    input:
        bam="results/dedup/{sample}.bam",
        ref=config["reference"],
        known_sites=config["known_sites"]
    output:
        recal_table="results/recal/{sample}.table"
    log:
        "logs/picard/recal_table/{sample}.log"
    conda:
        "envs/biotools.yml"
    shell:
        "gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known_sites} "
        "-O {output.recal_table}"

# Apply BQSR
rule apply_bqsr:
    input:
        bam="results/dedup/{sample}.bam",
        recal_table="results/recal/{sample}.table"
    output:
        bam="results/recal/{sample}.bam"
    log:
        "logs/picard/recal_apply/{sample}.log"
    conda:
        "envs/biotools.yml"
    shell:
        "gatk ApplyBQSR -I {input.bam} --bqsr-recal-file {input.recal_table} -O {output.bam}"

# Index recalibrated BAMs
rule index_recal_bam:
    input:
        bam="results/recal/{sample}.bam"
    output:
        bai="results/recal/{sample}.bam.bai"
    conda:
        "envs/biotools.yml"
    log:
        "logs/samtools/index_recal/{sample}.log"
    shell:
        "samtools index -b {input.bam} > {log} 2>&1"

rule deepvariant:
    input:
        bam="results/merged_bams/{sample}.bam",
        bai="results/merged_bams/{sample}.bam.bai",
        ref=config["reference"]
    output:
        vcf="results/deepvariant/{sample}.vcf.gz"
    params:
        model="WES"
    threads: 
        16
    log:
        log="snakemake/logs/{sample}/"
    singularity:
        "singularity/deepvariant_latest.sif"
    shell:
        """
        /opt/deepvariant/bin/run_deepvariant \
            --model_type {params.model} \
            --ref {input.ref} \
            --reads {input.bam} \
            --output_vcf {output.vcf} \
            --num_shards {threads} \
            --logging_dir {log.log} \
            --make_examples_extra_args='vsc_min_count_snps=3,vsc_min_fraction_snps=0.2,vsc_min_count_indels=3,vsc_min_fraction_indels=0.10'
        """

rule annovar_recode:
    input:
        "results/deepvariant/{sample}.vcf.gz"
    output:
        recode='results/annovar/{sample}.Info.recode.vcf',
    params:
        param='results/annovar/{sample}.Info'
    conda:
        'envs/biotools.yml'        
    shell:
        'vcftools --gzvcf {input} --recode --out {params.param}'

rule convert2annovar:
    input:
        rules.annovar_recode.output.recode
    output:
        'results/annovar/{sample}.avinput'
    params:
        out='results/annovar/{sample}',
        annovar=config["annovar"]
    shell:
        'perl {params.annovar}/convert2annovar.pl -format vcf4 {input} > {output} -allsample -includeinfo -withfreq'

rule annotate_variation:
    input:
        rules.convert2annovar.output
    output:
        'results/annovar/multianno_csv/{sample}.hg38_multianno.csv'
    params:
        out='results/annovar/multianno_csv/{sample}',
        annovar=config["annovar"]
    threads:
        32
    shell:
        "perl {params.annovar}/table_annovar.pl -thread {threads} {input} {params.annovar}/humandb/ -buildver hg38 -out {params.out} -otherinfo -remove -protocol refGene,dbnsfp47a_interpro,gnomad41_genome,dbnsfp47a,clinvar_20240611 -operation g,f,f,f,f -nastring . -csvout"