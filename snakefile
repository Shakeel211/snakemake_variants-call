# --- Snakefile ---

configfile: "config.yaml"

# Get all samples by finding R1 files

def get_samples():
    import glob
    import os
    
    r1_files = glob.glob("data/raw/*_R1_*.fastq.gz")
    samples = []
    
    for f in r1_files:
        filename = os.path.basename(f)
        
        # For files like: BC-GP-S46_S46_L001_R1_001.fastq.gz
        # Extract everything before '_R1_001.fastq.gz'
        if '_R1_001.fastq.gz' in filename:
            sample = filename.replace('_R1_001.fastq.gz', '')
        # Or if the pattern is different, use this more flexible approach:
        elif '_R1_' in filename:
            sample = filename.split('_R1_')[0]
        else:
            # Fallback to original method
            sample = filename.split('_')[0]
            
        samples.append(sample)
        print(f"File: {filename} → Sample: {sample}")  # Debug output
    
    unique_samples = list(set(samples))
    print(f"Final detected samples: {unique_samples}")
    return unique_samples

SAMPLES = get_samples()


# Rule to run all samples through the final variant calling step
rule all:
    input:
        expand("results/variants/{sample}.vcf.gz", sample=SAMPLES)


# 1. FASTQC - Quality Control
rule fastqc_raw:
    input:
        r1 = "data/raw/{sample}_R1_001.fastq.gz",
        r2 = "data/raw/{sample}_R2_001.fastq.gz"
    output:
        html = ["results/fastqc/{sample}_R1_fastqc.html", 
                "results/fastqc/{sample}_R2_fastqc.html"],
        zip = ["results/fastqc/{sample}_R1_fastqc.zip", 
               "results/fastqc/{sample}_R2_fastqc.zip"]
    log:
        "logs/fastqc/{sample}.log"
    threads: 2
    shell:
        """
        fastqc {input.r1} {input.r2} \
            -o results/fastqc/ \
            -t {threads} \
            > {log} 2>&1
        """


# 2. FASTP - Quality Trimming

rule fastp_trim:
    input:
        r1 = "data/raw/{sample}_R1_001.fastq.gz",
        r2 = "data/raw/{sample}_R2_001.fastq.gz"
    output:
        r1 = "results/trimmed/{sample}_R1_trimmed.fastq.gz",
        r2 = "results/trimmed/{sample}_R2_trimmed.fastq.gz",
        html = "results/reports/fastp/{sample}.html",
        json = "results/reports/fastp/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    threads: 4
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              -h {output.html} -j {output.json} \
              --thread {threads} \
              --detect_adapter_for_pe \
              --correction \
              --cut_right \
              --cut_window_size 4 \
              --cut_mean_quality 20 \
              --length_required 36 \
              2> {log}
        """

# 3. FASTQC on trimmed reads
rule fastqc_trimmed:
    input:
        r1 = "results/trimmed/{sample}_R1_trimmed.fastq.gz",
        r2 = "results/trimmed/{sample}_R2_trimmed.fastq.gz"
    output:
        html = ["results/fastqc/trimmed/{sample}_R1_trimmed_fastqc.html", 
                "results/fastqc/trimmed/{sample}_R2_trimmed_fastqc.html"],
        zip = ["results/fastqc/trimmed/{sample}_R1_trimmed_fastqc.zip", 
               "results/fastqc/trimmed/{sample}_R2_trimmed_fastqc.zip"]
    log:
        "logs/fastqc/trimmed/{sample}.log"
    threads: 2
    shell:
        """
        fastqc {input.r1} {input.r2} \
            -o results/fastqc/trimmed/ \
            -t {threads} \
            > {log} 2>&1
        """

# 4. BWA-MEM - Alignment
rule bwa_mem:
    input:
        r1 = "results/trimmed/{sample}_R1_trimmed.fastq.gz",
        r2 = "results/trimmed/{sample}_R2_trimmed.fastq.gz",
        ref = config["reference_fasta"]
    output:
        "results/aligned/{sample}.bam"
    log:
        "logs/bwa/{sample}.log"
    params:
        read_group = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"
    threads: 8
    shell:
        """
        bwa mem -t {threads} \
            -R '{params.read_group}' \
            {input.ref} {input.r1} {input.r2} \
            2> {log} | \
        samtools view -@ 4 -Sb - > {output} 2>> {log}
        """

# 5. Sort BAM file
rule samtools_sort:
    input:
        "results/aligned/{sample}.bam"
    output:
        "results/aligned/{sample}.sorted.bam"
    log:
        "logs/samtools_sort/{sample}.log"
    threads: 4
    shell:
        """
        samtools sort -@ {threads} \
            -o {output} {input} \
            2> {log}
        """


# 6. Mark Duplicates with Picard
rule picard_markduplicates:
    input:
        "results/aligned/{sample}.sorted.bam"
    output:
        bam = "results/dedup/{sample}.dedup.bam",
        metrics = "results/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/{sample}.log"
    conda:
        "envs/picard.yaml"
    shell:
        """
        picard MarkDuplicates \
            I={input} \
            O={output.bam} \
            M={output.metrics} \
            REMOVE_DUPLICATES=false \
            CREATE_INDEX=true \
            > {log} 2>&1
        """

# 7. GATK Base Quality Score Recalibration (BQSR)

rule gatk_baserecalibrator:
    input:
        bam = "results/dedup/{sample}.dedup.bam",
        ref = config["reference_fasta"],
        known_sites1 = config["known_sites1"],
        known_sites2 = config["known_sites2"]
    output:
        table = "results/bqsr/{sample}.recal.table"
    log:
        "logs/gatk/baserecalibrator/{sample}.log"
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {input.ref} \
            --known-sites {input.known_sites1} \
            --known-sites {input.known_sites2} \
            -O {output.table} \
            > {log} 2>&1
        """

rule gatk_applybqsr:
    input:
        bam = "results/dedup/{sample}.dedup.bam",
        table = "results/bqsr/{sample}.recal.table",
        ref = config["reference_fasta"]
    output:
        "results/bqsr/{sample}.recal.bam"
    log:
        "logs/gatk/applybqsr/{sample}.log"
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk ApplyBQSR \
            -I {input.bam} \
            -R {input.ref} \
            --bqsr-recal-file {input.table} \
            -O {output} \
            > {log} 2>&1
        """


# 8. GATK HaplotypeCaller - Variant Calling

rule gatk_haplotypecaller:
    input:
        bam = "results/bqsr/{sample}.recal.bam",
        ref = config["reference_fasta"]
    output:
        "results/variants/{sample}.vcf.gz"
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        extra = config.get("gatk_extra_params", "")
    conda:
        "envs/gatk.yaml"
    shell:
        """
        gatk HaplotypeCaller \
            -I {input.bam} \
            -R {input.ref} \
            -O {output} \
            {params.extra} \
            > {log} 2>&1
        """
