from pathlib import Path
import random

MASON_BIN = Path(workflow.basedir, "mason2-2.0.9-Linux-x86_64/bin/mason_simulator")
CLAM = Path(workflow.basedir, "clam/target/release/clam")
SAMPLES = [f"tsk_{i}" for i in range(10)]

SEQLEN = 1000
COVERAGE = 10
NUM_READS = (SEQLEN * COVERAGE) // 300
SEED = int(config.get("seed", 1234))
logger.warning(f"seed = {SEED}")
NE=1e4
MU=2.5e-8
MIN_DEPTH = 2
WINDOW_SIZE = 1000

MAX_MISSING = 0.8
CLAM_LOCI_THREADS = 4
CLAM_STAT_THREADS = 4


def read_seed(wc):
    sample_number = int(wc.sample.split("_")[1])
    return SEED + sample_number



rule all:
    input:
        # expand("depths/{sample}.per-base.d4", sample=SAMPLES)
        "pixy_pi.txt",
        "clam_pi.tsv",
        # "msprime_pi_windows.tsv",
        # "vcftools_pi.txt"
        
rule simulate:
    output:
        "sim/all.vcf.gz",
        "sim/all.vcf.gz.tbi",
        "sim/sim_results.txt",
        "msprime_pi_windows.tsv",
        "sim/ref.fa"
    conda:
        "clam_benchmarking"
    threads: 4
    params:
        ne = NE,
        mu = MU,
        sample_size = len(SAMPLES),
        seed = SEED,
        outdir = "sim",
        seqlen = SEQLEN,
        windows=1000,
        indv_vcfs=False
    script:
        "scripts/simulate.py"

rule split_vcf:
    input:
        vcf = "sim/all.vcf.gz",
        idx = "sim/all.vcf.gz.tbi",
    output:
        "sim/{sample}.vcf"
    conda:
        "clam_benchmarking"
    shell:
        """
        bcftools view -s {wildcards.sample} {input.vcf} > {output}
        """

rule bwa_index:
    input:
        "sim/ref.fa"
    output:
        expand("sim/ref.fa.{index}", index=["bwt", "amb", "ann", "pac", "sa"])
    conda:
        "clam_benchmarking"
    shell:
        "bwa index {input}"

rule samtools_index:
    input:
        "sim/ref.fa"
    output:
        fai = "sim/ref.fa.fai",
        dictf = "sim/ref.dict"
    conda:
        "clam_benchmarking"
    shell:
        """
        samtools faidx {input}
        samtools dict {input} > {output.dictf}
        """

rule simulate_reads:
    input:
        mason = MASON_BIN,
        ref = "sim/ref.fa",
        vcf = "sim/{sample}.vcf"
    params:
        num_reads = NUM_READS,
        seed = read_seed,
        
    output:
        r1 = temp("reads/{sample}_R1.fq"),
        r2 = temp("reads/{sample}_R2.fq"),
        sam=temp("alignments/{sample}.sam")
    log:
        "logs/sim_reads/{sample}.txt"
    threads: 4
    shell:
        """
        {input.mason} \
        --seed {params.seed} \
        -ir {input.ref} \
        -iv {input.vcf} \
        -n {params.num_reads} \
        --illumina-read-length 150 \
        --num-threads {threads} \
        -o {output.r1} -or {output.r2} -oa {output.sam} &> {log}
        """

rule samtools_sort:
    input:
        sam = "alignments/{sample}.sam"
    output:
        bam="bams/{sample}.bam",
        bai="bams/{sample}.bam.bai"
    params:
        rg = r"'@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA'"
    conda:
        "clam_benchmarking"
    log:
        "logs/bwa/{sample}.txt"
    shell:
        "samtools addreplacerg -r {params.rg} {input.sam} 2> {log} | samtools sort - 2>> {log}| samtools view -b > {output.bam} 2>> {log} && samtools index {output.bam} 2>> {log}"


# rule bwa:
#     input:
#         ref = "sim/ref.fa",
#         ref_idx = rules.bwa_index.output,
#         r1= "reads/{sample}_R1.fq",
#         r2= "reads/{sample}_R2.fq"
#     params:
#         rg = r"'@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA'"
#     conda:
#         "clam_benchmarking"
#     log:
#         "logs/bwa/{sample}.txt"
#     threads: 4
#     output:
#         bam="bams/{sample}.bam",
#         bai="bams/{sample}.bam.bai"
#     shadow: "minimal"
#     shell:
#         "bwa mem -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort - 2>> {log}| samtools view -b > {output.bam} 2>> {log} && samtools index {output.bam} 2>> {log}"

rule mosdepth:
    input:
        bam = "bams/{sample}.bam",
        index = "bams/{sample}.bam.bai"
    output:
        "depths/{sample}.per-base.d4"
    params:
        prefix=lambda wc, output: output[0].replace(".per-base.d4", ""),
    conda:
        "clam_benchmarking"
    shadow: "minimal"
    
    shell:
        """
        mosdepth --d4 {params.prefix} {input.bam}
        """

rule bgzip_d4:
    input:
        "depths/{sample}.per-base.d4"
    output:
        "depths/{sample}.per-base.d4.gz",
        "depths/{sample}.per-base.d4.gz.gzi"
    shell:
        """
        bgzip --index {input}
        """

rule d4_fof:
    input:
        expand("depths/{sample}.per-base.d4.gz", sample=SAMPLES)
    output:
        "depths/fof.txt"
    shell:
        """
        ls depths/*.per-base.d4.gz > {output}
        """


rule clam_loci:
    input:
        clam_binary = CLAM,
        d4 ="depths/fof.txt",
        pops = "vcf/clam_pops.txt"
    output:
        "clam_loci/callable_sites.d4"
    params:
        depth = f"-m {MIN_DEPTH}",
        prefix = "clam_loci/callable_sites"
    threads: 1
    log: "logs/clam_loci/log.txt"
    benchmark:
        "benchmarks/clam_loci/benchmark.txt"
    shell:
        """
        {input.clam_binary} loci -p {input.pops} {params.depth} -t {threads} {input.d4} {params.prefix} 2> {log}
        """


rule bcftools_call:
    input:
        ref = "sim/ref.fa",
        fai = "sim/ref.fa.fai",
        dictf = "sim/ref.dict",
        bams = expand("bams/{sample}.bam", sample=SAMPLES)
    output:
        vcf = "vars/allsites.vcf", # for debugging
        vcfgz = "vars/allsites.vcf.gz",
        vcftbi = "vars/allsites.vcf.gz.tbi"
    conda:
        "clam_benchmarking"
    shell:
        """
        bcftools mpileup -Ou --annotate 'FMT/DP' -f {input.ref} {input.bams} | bcftools call -m -Ov -o {output.vcf}
        bgzip -c {output.vcf} > {output.vcfgz}
        tabix -p vcf {output.vcfgz}
        """

rule bcftools_filter:
    input:
        vcfgz = "vars/allsites.vcf.gz",
        vcftbi = "vars/allsites.vcf.gz.tbi"
    output:
        vcf = "vars/allsites.filtered.vcf", # for debugging
        vcfgz = "vars/allsites.filtered.vcf.gz",
        vcftbi = "vars/allsites.filtered.vcf.gz.tbi"
    conda:
        "clam_benchmarking"
    params:
        min_depth = MIN_DEPTH,
    shell:
        """
        bcftools view --exclude-types indels --max-alleles 2 {input.vcfgz} | bcftools +setGT - -- -t q -n . -e 'FMT/DP>{params.min_depth}' > {output.vcf}
        bgzip -c {output.vcf} > {output.vcfgz}
        tabix -p vcf {output.vcfgz}
        """

rule bcftools_vars_only:
    input:
        vcfgz = "vars/allsites.filtered.vcf.gz",
        vcftbi = "vars/allsites.filtered.vcf.gz.tbi"
    output:
        vcf = "vars/vars_only.filtered.vcf", # for debugging
        vcfgz = "vars/vars_only.filtered.vcf.gz",
        vcftbi = "vars/vars_only.filtered.vcf.gz.tbi"
    conda:
        "clam_benchmarking"
    shell:
        """
        bcftools view -e 'ALT="."' {input.vcfgz} > {output.vcf}
        bgzip -c {output.vcf} > {output.vcfgz}
        tabix -p vcf {output.vcfgz}
        """

rule write_pops_file:
    output: "vcf/pops.txt"
    run:
        lines = [f"{s}\tpop1\n" for s in SAMPLES]
        
        # Split the list into two populations
        mid = len(SAMPLES) // 2
        pop1 = lines[:mid]
        pop2 = [line.replace("\tpop1", "\tpop2") for line in lines[mid:]]
        
        # Write to the output file
        with open(output[0], "w") as f:
            f.writelines(pop1 + pop2)


# rule vcftools_pi:
#     input:
#         vcf = "vars/allsites.filtered.vcf.gz",
#         vcftbi = "vars/allsites.filtered.vcf.gz.tbi",
#         pops = "vcf/pops.txt"
#     params:
#         windows=1000
#     output:
#         "vcftools_pi.txt"
#     conda: "clam_benchmarking"
#     shadow: "minimal"
#     benchmark:
#         "benchmarks/vcftools_pi/benchmark.txt"
#     threads: 1
#     shell:
#         """
#         vcftools --gzvcf {input.vcf} --window-pi {params.windows} --stdout > {output}
#         """

rule pixy:
    input:
        vcf = "vars/allsites.filtered.vcf.gz",
        vcftbi = "vars/allsites.filtered.vcf.gz.tbi",
        pops = "vcf/pops.txt"
    params:
        windows=WINDOW_SIZE
    output:
        "pixy_pi.txt"
    conda: "pixy"
    # shadow: "minimal"
    benchmark:
        "benchmarks/pixy/benchmark.txt"
    threads: 1
    shell:
        "pixy --n_cores {threads} --stats pi dxy fst --fst_type hudson --vcf {input.vcf} --window_size {params.windows} --populations {input.pops}"

rule clam_stat_bgzf:
    input:
        clam_binary = CLAM,
        vcf="vars/vars_only.filtered.vcf.gz",
        vcftbi = "vars/vars_only.filtered.vcf.gz.tbi",
        callable_sites= "clam_loci/callable_sites.d4",
        pops = "vcf/clam_pops.txt"
    output:
        pi="clam_pi.tsv",
    params:
        outdir="clam_stat_bgzf",
        windows=WINDOW_SIZE
    log: "logs/clam_stat_bgzf/log.txt"
    benchmark:
        "benchmarks/clam_stat/benchmark.txt"
    threads: 1
    shell:
        """
        {input.clam_binary} stat -t {threads} -p {input.pops} -w {params.windows} {input.vcf} {input.callable_sites}
        """
