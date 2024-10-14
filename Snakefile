MASON_BIN = "../mason2-2.0.9-Linux-x86_64/bin/mason_simulator"
MOPRS = "../mop-rs/target/release/moprs"
SAMPLES = [f"sample_{i}" for i in range(100)]
SEQLEN = 10_000_000
COVERAGE = 100
NUM_READS = (SEQLEN * COVERAGE) // 300
SEED = 1234
NE=1e6
MU=1e-8
MIN_DEPTH = 10
MAX_DEPTH = 10000000
MIN_DEPTH_MEAN = 20
MAX_DEPTH_MEAN = 500
MAX_MISSING = 0.8
MOPRS_THREADS = 4

rule all:
    input:
        "pixy_pi.txt",
        "callable_sites/moprs.bed",
        
rule simulate:
    output:
        temp(expand("sim/{sample}.vcf", sample=SAMPLES)),
        "sim/sim_results.txt",
        "sim/pi_windows.txt",
        "sim/ref.fa"
    conda:
        "envs/env.yaml"
    params:
        ne = NE,
        mu = MU,
        sample_size = len(SAMPLES),
        seed = SEED,
        outdir = "sim",
        seqlen = SEQLEN,
        windows=1000
    script:
        "scripts/simulate.py"

rule bwa_index:
    input:
        "sim/ref.fa"
    output:
        "sim/ref.fa.bwt"
    conda:
        "envs/env.yaml"
    shell:
        "bwa index {input}"


rule simulate_reads:
    input:
        mason = MASON_BIN,
        ref = "sim/ref.fa",
        vcf = "sim/{sample}.vcf"
    params:
        num_reads = NUM_READS,
        seed = SEED,
        
    output:
        r1 = temp("reads/{sample}_R1.fq"),
        r2 = temp("reads/{sample}_R2.fq")
    log:
        "logs/sim_reads/{sample}.txt"
    threads: 6
    shell:
        """
        {input.mason} \
        --seed {params.seed} \
        -ir {input.ref} \
        -iv {input.vcf} \
        -n {params.num_reads} \
        --illumina-read-length 150 \
        --num-threads {threads} \
        -o {output.r1} -or {output.r2} &> {log}
        """
        

rule bwa:
    input:
        ref = "sim/ref.fa",
        ref_idx = rules.bwa_index.output,
        r1= "reads/{sample}_R1.fq",
        r2= "reads/{sample}_R2.fq"
    params:
        rg = r"'@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA'"
    conda:
        "envs/env.yaml"
    log:
        "logs/bwa/{sample.txt}"
    output:
        bam=temp("bams/{sample}.bam"),
        bai=temp("bams/{sample}.bam.bai")
    shell:
        "bwa mem -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort - 2>> {log}| samtools view -b > {output.bam} 2>> {log} && samtools index {output.bam} 2>> {log}"

rule mosdepth:
    input:
        bam = "bams/{sample}.bam",
        index = "bams/{sample}.bam.bai"
    output:
        "depths/{sample}.per-base.d4"
    params:
        prefix=lambda wc, output: output[0].replace(".per-base.d4", ""),
    conda:
        "envs/mosdepth.yaml"
    shadow: "minimal"
    shell:
        """
        mosdepth --d4 {params.prefix} {input.bam}
        """
rule merge_d4:
    input:
        expand("depths/{sample}.per-base.d4", sample=SAMPLES)
    output:
        "d4/merged.d4"
    conda:
        "envs/mosdepth.yaml"
    shell:
        """
        d4tools merge {input} {output}
        """

rule moprs:
    input:
        moprs_binary = MOPRS,
        d4 = "d4/merged.d4"
    output:
        "callable_sites/moprs.bed"
    params:
        min_mean_depth = f"-u {MIN_DEPTH_MEAN}",
        max_mean_depth = f"-U {MAX_DEPTH_MEAN}",
        proportion = f"-d {1-MAX_MISSING}"
    threads: MOPRS_THREADS
    log: "logs/moprs.txt"
    benchmark:
        "benchmarks/moprs.txt"
    shell:
        """
        {input.moprs_binary} {params} -c -t {threads} --d4 {input.d4} > {output} 2> {log}
        """


rule bcftools:
    input:
        ref = "sim/ref.fa",
        bams = expand("bams/{sample}.bam", sample=SAMPLES)
    output:
        vcf="vcf/vars.vcf.gz",
        tbi="vcf/vars.vcf.gz.tbi"
    conda:
        "envs/env.yaml"
    shell:
        "(bcftools mpileup -f {input.ref} {input.bams} | bcftools call -m -f GQ | bgzip -c > {output.vcf}) && tabix -p vcf {output.vcf}"


rule create_invar_only:
    input:
        vcf="vcf/vars.vcf.gz",
        tbi="vcf/vars.vcf.gz.tbi"
    output:
        vcf="vcf/invars.vcf.gz",
        tbi="vcf/invars.vcf.gz.tbi"
    params:
        min_mean_depth = MIN_DEPTH_MEAN,
        max_mean_depth = MAX_DEPTH_MEAN,
        max_missing = MAX_MISSING,
    conda:
        "envs/env.yaml"
    shell:
        """
        (vcftools --gzvcf {input.vcf} \
        --max-maf 0 \
        --max-missing {params.max_missing} \
        --min-meanDP {params.min_mean_depth} \
        --max-meanDP {params.max_mean_depth} \
        --recode --stdout | bgzip -c > {output.vcf}) && tabix -p vcf {output.vcf}
        """

    
rule vcftools_filter:
    input:
        vcf="vcf/vars.vcf.gz",
        tbi="vcf/vars.vcf.gz.tbi"
    output:
        vcf="vcf/filtered_vars.vcf.gz",
        tbi="vcf/filtered_vars.vcf.gz.tbi"
    params:
        min_mean_depth = MIN_DEPTH_MEAN,
        max_mean_depth = MAX_DEPTH_MEAN,
        max_missing = MAX_MISSING,
    conda:
        "envs/env.yaml"
    shell:
        """
        (vcftools --gzvcf {input.vcf} \
        --remove-indels \
        --max-missing {params.max_missing} \
        --min-meanDP {params.min_mean_depth} \
        --max-meanDP {params.max_mean_depth} \
        --recode --stdout | bgzip -c > {output.vcf}) && tabix -p vcf {output.vcf}
        """

rule bcftools_concat:
    input:
        varvcf="vcf/filtered_vars.vcf.gz",
        vartbi="vcf/filtered_vars.vcf.gz.tbi",
        invarvcf="vcf/invars.vcf.gz",
        invartbi="vcf/invars.vcf.gz.tbi"
    output:
        vcf="vcf/allsites_filtered.vcf.gz",
        tbi="vcf/allsites_filtered.vcf.gz.tbi"
    conda:
        "envs/env.yaml"
    shell:
        """
        (bcftools concat \
        --allow-overlaps \
        {input.varvcf} {input.invarvcf} \
        | bgzip -c {output.vcf}) && tabix -p vcf {output.vcf}
        """



rule write_pops_file:
    output: "vcf/pops.txt"
    run:
        lines = "\n".join([f"{s}\tpop1" for s in SAMPLES])
        with open(output[0], "w") as f:
            f.writelines(lines)

rule pixy:
    input:
        vcf="vcf/allsites_filtered.vcf.gz",
        tbi="vcf/allsites_filtered.vcf.gz.tbi",
        pops = "vcf/pops.txt"
    params:
        windows=1000
    output:
        "pixy_pi.txt"
    conda: "envs/pixy.yaml"
    shadow: "minimal"
    shell:
        "pixy --stats pi --vcf {input.vcf} --window_size {params.windows} --populations {input.pops}"