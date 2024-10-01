MASON_BIN = "/Users/cade/dev/vcf-stuff/sim/mason2-2.0.9-Mac-x86_64/bin/mason_simulator"
SAMPLES = [f"sample_{i}" for i in range(20)]
SEQLEN = 10000
COVERAGE = 60
NUM_READS = (SEQLEN * COVERAGE) // 300
SEED = 1234
NE=1e6
MU=1e-8

rule all:
    input:
        expand("bams/{sample}.bam", sample=SAMPLES),
        "vars.vcf.gz",
        "pops.txt"

rule simulate:
    output:
        expand("sim/{sample}.fa", sample=SAMPLES),
        "sim/sim.vcf",
        "sim/sim_results.txt",
        "sim/pi_windows.txt",
        "sim/ref.fa"
    params:
        ne = NE,
        mu = MU,
        sample_size = len(SAMPLES),
        seed = SEED,
        outdir = "sim",
        seqlen = SEQLEN,
        windows=1000
    script:
        "msprime_test.py"
rule bwa_index:
    input:
        "sim/ref.fa"
    output:
        "sim/ref.fa.bwt"
    shell:
        "bwa index {input}"

rule simulate_reads:
    input:
        mason = MASON_BIN,
        ref = "sim/{sample}.fa"
    params:
        num_reads = NUM_READS,
        seed = SEED
    output:
        r1= "reads/{sample}_R1.fq",
        r2= "reads/{sample}_R2.fq"
    shell:
        "{input.mason} --seed {params.seed} -ir {input.ref} -n {params.num_reads} --illumina-read-length 150 -o {output.r1} -or {output.r2}"


rule bwa:
    input:
        ref = "sim/ref.fa",
        ref_idx = rules.bwa_index.output,
        r1= "reads/{sample}_R1.fq",
        r2= "reads/{sample}_R2.fq"
    params:
        rg = r"'@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA'"
    output:
        "bams/{sample}.bam"
    shell:
        "bwa mem -R {params.rg} {input.ref} {input.r1} {input.r2} | samtools sort - | samtools view -b > {output}"

rule bcftools:
    input:
        ref = "sim/ref.fa",
        bams = expand("bams/{sample}.bam", sample=SAMPLES)
    output:
        "vars.vcf.gz"
    shell:
        "bcftools mpileup -f {input.ref} {input.bams} | bcftools call -m -Oz -f GQ -o {output} && tabix -p vcf {output}"

rule write_pops_file:
    output: "pops.txt"
    run:
        lines = "\n".join([f"{s}\tpop1" for s in SAMPLES])
        with open(output[0], "w") as f:
            f.writelines(lines)

rule pixy:
    input:
        vcf = "vars.vcf.gz",
        pops = "pops.txt"
    params:
        windows=1000
    output:
        "pixy_pi.txt"
    shell:
        "pixy --stats pi --vcf {input.vcf} --window_size {params.windows} --populations {input.pops}"