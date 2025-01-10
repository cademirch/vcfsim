from pathlib import Path
import numpy as np
import random
import pandas as pd
# wildcard_constraints:
#     missingtype="[^\W_]"

CLAM = Path(workflow.basedir, "clam/target/release/clam")
SAMPLES = [f"sample_{i}" for i in range(10)]
SEQLEN = 10_000
SEED = int(config.get("seed", 1234))
logger.warning(f"{SEED=}")
NE=1e4
MU=2.5e-8
PROPS = list(np.round(np.arange(0, 1, 0.01), 2))
MISSINGTYPES = ["sites", "gts"]
WINDOW_SIZE=SEQLEN


rule all:
    input:
        "vcf/allsites.vcf",
        "pixy_allsites/pixy_pi.txt",
        "clam_allsites/clam_pi.tsv",
        "pixy_missing/prop_{prop}/{missingtype}/pixy_pi.txt"
        "clam_missing/prop_{prop}/{missingtype}/clam_pi.tsv"

rule write_pops_file:
    output: "vcf/pops.txt"
    run:
        lines = "\n".join([f"tsk_{s}\tpop1" for s in range(len(SAMPLES))])
        with open(output[0], "w") as f:
            f.writelines(lines)

rule simulate:
    output:
        "sim/sim_results.txt",
        "msprime_pi_windows.tsv",
        "sim/ref.fa",
        "sim/all.vcf.gz",
        "sim/all.vcf.gz.tbi"
    conda:
        "clam_benchmarking"
    threads: 1
    params:
        ne = NE,
        mu = MU,
        sample_size = len(SAMPLES),
        seed = SEED,
        outdir = "sim",
        seqlen = SEQLEN,
        windows=WINDOW_SIZE,
        indv_vcfs=False
    script:
        "scripts/simulate.py"

rule samtools_faidx:
    input:
        "sim/ref.fa",
    output:
        "sim/ref.fa.fai",
    shell:
        "samtools faidx {input}"
rule add_invar:
    input:
        vcf = "sim/all.vcf.gz",
        vcfidx = "sim/all.vcf.gz.tbi",
        ref = "sim/ref.fa",
        refidx = "sim/ref.fa.fai",
    output:
        vcf = "vcf/allsites.vcf",
        vcfgz = "vcf/allsites.vcf.gz",
        vcfidx = "vcf/allsites.vcf.gz.tbi"
        
    conda:
        "clam_benchmarking"
    script:
        "scripts/addinvar.py"
rule add_missing:
    input:
        vcf = "vcf/allsites.vcf"
    output:
        gtsvcf = "vcf/prop_{prop}/allsites_missing_gts.vcf",
        gtsbed = "callable/prop_{prop}/callable_gts.bedgraph",
        gtsgenome = "callable/prop_{prop}/genome.txt",
        sitesvcf = "vcf/prop_{prop}/allsites_missing_sites.vcf",
        sitesbed = "callable/prop_{prop}/callable_sites.bedgraph",
    threads: 2
    params:
        prop="{prop}"
    script:
        "scripts/add_missing.py"


rule variants_only_missing:
    input:
        allsitesvcf = "vcf/prop_{prop}/allsites_missing_{missingtype}.vcf",
    output:
        allsitesvcf = "vcf/prop_{prop}/allsites_missing_{missingtype}.vcf.gz",
        allsitesvcfidx = "vcf/prop_{prop}/allsites_missing_{missingtype}.vcf.gz.tbi",
        varsonly = "vcf/prop_{prop}/varsonly_{missingtype}.vcf",
        varsonlygz = "vcf/prop_{prop}/varsonly_{missingtype}.vcf.gz",
        varsonlyidx = "vcf/prop_{prop}/varsonly_{missingtype}.vcf.gz.tbi",
        
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} /^#/ {{print $0}} !/^#/ && $5 != "." {{print $0}}' {input.allsitesvcf} > {output.varsonly}
        bgzip -c {input} > {output.allsitesvcf}
        tabix -p vcf {output.allsitesvcf}
        bgzip -c {output.varsonly} > {output.varsonlygz}
        tabix -p vcf {output.varsonlygz}
        """

rule d4tools_create:
    input:
        bed = "callable/prop_{prop}/callable_{missingtype}.bedgraph",
        genome = "callable/prop_{prop}/genome.txt" # for d4tools create
    output:
        "callable/prop_{prop}/callable_{missingtype}.d4"
    conda:
        "mosdepth"
    shell:
        """
        d4tools create -z -g {input.genome} {input.bed} {output}
        d4tools index build -s {output}
        """

    
rule clam_allsites:
    input:
        clam_binary = CLAM,
        allsites = "vcf/allsites.vcf.gz",
        allsitesidx = "vcf/allsites.vcf.gz.tbi",
        
    output:
        allsites = "clam_allsites/clam_pi.tsv",
        
    benchmark:
        "benchmarks/clam_allsites.txt"
    params:
        windows=WINDOW_SIZE
    shell:
        """
        {input.clam_binary} --quiet stat -w {params.windows} -o clam_allsites {input.allsites}
        """
rule clam_missing:
    input:
        clam_binary = CLAM,
        varmissingvcf = "vcf/prop_{prop}/varsonly_{missingtype}.vcf.gz",

        varmissingvcf_idx = "vcf/prop_{prop}/varsonly_{missingtype}.vcf.gz.tbi",
        callable_sites = "callable/prop_{prop}/callable_{missingtype}.d4",
    output:
        missing = "clam_missing/prop_{prop}/{missingtype}/clam_pi.tsv"
    params:
        windows=WINDOW_SIZE,
        outdir= "clam_missing/prop_{prop}/{missingtype}"
    log:
        "logs/clam_{missingtype}_{prop}.txt"
    shell:
        """
        {input.clam_binary} --quiet stat -w {params.windows} -o {params.outdir} {input.varmissingvcf} {input.callable_sites} || touch {output}
        """

rule pixy_allsites:
    input:
        allsites = "vcf/allsites.vcf.gz",
        
        pops = "vcf/pops.txt"
    output:
        allsites = "pixy_allsites/pixy_pi.txt",
        
    params:
        windows=WINDOW_SIZE
    benchmark:
        "benchmarks/pixy_allsites.txt"
    conda:
        "pixy"
    shell:
        """
        pixy --stats pi --output_folder pixy_allsites --vcf {input.allsites} --window_size {params.windows} --populations {input.pops}
        """

rule pixy_missing:
    input:
        varmissingvcf = "vcf/prop_{prop}/allsites_missing_{missingtype}.vcf.gz",
        pops = "vcf/pops.txt"
    output:
        missing = "pixy_missing/prop_{prop}/{missingtype}/pixy_pi.txt"
    params:
        windows=WINDOW_SIZE,
        outdir= "pixy_missing/prop_{prop}/{missingtype}"
    
    conda:
        "pixy"
    shell:
        """
        pixy --stats pi --output_folder {params.outdir} --vcf {input.varmissingvcf} --window_size {params.windows} --populations {input.pops} || touch {output}
        """