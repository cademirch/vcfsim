"""
Runs sim_no_reads.smk. This simulation just simulates vcf from msprime, adds invariant sites, then randomly removes sites/genotypes for missing data
For clam, since we dont have bams in this simulation, we just record the number of callable individuals when we remove sites/genotypes.
"""

localrules: setup_clam, setup_mason, setup_mosdepth, all
import pandas as pd

NUM_REPS=range(1,1001) #seed for msprime gotta be greater than 0
THREADS_PER_RUN=32
SNAKEFILE_PATH = Path(workflow.basedir,"sim_no_reads.smk")
CONDA_PREFIX = Path(workflow.basedir,"")
CLAM_URL = "https://github.com/cademirch/clam.git"
MASON_URL = "http://packages.seqan.de/mason2/mason2-2.0.9-Linux-x86_64.tar.xz"
CLAM_PATH = Path(workflow.basedir,"clam/target/release/clam")
MASON_PATH = Path(workflow.basedir,"mason2-2.0.9-Linux-x86_64/bin/mason_simulator")



rule all:
    input:
        expand("results/no_reads/{i}/all.tsv",i=NUM_REPS),
    output:
        "results/no_reads/complete.tsv"
    run:
        dfs = []
        for i, f in enumerate(input):
            d = pd.read_csv(str(f), sep="\t")
            d["replicate"] = i
            dfs.append(d)
        df = pd.concat(dfs)
        df.to_csv(str(output), index=False, sep="\t")
        



rule setup_clam:
    output:
        CLAM_PATH
    params:
        url = CLAM_URL
    shell:
        """
        rm -rf clam
        git clone {params.url}
        cd clam
        cargo build --release
        """


rule setup_mason:
    output:
        MASON_PATH
    params:
        url = MASON_URL
    shell:
        """
        wget {params.url}
        tar -xf mason2-2.0.9-Linux-x86_64.tar.xz
        """
    
rule run:
    input:
        CLAM_PATH,
        MASON_PATH,
        snk = SNAKEFILE_PATH,
    output:
        temp("results/no_reads/{i}/all.tsv")
    params:
        conda_prefix = CONDA_PREFIX
    conda:
        "snakemake"
    threads: THREADS_PER_RUN
    log: "logs/meta/{i}.log"
    shadow: "minimal"
    shell:
        """
        snakemake -s {input.snk} -d results/no_reads/{wildcards.i} --config seed={wildcards.i} --use-conda --conda-prefix {params.conda_prefix} --cores {threads} --nolock --show-failed-logs --retries 3 &> {log}
        """
        
        