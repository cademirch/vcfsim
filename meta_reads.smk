"""
Runs sim_with_reads.smk. For this simulation, we simulate vcf with msprime. Then simulate reads for each sample in the vcf using mason. 
Then align with bwa and call variants. Get depth with mosdepth. Then do clam loci/stat and pixy.
"""
localrules: setup_clam, setup_mason, setup_mosdepth, all
import pandas as pd

NUM_REPS=range(1,101) #seed for msprime gotta be greater than 0
THREADS_PER_RUN=32
SNAKEFILE_PATH = Path(workflow.basedir,"sim_with_reads.smk")
CONDA_PREFIX = Path(workflow.basedir,"")
CLAM_URL = "https://github.com/cademirch/clam.git"
MASON_URL = "http://packages.seqan.de/mason2/mason2-2.0.9-Linux-x86_64.tar.xz"
CLAM_PATH = Path(workflow.basedir,"clam/target/release/clam")
MASON_PATH = Path(workflow.basedir,"mason2-2.0.9-Linux-x86_64/bin/mason_simulator")



rule all:
    input:
        pixy = expand("intermediate/{i}/pixy_pi.txt", i = NUM_REPS),
        clam = expand("intermediate/{i}/clam_pi.tsv", i = NUM_REPS),
        prime = expand("intermediate/{i}/msprime_pi_windows.tsv", i = NUM_REPS),
        vcftools = expand("intermediate/{i}/vcftools_pi.txt", i = NUM_REPS),
        pixy_benchmark = expand("intermediate/{i}/benchmarks/pixy/benchmark.txt", i = NUM_REPS),
        clam_benchmark = expand("intermediate/{i}/benchmarks/clam_stat/benchmark.txt", i = NUM_REPS)
    output:
        pixy = "results_with_reads/pixy.csv",
        clam = "results_with_reads/clam.csv",
        prime = "results_with_reads/msprime.csv",
        vcftools = "results_with_reads/vcftools.csv",
        bench = "results_with_reads/bench.csv"
    run:
        clam_dfs = [pd.read_csv(f, sep="\t") for f in input.clam]
        pixy_dfs = [pd.read_csv(f, sep="\t") for f in input.pixy]
        prime_dfs = [pd.read_csv(f, sep="\t") for f in input.prime]
        vcftools_dfs = [pd.read_csv(f, sep="\t") for f in input.vcftools]
        clam_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(clam_dfs)]
        pixy_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(pixy_dfs)]
        prime_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(prime_dfs)]
        vcftools_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(vcftools_dfs)]

        pd.concat(clam_dfs_with_index).to_csv(output.clam, index=False)
        pd.concat(pixy_dfs_with_index).to_csv(output.pixy, index=False)
        pd.concat(prime_dfs_with_index).to_csv(output.prime, index=False)
        pd.concat(vcftools_dfs_with_index).to_csv(output.vcftools, index=False)

        clam_benches = [pd.read_csv(f, sep="\t") for f in input.clam_benchmark]
        pixy_benches = [pd.read_csv(f, sep="\t") for f in input.pixy_benchmark]

        clam_bench = pd.concat(clam_benches)
        pixy_bench = pd.concat(pixy_benches)

        clam_bench["tool"] = "clam"
        pixy_bench["tool"] = "pixy"

        pd.concat([clam_bench, pixy_bench]).to_csv(output.bench, index=False)




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
        pixy = temp("intermediate/{i}/pixy_pi.txt"),
        clam = temp("intermediate/{i}/clam_pi.tsv"),
        vcftools = temp("intermediate/{i}/vcftools_pi.txt"),
        prime =temp("intermediate/{i}/msprime_pi_windows.tsv"),
        pixy_benchmark = temp("intermediate/{i}/benchmarks/pixy/benchmark.txt"),
        clam_benchmark = temp("intermediate/{i}/benchmarks/clam_stat/benchmark.txt")
    params:
        conda_prefix = CONDA_PREFIX
    conda:
        "snakemake"
    threads: THREADS_PER_RUN
    log: "logs/meta/{i}.log"
    shadow: "minimal"
    shell:
        """
        snakemake -s {input.snk} -d intermediate/{wildcards.i} --config seed={wildcards.i} --use-conda --conda-prefix {params.conda_prefix} --cores {threads} --nolock --show-failed-logs --retries 3 &> {log}
        """
        
        