import pandas as pd

NUM_REPS=range(1,3) #seed for msprime gotta be greater than 0
THREADS_PER_RUN=16
SNAKEFILE_PATH = Path(workflow.basedir,"Snakefile")
CONDA_PREFIX = Path(workflow.basedir,"")
CLAM_URL = "https://github.com/cademirch/clam.git"
MOSDEPTH_URL = "https://github.com/brentp/mosdepth/releases/download/v0.3.10/mosdepth_d4" # get latest mosdepth that builds d4 index
MASON_URL = "http://packages.seqan.de/mason2/mason2-2.0.9-Linux-x86_64.tar.xz"




rule all:
    input:
        pixy = expand("intermediate/{i}/pixy_pi.txt", i = NUM_REPS),
        clam = expand("intermediate/{i}/clam_pi.tsv", i = NUM_REPS),
        prime = expand("intermediate/{i}/msprime_pi_windows.tsv", i = NUM_REPS),
    output:
        pixy = "results/pixy.csv",
        clam = "results/clam.csv",
        prime = "results/msprime.csv",
    run:
        clam_dfs = [pd.read_csv(f, sep="\t") for f in input.clam]
        pixy_dfs = [pd.read_csv(f, sep="\t") for f in input.pixy]
        prime_dfs = [pd.read_csv(f, sep="\t") for f in input.prime]
        clam_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(clam_dfs)]
        pixy_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(pixy_dfs)]
        prime_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(prime_dfs)]

        clam = pd.concat(clam_dfs_with_index).to_csv(output.clam)
        pixy = pd.concat(pixy_dfs_with_index).to_csv(output.pixy)
        prime = pd.concat(prime_dfs_with_index).to_csv(output.prime)

rule setup_clam:
    output:
        "clam/target/release/clam"
    params:
        url = CLAM_URL
    shell:
        """
        rm -rf clam
        git clone {params.url}
        cd clam
        cargo build --release
        """

rule setup_mosdepth:
    output:
        "mosdepth_d4"
    params:
        url = MOSDEPTH_URL
    shell:
        """
        wget {params.url}
        chmod +x mosdepth_d4
        """

rule setup_mason:
    output:
        "mason2-2.0.9-Linux-x86_64/bin/mason_simulator"
    params:
        url = MASON_URL
    shell:
        """
        wget {params.url}
        tar -xf mason2-2.0.9-Linux-x86_64.tar.xz
        """
    
rule run:
    input:
        "clam/target/release/clam",
        "mosdepth_d4",
        "mason2-2.0.9-Linux-x86_64/bin/mason_simulator",
        snk = SNAKEFILE_PATH,
    output:
        pixy = temp("intermediate/{i}/pixy_pi.txt"),
        clam = temp("intermediate/{i}/clam_pi.tsv"),
        prime =temp("intermediate/{i}/msprime_pi_windows.tsv")
    params:
        conda_prefix = CONDA_PREFIX
    conda:
        "snakemake"
    threads: THREADS_PER_RUN
    log: "logs/meta/{i}.log"
    shadow: "minimal"
    shell:
        """
        snakemake -s {input.snk} -d intermediate/{wildcards.i} --config seed={wildcards.i} --use-conda --conda-prefix {params.conda_prefix} --cores {threads} --nolock &> {log}
        """
        
        