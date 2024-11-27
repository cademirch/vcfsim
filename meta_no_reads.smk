import pandas as pd

NUM_REPS=range(1,3) #seed for msprime gotta be greater than 0
THREADS_PER_RUN=16
SNAKEFILE_PATH = Path(workflow.basedir,"sim_no_reads.smk")
CONDA_PREFIX = Path(workflow.basedir,"")
CLAM_URL = "https://github.com/cademirch/clam.git"



rule all:
    input:
        pixy = expand("intermediate_missing/{i}/pixy_pi.txt", i = NUM_REPS),
        clam = expand("intermediate_missing/{i}/clam_pi.tsv", i = NUM_REPS),
        prime = expand("intermediate_missing/{i}/msprime_pi_windows.tsv", i = NUM_REPS),
    output:
        pixy = "results/pixy.csv",
        clam = "results/clam.csv",
        prime = "results/msprime.csv",
    run:
        def calculate_proportion(idx):
            if idx <= 10000:
                return (idx // 10) / 1000.0, "gts"  # Proportion of missing gts
            else:
                return ((idx - 10000) // 10) / 1000.0, "sites"  # Proportion of missing sites

        clam_dfs = [pd.read_csv(f, sep="\t") for f in input.clam]
        pixy_dfs = [pd.read_csv(f, sep="\t") for f in input.pixy]
        prime_dfs = [pd.read_csv(f, sep="\t") for f in input.prime]

        # Add source index and proportion columns
        clam_dfs_with_index = [
            df.assign(source_idx=i, proportion=calculate_proportion(i)[0], type=calculate_proportion(i)[1])
            for i, df in enumerate(clam_dfs)
        ]
        pixy_dfs_with_index = [
            df.assign(source_idx=i, proportion=calculate_proportion(i)[0], type=calculate_proportion(i)[1])
            for i, df in enumerate(pixy_dfs)
        ]
        prime_dfs_with_index = [
            df.assign(source_idx=i, proportion=calculate_proportion(i)[0], type=calculate_proportion(i)[1])
            for i, df in enumerate(prime_dfs)
        ]

        # Concatenate and save to CSV
        pd.concat(clam_dfs_with_index).to_csv(output.clam, index=False)
        pd.concat(pixy_dfs_with_index).to_csv(output.pixy, index=False)
        pd.concat(prime_dfs_with_index).to_csv(output.prime, index=False)

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
def set_missing_gts(wildcards):
    if int(wildcards.i) <= 10000:
        return (int(wildcards.i) // 10) / 10000.0
    return 0

def set_missing_sites(wildcards):
    if int(wildcards.i) > 10000:
        return ((int(wildcards.i) - 10_000) // 10) / 10000.0
    return 0


rule run:
    input:
        "clam/target/release/clam",
        snk = SNAKEFILE_PATH,
    output:
        pixy = temp("intermediate_missing/{i}/pixy_pi.txt"),
        clam = temp("intermediate_missing/{i}/clam_pi.tsv"),
        prime =temp("intermediate_missing/{i}/msprime_pi_windows.tsv")
    params:
        conda_prefix = CONDA_PREFIX
        missing_gts = set_missing_gts,
        missing_sites = set_missing_sites
    conda:
        "snakemake"
    threads: THREADS_PER_RUN
    log: "logs/meta/{i}.log"
    shadow: "minimal"
    shell:
        """
        snakemake -s {input.snk} -d intermediate_missing/{wildcards.i} --config seed={wildcards.i} missing_sites={params.missing_sites} missing_gts={params.missing_gts} --use-conda --conda-prefix {params.conda_prefix} --cores {threads} --nolock &> {log}
        """
        
        