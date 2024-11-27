from pathlib import Path
import random


CLAM = Path(workflow.basedir, "clam/target/release/clam")
SAMPLES = [f"sample_{i}" for i in range(10)]
SEQLEN = 10_000
SEED = int(config.get("seed", 1234))
logger.warning(f"{SEED=}")
NE=1e6
MU=1e-8
PROP_MISSING_SITES = float(config.get("missing_sites", 0.0))
PROP_MISSING_GENOTYPES = float(config.get("missing_gts", 0.0))

rule all:
    input:"clam_pi.tsv", "pixy_pi.tsv"

rule simulate:
    output:
        temp(expand("sim/{sample}.vcf", sample=SAMPLES)),
        "sim/sim_results.txt",
        "msprime_pi_windows.tsv",
        "sim/ref.fa",
        "sim/all.vcf"
    conda:
        "envs/env.yaml"
    threads: 4
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

rule add_invar:
    input:
        vcf = "sim/all.vcf",
        ref = "sim/ref.fa",
    output:
        vcf = "vcf/allsites.vcf"
    conda:
        "envs/env.yaml"
    script:
        "scripts/addinvar.py"

rule add_missing:
    input:
        vcf = "vcf/allsites.vcf"
    output:
        vcf = "vcf/allsites_missing.vcf",
        bed = "callable/callable.bedgraph",
        genome = "callable/genome.txt" 
    params:
        prop_missing_genotypes = PROP_MISSING_GENOTYPES,
        prop_missing_sites = PROP_MISSING_SITES
    script:
        "scripts/add_missing.py"

rule d4tools_create:
    input:
        bed = "callable/callable.bedgraph",
        genome = "callable/genome.txt" # for d4tools create
    output:
        "callable/callable_sites.d4"
    shell:
        """
        d4tools create -z -g {input.genome} {input.bed} {output}
        """

rule variants_only:
    input:
        allsitesvcf = "vcf/allsites.vcf",
        missingvcf = "vcf/allsites_missing.vcf",
    output:
        allsitesvcf = "vcf/varsonly_allsites.vcf",
        missingvcf = "vcf/varsonly_allsites_missing.vcf",
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} /^#/ {{print $0}} !/^#/ && $5 != "." {{print $0}}' {input.allsitesvcf} > {output.allsitesvcf}
        awk 'BEGIN {{OFS="\\t"}} /^#/ {{print $0}} !/^#/ && $5 != "." {{print $0}}' {input.missingvcf} > {output.missingvcf}
        """

rule bgzip_all:
    input:
        allsitesvcf = "vcf/allsites.vcf",
        missingvcf = "vcf/allsites_missing.vcf",
        varallsitesvcf = "vcf/varsonly_allsites.vcf",
        varmissingvcf = "vcf/varsonly_allsites_missing.vcf",
    output:
        vcfs = ["vcf/allsites.vcf.gz","vcf/allsites_missing.vcf.gz","vcf/varsonly_allsites.vcf.gz","vcf/varsonly_allsites_missing.vcf.gz",],
        idxs = [ "vcf/allsites.vcf.gz.tbi","vcf/allsites_missing.vcf.gz.tbi","vcf/varsonly_allsites.vcf.gz.tbi","vcf/varsonly_allsites_missing.vcf.gz.tbi",]
    run:
        for i,o in zip(input,output.vcfs):
            shell(f"bgzip -c {i} > {o}")
            shell(f"tabix -p vcf {o}")
        

rule clam:
    input:
        clam_binary = CLAM,
        allsites = "vcf/allsites.vcf.gz",
        varmissingvcf = "vcf/varsonly_allsites_missing.vcf.gz",
        callable_sites = "callable/callable_sites.d4",
    output:
        allsites = "clam_allsites/clam_pi.tsv",
        missing = "clam_missing/clam_pi.tsv"
    shell:
        """
        {input.clam_binary} stat -w 1000 -o clam_allsites {input.allsites}
        {input.clam_binary} stat -w 1000 -o clam_missing {input.varmissingvcf} {input.callable_sites}
        """

rule combine_clam:
    input:
        allsites = "clam_allsites/clam_pi.tsv",
        missing = "clam_missing/clam_pi.tsv"
    output:
        "clam_pi.tsv"
    run:
        import pandas as pd

        # Read both input TSV files
        allsites_df = pd.read_csv(input.allsites, sep="\t")
        missing_df = pd.read_csv(input.missing, sep="\t")

        # Add the 'vcf_type' column
        allsites_df['vcf_type'] = 'allsites'
        missing_df['vcf_type'] = 'missing'

        # Concatenate the two DataFrames
        combined_df = pd.concat([allsites_df, missing_df], ignore_index=True)

        # Save the combined DataFrame to the output file
        combined_df.to_csv(output[0], sep="\t", index=False)
rule write_pops_file:
    output: "vcf/pops.txt"
    run:
        lines = "\n".join([f"tsk_{s}\tpop1" for s in range(len(SAMPLES))])
        with open(output[0], "w") as f:
            f.writelines(lines)

rule pixy:
    input:
        allsites = "vcf/allsites.vcf.gz",
        varmissingvcf = "vcf/allsites_missing.vcf.gz",
        pops = "vcf/pops.txt"
    output:
        allsites = "pixy_allsites/pixy_pi.txt",
        missing = "pixy_missing/pixy_pi.txt"
    conda:
        "envs/pixy.yaml"
    shell:
        """
        pixy --stats pi --output_folder pixy_allsites --vcf {input.allsites} --window_size 1000 --populations {input.pops}
        pixy --stats pi --output_folder pixy_missing --vcf {input.varmissingvcf} --window_size 1000 --populations {input.pops}
        """


rule combine_pixy:
    input:
        allsites = "pixy_allsites/pixy_pi.txt",
        missing = "pixy_missing/pixy_pi.txt"
    output:
        "pixy_pi.tsv"
    run:
        import pandas as pd

        # Read both input TSV files
        allsites_df = pd.read_csv(input.allsites, sep="\t")
        missing_df = pd.read_csv(input.missing, sep="\t")

        # Add the 'vcf_type' column
        allsites_df['vcf_type'] = 'allsites'
        missing_df['vcf_type'] = 'missing'

        # Concatenate the two DataFrames
        combined_df = pd.concat([allsites_df, missing_df], ignore_index=True)

        # Save the combined DataFrame to the output file
        combined_df.to_csv(output[0], sep="\t", index=False)
