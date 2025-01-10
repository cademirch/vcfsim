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
        clam = "clam_ratios.tsv",
        pixy = "pixy_ratios.tsv",
        vcftools = "vcftools_ratios.tsv",
    output:
        "all.tsv"
    run:
        dfs = []
        for f in input:
            d = pd.read_csv(str(f), sep="\t")
            dfs.append(d)
        df = pd.concat(dfs)
        df.to_csv(str(output), sep="\t", index=False)

        

rule simulate:
    output:
        "sim/sim_results.txt",
        "msprime_pi_windows.tsv",
        "sim/ref.fa",
        "sim/all.vcf.gz",
        "sim/all.vcf.gz.tbi"
    conda:
        "clam_benchmarking"
    threads: 2
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
        vcf = "vcf/allsites.vcf"
    params:
        chunk_size = 1_000_000
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

rule variants_only_allsites:
    input:
        allsitesvcf = "vcf/allsites.vcf",
        
    output:
        allsitesvcf = "vcf/varsonly_allsites.vcf",
        
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} /^#/ {{print $0}} !/^#/ && $5 != "." {{print $0}}' {input.allsitesvcf} > {output.allsitesvcf}
        """

rule variants_only_missing:
    input:
        missingvcf = "vcf/prop_{prop}/allsites_missing_{missingtype}.vcf",
    output:
        missingvcf = "vcf/prop_{prop}/varsonly_{missingtype}.vcf",
        
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} /^#/ {{print $0}} !/^#/ && $5 != "." {{print $0}}' {input.missingvcf} > {output.missingvcf}
        """

rule bgzip_allsites:
    input:
        allsitesvcf = "vcf/allsites.vcf",
        varallsitesvcf = "vcf/varsonly_allsites.vcf",
    output:
        vcfs = ["vcf/allsites.vcf.gz","vcf/varsonly_allsites_missing.vcf.gz",],
        idxs = [ "vcf/allsites.vcf.gz.tbi","vcf/varsonly_allsites_missing.vcf.gz.tbi",]
    run:
        for i,o in zip(input,output.vcfs):
            shell(f"bgzip -c {i} > {o}")
            shell(f"tabix -p vcf {o}")

rule bgzip_missing:
    input:
        missingvcf = "vcf/prop_{prop}/allsites_missing_{missingtype}.vcf",
        varmissingvcf = "vcf/prop_{prop}/varsonly_{missingtype}.vcf",
    output:
        vcfs = ["vcf/prop_{prop}/allsites_missing_{missingtype}.vcf.gz","vcf/prop_{prop}/varsonly_{missingtype}.vcf.gz"],
        idxs = ["vcf/prop_{prop}/allsites_missing_{missingtype}.vcf.gz.tbi","vcf/prop_{prop}/varsonly_{missingtype}.vcf.gz.tbi"]
    run:
        for i,o in zip(input,output.vcfs):
            shell(f"bgzip -c {i} > {o}")
            shell(f"tabix -p vcf {o}")
        

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


rule combine_clam:
    input:
        allsites = "clam_allsites/clam_pi.tsv",
        missing = expand("clam_missing/prop_{prop}/{missingtype}/clam_pi.tsv", prop=PROPS, missingtype=MISSINGTYPES)
    output:
        raw = "clam_pi.tsv",
        ratio = "clam_ratios.tsv"
    run:
        import pandas as pd
        from pathlib import Path
        missing_dfs = []
        for tsv in input.missing:
            missingtype = Path(tsv).parent.name
            prop = float(Path(tsv).parent.parent.name.split("_")[1])
            try:
                df = pd.read_csv(tsv, sep="\t")
                df["missingtype"] = missingtype
                df["prop"] = prop
                
                missing_dfs.append(df)
            except pd.errors.EmptyDataError:
                continue
        missing_df = pd.concat(missing_dfs)

        allsites_df = pd.read_csv(input.allsites, sep="\t")
        
        # # Concatenate the two DataFrames
        combined_df = pd.concat([allsites_df, missing_df], ignore_index=True)
        combined_df.to_csv(output[0], sep="\t", index=False)
        
        merged_df = allsites_df.merge(
            missing_df,
            on=["chrom", "start", "end", "population_name"],
            suffixes=("_allsites", "_missing")
        )

        # # Divide the 'pi' values: allsites / missing
        merged_df['pi_ratio'] =  merged_df['pi_missing'] / merged_df['pi_allsites']
        merged_df["tool"] = "clam"
        merged_df.to_csv(output[1], sep="\t", index=False)
        
rule write_pops_file:
    output: "vcf/pops.txt"
    run:
        lines = "\n".join([f"tsk_{s}\tpop1" for s in range(len(SAMPLES))])
        with open(output[0], "w") as f:
            f.writelines(lines)

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


rule combine_pixy:
    input:
        allsites = "pixy_allsites/pixy_pi.txt",
        missing = expand("pixy_missing/prop_{prop}/{missingtype}/pixy_pi.txt", prop=PROPS, missingtype=MISSINGTYPES)
    output:
        raw = "pixy_pi.tsv",
        ratio = "pixy_ratios.tsv"
    run:
        import pandas as pd
        from pathlib import Path
        missing_dfs = []
        for tsv in input.missing:
            missingtype = Path(tsv).parent.name
            prop = float(Path(tsv).parent.parent.name.split("_")[1])
            try:
                df = pd.read_csv(tsv, sep="\t")
                df["missingtype"] = missingtype
                df["prop"] = prop
                
                missing_dfs.append(df)
            except pd.errors.EmptyDataError:
                continue
        missing_df = pd.concat(missing_dfs)

        allsites_df = pd.read_csv(input.allsites, sep="\t")
        
        # # Concatenate the two DataFrames
        combined_df = pd.concat([allsites_df, missing_df], ignore_index=True)
        combined_df.to_csv(output[0], sep="\t", index=False)
        
        merged_df = allsites_df.merge(
            missing_df,
            on=["chromosome", "window_pos_1", "window_pos_2", "pop"],
            suffixes=("_allsites", "_missing")
        )

        # # Divide the 'pi' values: allsites / missing
        merged_df['pi_ratio'] =  merged_df['avg_pi_missing'] / merged_df['avg_pi_allsites']
        merged_df["tool"] = "pixy"
        merged_df.to_csv(output[1], sep="\t", index=False)



rule vcftools_allsites:
    input:
        vcf = "vcf/allsites.vcf.gz",
      
    output:
        "vcftools_allsites/vcftools_pi.txt",
    params:
        windows=WINDOW_SIZE
    conda:
        "clam_benchmarking"
    shell:
        """
        vcftools --gzvcf {input.vcf} --window-pi {params.windows} --stdout > {output}
        """

rule vcftools_missing:
    input:
        vcf = "vcf/prop_{prop}/allsites_missing_{missingtype}.vcf.gz",
      
    output:
       "vcftools_missing/prop_{prop}/{missingtype}/vcftools_pi.txt"
    params:
        windows=WINDOW_SIZE,
    
    conda:
        "clam_benchmarking"
    shell:
        """
        vcftools --gzvcf {input.vcf} --window-pi {params.windows} --stdout > {output}
        """


rule combine_vcftools:
    input:
        allsites = "vcftools_allsites/vcftools_pi.txt",
        missing = expand("vcftools_missing/prop_{prop}/{missingtype}/vcftools_pi.txt", prop=PROPS, missingtype=MISSINGTYPES)
    output:
        raw = "vcftools_pi.tsv",
        ratio = "vcftools_ratios.tsv"
    run:
        import pandas as pd
        from pathlib import Path
        missing_dfs = []
        for tsv in input.missing:
            missingtype = Path(tsv).parent.name
            prop = float(Path(tsv).parent.parent.name.split("_")[1])
            try:
                df = pd.read_csv(tsv, sep="\t")
                df["missingtype"] = missingtype
                df["prop"] = prop
                
                missing_dfs.append(df)
            except pd.errors.EmptyDataError:
                continue
        missing_df = pd.concat(missing_dfs)

        allsites_df = pd.read_csv(input.allsites, sep="\t")
        
        # # Concatenate the two DataFrames
        combined_df = pd.concat([allsites_df, missing_df], ignore_index=True)
        combined_df.to_csv(output[0], sep="\t", index=False)
        
        merged_df = allsites_df.merge(
            missing_df,
            on=["CHROM", "BIN_START", "BIN_END"],
            suffixes=("_allsites", "_missing")
        )

        # # Divide the 'pi' values: allsites / missing
        merged_df['pi_ratio'] =  merged_df['PI_missing'] / merged_df['PI_allsites']
        merged_df["tool"] = "vcftools"
        merged_df.to_csv(output[1], sep="\t", index=False)