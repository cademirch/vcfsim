from pathlib import Path
import numpy as np
import random
import pandas as pd
# wildcard_constraints:
#     missingtype="[^\W_]"

CLAM = Path(workflow.basedir, "clam/target/release/clam")
SAMPLES = [f"tsk_{i}" for i in range(10)]
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
        vcf = "vcf/allsites.vcf",
        populations = "vcf/pops.txt"
    output:
        gtsvcf = "vcf/prop_{prop}/allsites_missing_gts.vcf",
        gtsbed1 = "callable/prop_{prop}/callable_gts.pop1.bedgraph",
        gtsbed2 = "callable/prop_{prop}/callable_gts.pop2.bedgraph",
        gtsgenome = "callable/prop_{prop}/genome.txt",
        sitesvcf = "vcf/prop_{prop}/allsites_missing_sites.vcf",
        sitesbed1 = "callable/prop_{prop}/callable_sites.pop1.bedgraph",
        sitesbed2 = "callable/prop_{prop}/callable_sites.pop2.bedgraph",
    threads: 2
    params:
        prop="{prop}"
    script:
        "scripts/add_missing.py"

rule d4tools_create:
    input:
        bed = "callable/prop_{prop}/callable_{missingtype}.{population}.bedgraph",
        genome = "callable/prop_{prop}/genome.txt" # for d4tools create
    output:
        "callable/prop_{prop}/callable_{missingtype}.{population}.d4"
    conda:
        "mosdepth"
    shell:
        """
        d4tools create -z -g {input.genome} {input.bed} {output}
        d4tools index build -s {output}
        """

rule d4tools_merge:
    input:
        expand("callable/prop_{{prop}}/callable_{{missingtype}}.{population}.d4", population=["pop1", "pop2"])
    output:
        "callable/prop_{prop}/callable_{missingtype}.d4"
    conda:
        "mosdepth"
    shell:
        """
        d4tools merge {input} {output}
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
        p = "vcf/pops.txt"
        
    output:
        allsites = "clam_allsites/clam_pi.tsv",
        dxy = "clam_allsites/clam_dxy.tsv"
        
    benchmark:
        "benchmarks/clam_allsites.txt"
    params:
        windows=WINDOW_SIZE
    shell:
        """
        {input.clam_binary} --quiet stat -p {input.p} -w {params.windows} -o clam_allsites {input.allsites}
        """

rule clam_missing:
    input:
        clam_binary = CLAM,
        varmissingvcf = "vcf/prop_{prop}/varsonly_{missingtype}.vcf.gz",

        varmissingvcf_idx = "vcf/prop_{prop}/varsonly_{missingtype}.vcf.gz.tbi",
        callable_sites = "callable/prop_{prop}/callable_{missingtype}.d4",
        p = "vcf/pops.txt",
    output:
        missing = "clam_missing/prop_{prop}/{missingtype}/clam_pi.tsv",
        dxy = "clam_missing/prop_{prop}/{missingtype}/clam_dxy.tsv"
    params:
        windows=WINDOW_SIZE,
        outdir= "clam_missing/prop_{prop}/{missingtype}"
    log:
        "logs/clam_{missingtype}_{prop}.txt"
    shell:
        """
        {input.clam_binary} --quiet stat -p {input.p} -w {params.windows} -o {params.outdir} {input.varmissingvcf} {input.callable_sites} || touch {output}
        """


def process_stats(pi_tsv, dxy_tsv, allsites_pi, allsites_dxy, tool):
    """
    Process pi and dxy statistics for a given tool and missing data files
    
    Args:
        pi_tsv: Path to pi statistics file
        dxy_tsv: Path to dxy statistics file
        allsites_pi: DataFrame with pi statistics for all sites
        allsites_dxy: DataFrame with dxy statistics for all sites
        tool: String indicating the tool used ('pixy' or 'clam')
    
    Returns:
        List of dictionaries with processed statistics
    """
    results = []
    missingtype = Path(pi_tsv).parent.name
    prop = float(Path(pi_tsv).parent.parent.name.split("_")[1])
    
    try:
        # Handle pi
        pi_df = pd.read_csv(pi_tsv, sep="\t")
        pi_dict = {
            "tool": tool,
            "stat": "pi",
            "allsites_value": allsites_pi['avg_pi' if tool == 'pixy' else 'pi'].iloc[0],
            "missing_value": pi_df['avg_pi' if tool == 'pixy' else 'pi'].iloc[0],
            "missingtype": missingtype,
            "prop_missing": prop
        }
        results.append(pi_dict)
        
        # Handle dxy
        dxy_df = pd.read_csv(dxy_tsv, sep="\t")
        dxy_dict = {
            "tool": tool,
            "stat": "dxy",
            "allsites_value": allsites_dxy['avg_dxy' if tool == 'pixy' else 'dxy'].iloc[0],
            "missing_value": dxy_df['avg_dxy' if tool == 'pixy' else 'dxy'].iloc[0],
            "missingtype": missingtype,
            "prop_missing": prop
        }
        results.append(dxy_dict)
        
    except pd.errors.EmptyDataError:
        pass
        
    return results


rule combine_clam:
    input:
        allsites_pi = "clam_allsites/clam_pi.tsv",
        allsites_dxy = "clam_allsites/clam_dxy.tsv",
        missing_pi = expand("clam_missing/prop_{prop}/{missingtype}/clam_pi.tsv", prop=PROPS, missingtype=MISSINGTYPES),
        missing_dxy = expand("clam_missing/prop_{prop}/{missingtype}/clam_dxy.tsv", prop=PROPS, missingtype=MISSINGTYPES)
    output:
        # raw = "clam_pi_dxy.tsv",
        ratio = "clam_ratios.tsv"
    run:
        import pandas as pd
        from pathlib import Path
        
        # Read allsites data first
        allsites_pi = pd.read_csv(input.allsites_pi, sep="\t")
        allsites_dxy = pd.read_csv(input.allsites_dxy, sep="\t")
        
        results = []
        for pi_tsv, dxy_tsv in zip(input.missing_pi, input.missing_dxy):
            results.extend(process_stats(pi_tsv, dxy_tsv, allsites_pi, allsites_dxy, "clam"))
        
        result_df = pd.DataFrame(results)
        result_df.to_csv(output.ratio, sep="\t", index=False)
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

rule pixy_allsites:
    input:
        allsites = "vcf/allsites.vcf.gz",
        
        pops = "vcf/pops.txt"
    output:
        allsites = "pixy_allsites/pixy_pi.txt",
        dxy = "pixy_allsites/pixy_dxy.txt",
        
    params:
        windows=WINDOW_SIZE
    benchmark:
        "benchmarks/pixy_allsites.txt"
    conda:
        "pixy"
    shell:
        """
        pixy --stats pi dxy fst --output_folder pixy_allsites --vcf {input.allsites} --window_size {params.windows} --populations {input.pops}
        """

rule pixy_missing:
    input:
        varmissingvcf = "vcf/prop_{prop}/allsites_missing_{missingtype}.vcf.gz",
        pops = "vcf/pops.txt"
    output:
        missing = "pixy_missing/prop_{prop}/{missingtype}/pixy_pi.txt",
        dxy = "pixy_missing/prop_{prop}/{missingtype}/pixy_dxy.txt"
    params:
        windows=WINDOW_SIZE,
        outdir= "pixy_missing/prop_{prop}/{missingtype}"
    
    conda:
        "pixy"
    shell:
        """
        pixy --stats pi dxy fst --output_folder {params.outdir} --vcf {input.varmissingvcf} --window_size {params.windows} --populations {input.pops} || touch {output}
        """


rule combine_pixy:
    input:
        allsites_pi = "pixy_allsites/pixy_pi.txt",
        allsites_dxy = "pixy_allsites/pixy_dxy.txt",
        missing_pi = expand("pixy_missing/prop_{prop}/{missingtype}/pixy_pi.txt", prop=PROPS, missingtype=MISSINGTYPES),
        missing_dxy = expand("pixy_missing/prop_{prop}/{missingtype}/pixy_dxy.txt", prop=PROPS, missingtype=MISSINGTYPES)
    output:
        ratio = "pixy_ratios.tsv"
    run:
        import pandas as pd
        from pathlib import Path
        
        # Read allsites data first
        allsites_pi = pd.read_csv(input.allsites_pi, sep="\t")
        allsites_dxy = pd.read_csv(input.allsites_dxy, sep="\t")
        
        results = []
        for pi_tsv, dxy_tsv in zip(input.missing_pi, input.missing_dxy):
            results.extend(process_stats(pi_tsv, dxy_tsv, allsites_pi, allsites_dxy, "pixy"))
        
        result_df = pd.DataFrame(results)
        result_df.to_csv(output.ratio, sep="\t", index=False)



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
        ratio = "vcftools_ratios.tsv"
    run:
        import pandas as pd
        from pathlib import Path
        
        allsites_df = pd.read_csv(input.allsites, sep="\t")
        print("Allsites shape:", allsites_df.shape)
        print("Allsites head:", allsites_df.head())
        
        results = []
        
        for tsv in input.missing:
            missingtype = Path(tsv).parent.name
            prop = float(Path(tsv).parent.parent.name.split("_")[1])
            try:
                df = pd.read_csv(tsv, sep="\t")
            
                
                if not df.empty and not allsites_df.empty:
                    result_dict = {
                        "tool": "vcftools",
                        "stat": "pi",
                        "allsites_value": allsites_df['PI'].iloc[0] if len(allsites_df) > 0 else None,
                        "missing_value": df['PI'].iloc[0] if len(df) > 0 else None,
                        "missingtype": missingtype,
                        "prop_missing": prop
                    }
                    results.append(result_dict)
            except pd.errors.EmptyDataError:
                print(f"Empty file: {tsv}")
                continue
            except Exception as e:
                print(f"Error processing {tsv}: {str(e)}")
                continue
                
        result_df = pd.DataFrame(results)
        result_df.to_csv(output.ratio, sep="\t", index=False)