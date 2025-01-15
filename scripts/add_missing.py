import random
import argparse
import gzip
from pathlib import Path
from multiprocessing import Process
import numpy as np
import pandas as pd


def parse_populations(pop_file):
    """
    Parse the population assignment file.

    Args:
        pop_file (str): Path to the tab-separated file containing sample_name and population

    Returns:
        dict: Mapping of sample names to populations
        set: Set of unique populations
    """
    pop_df = pd.read_csv(
        pop_file, sep="\t", header=None, names=["sample_name", "population"]
    )
    pop_map = dict(zip(pop_df["sample_name"], pop_df["population"]))
    unique_pops = set(pop_df["population"])
    return pop_map, unique_pops


def parse_vcf(input_vcf):
    with (
        gzip.open(input_vcf, "rt")
        if input_vcf.endswith(".gz")
        else open(input_vcf, "r")
    ) as vcf:
        header = []
        body = []
        samples = None
        for line in vcf:
            if line.startswith("#"):
                header.append(line.strip())
                if line.startswith("#CHROM"):
                    # Extract sample names from the header
                    samples = line.strip().split("\t")[9:]
            else:
                body.append(line.strip().split("\t"))
        return header, body, samples


def write_vcf(header, body, output_vcf):
    with open(output_vcf, "w") as vcf:
        vcf.write("\n".join(header) + "\n")
        for row in body:
            vcf.write("\t".join(row) + "\n")


def create_bedgraph(body, samples, pop_map, bedgraph_base):
    """
    Create separate bedgraph files for each population.

    Args:
        body (list): VCF body content
        samples (list): List of sample names from VCF header
        pop_map (dict): Mapping of sample names to populations
        bedgraph_base (str): Base path for output bedgraph files
    """
    # Create a mapping of sample indices to populations
    sample_indices = {pop: [] for pop in set(pop_map.values())}
    for i, sample in enumerate(samples):
        if sample in pop_map:
            pop = pop_map[sample]
            sample_indices[pop].append(i + 9)  # Add 9 to account for VCF fixed fields

    # Create a bedgraph file for each population
    bedgraph_files = {}
    for pop in sample_indices:
        bedgraph_path = f"{bedgraph_base}.{pop}.bedgraph"
        bedgraph_files[pop] = open(bedgraph_path, "w")
        bedgraph_files[pop].write("#chrom\tstart\tend\tindividuals_non_missing\n")

    # Process each row in the VCF
    max_pos = 0
    chrom = None
    for row in body:
        chrom = row[0]
        pos = int(row[1])
        max_pos = max(max_pos, pos)

        # Count non-missing genotypes for each population
        for pop, indices in sample_indices.items():
            genotypes = [row[i] for i in indices]
            non_missing_count = sum(1 for gt in genotypes if gt != ".|.")
            bedgraph_files[pop].write(
                f"{chrom}\t{pos - 1}\t{pos}\t{non_missing_count}\n"
            )

    # Close all bedgraph files
    for f in bedgraph_files.values():
        f.close()

    # Write genome file
    with open(Path(bedgraph_base).parent / "genome.txt", "w") as f:
        f.write(f"{chrom}\t{max_pos}")


def modify_genotypes(body, prop_missing_genotypes):
    rows = len(body)
    cols = len(body[0][9:])  # Assuming genotype data starts at index 9
    n_genotypes = rows * cols
    n_missing = int(prop_missing_genotypes * n_genotypes)

    mask = np.array([True] * n_missing + [False] * (n_genotypes - n_missing))
    np.random.shuffle(mask)

    mask = mask.reshape(rows, cols)
    for i, row in enumerate(body):
        genotypes = row[9:]
        for j in range(len(genotypes)):
            if mask[i, j]:
                genotypes[j] = ".|."
        row[9:] = genotypes

    return body


def modify_sites(body, prop_missing_sites):
    body_array = np.array(body)
    num_missing_sites = int(len(body) * prop_missing_sites)
    indices_to_remove = np.random.choice(len(body), num_missing_sites, replace=False)
    mask = np.ones(len(body), dtype=bool)
    mask[indices_to_remove] = False
    modified_body = body_array[mask].tolist()
    return modified_body


def process_genotypes(
    header, body, samples, pop_map, prop_missing_genotypes, output_vcf, bedgraph_base
):
    modified_body = modify_genotypes(body, prop_missing_genotypes)
    write_vcf(header, modified_body, output_vcf)
    create_bedgraph(modified_body, samples, pop_map, bedgraph_base)


def process_sites(
    header, body, samples, pop_map, prop_missing_sites, output_vcf, bedgraph_base
):
    modified_body = modify_sites(body, prop_missing_sites)
    write_vcf(header, modified_body, output_vcf)
    create_bedgraph(modified_body, samples, pop_map, bedgraph_base)


def main():
    from snakemake.script import snakemake

    input_vcf = snakemake.input["vcf"]
    pop_file = snakemake.input["populations"]
    output_gtsvcf = snakemake.output["gtsvcf"]
    bedgraph_gtsbase = str(snakemake.output["gtsbed1"]).rsplit(".", 2)[0]
    output_sitesvcf = snakemake.output["sitesvcf"]
    bedgraph_sitesbase = str(snakemake.output["sitesbed1"]).rsplit(".", 2)[0]
    prop = float(snakemake.params["prop"])

    # Parse the population file
    pop_map, unique_pops = parse_populations(pop_file)

    # Parse the input VCF
    header, body, samples = parse_vcf(input_vcf)

    # Start multiprocessing tasks
    processes = []

    # Process missing genotypes
    p1 = Process(
        target=process_genotypes,
        args=(header, body, samples, pop_map, prop, output_gtsvcf, bedgraph_gtsbase),
    )
    processes.append(p1)

    # Process missing sites
    p2 = Process(
        target=process_sites,
        args=(
            header,
            body,
            samples,
            pop_map,
            prop,
            output_sitesvcf,
            bedgraph_sitesbase,
        ),
    )
    processes.append(p2)

    # Start processes
    for p in processes:
        p.start()

    # Wait for processes to finish
    for p in processes:
        p.join()


if __name__ == "__main__":
    main()