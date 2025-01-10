import random
import argparse
import gzip
from pathlib import Path
from multiprocessing import Process
import numpy as np


def parse_vcf(input_vcf):
    with (
        gzip.open(input_vcf, "rt")
        if input_vcf.endswith(".gz")
        else open(input_vcf, "r")
    ) as vcf:
        header = []
        body = []
        for line in vcf:
            if line.startswith("#"):
                header.append(line.strip())
            else:
                body.append(line.strip().split("\t"))
        return header, body


# Function to write VCF
def write_vcf(header, body, output_vcf):
    with open(output_vcf, "w") as vcf:
        vcf.write("\n".join(header) + "\n")
        for row in body:
            vcf.write("\t".join(row) + "\n")


# Function to create the BEDGRAPH file
def create_bedgraph(body, bedgraph_file):
    with open(bedgraph_file, "w") as bed:
        bed.write("#chrom\tstart\tend\tindividuals_non_missing\n")
        for row in body:
            chrom = row[0]
            pos = int(row[1])
            genotypes = row[9:]
            non_missing_count = sum(1 for gt in genotypes if gt != ".|.")
            bed.write(f"{chrom}\t{pos - 1}\t{pos}\t{non_missing_count}\n")
    with open(
        Path(
            Path(bedgraph_file).parents[0],
            "genome.txt",
        ),
        "w",
    ) as f:
        f.write(f"{chrom}\t{pos}")


# Function to modify genotypes in VCF
def modify_genotypes(body, prop_missing_genotypes):
    rows = len(body)
    cols = len(body[0][9:])  # Assuming genotype data starts at index 9
    n_genotypes = rows * cols
    n_missing = int(prop_missing_genotypes * n_genotypes)

    # Create a mask with `n_missing` True values (missing data)
    mask = np.array([True] * n_missing + [False] * (n_genotypes - n_missing))
    np.random.shuffle(mask)

    # Apply mask across the entire matrix
    mask = mask.reshape(rows, cols)
    for i, row in enumerate(body):
        genotypes = row[9:]
        for j in range(len(genotypes)):
            if mask[i, j]:
                genotypes[j] = ".|."
        row[9:] = genotypes

    return body


def modify_sites(body, prop_missing_sites):
    """
    Removes a proportion of sites (rows) from the VCF body using NumPy.

    Args:
        body (list): The VCF file contents split into lines, with header lines excluded.
        prop_missing_sites (float): Proportion of sites to remove (between 0 and 1).

    Returns:
        list: Modified VCF body with the specified proportion of rows removed.
    """
    # Convert body to NumPy array for efficient indexing
    body_array = np.array(body)

    # Calculate the number of rows to remove
    num_missing_sites = int(len(body) * prop_missing_sites)

    # Randomly select indices to remove
    indices_to_remove = np.random.choice(len(body), num_missing_sites, replace=False)

    # Create a mask for rows to keep
    mask = np.ones(len(body), dtype=bool)
    mask[indices_to_remove] = False

    # Apply mask to retain only the rows not selected for removal
    modified_body = body_array[mask].tolist()

    return modified_body


# Worker function for processing missing genotypes
def process_genotypes(header, body, prop_missing_genotypes, output_vcf, bedgraph_file):
    modified_body = modify_genotypes(body, prop_missing_genotypes)
    write_vcf(header, modified_body, output_vcf)
    create_bedgraph(modified_body, bedgraph_file)


# Worker function for processing missing sites
def process_sites(header, body, prop_missing_sites, output_vcf, bedgraph_file):
    modified_body = modify_sites(body, prop_missing_sites)
    write_vcf(header, modified_body, output_vcf)
    create_bedgraph(modified_body, bedgraph_file)


def main():
    from snakemake.script import snakemake

    input_vcf = snakemake.input["vcf"]
    output_gtsvcf = snakemake.output["gtsvcf"]
    bedgraph_gtsfile = snakemake.output["gtsbed"]
    output_sitesvcf = snakemake.output["sitesvcf"]
    bedgraph_sitesfile = snakemake.output["sitesbed"]
    prop = float(snakemake.params["prop"])

    # Parse the input VCF
    header, body = parse_vcf(input_vcf)

    # Start multiprocessing tasks
    processes = []

    # Process missing genotypes
    p1 = Process(
        target=process_genotypes,
        args=(header, body, prop, output_gtsvcf, bedgraph_gtsfile),
    )
    processes.append(p1)

    # Process missing sites
    p2 = Process(
        target=process_sites,
        args=(header, body, prop, output_sitesvcf, bedgraph_sitesfile),
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
