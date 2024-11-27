import random
import argparse
import gzip
from pathlib import Path


# Function to parse VCF files
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
            bed.write(f"{chrom}\t{pos-1}\t{pos}\t{non_missing_count}\n")
    with open(
        Path(
            Path(bedgraph_file).parents[0],
            "genome.txt",
        ),
        "w",
    ) as f:
        f.write(f"{chrom}\t{pos}")


# Function to randomly modify VCF
def modify_vcf(body, prop_missing_sites, prop_missing_genotypes):
    modified_body = []

    for row in body:
        if random.random() < prop_missing_sites:
            continue  # Skip the site entirely

        genotypes = row[9:]
        for i in range(len(genotypes)):
            if random.random() < prop_missing_genotypes:
                genotypes[i] = ".|."

        row[9:] = genotypes
        modified_body.append(row)

    return modified_body


def main():
    try:
        from snakemake.script import snakemake

        input_vcf = snakemake.input["vcf"]
        prop_missing_sites = snakemake.params["prop_missing_sites"]
        prop_missing_genotypes = snakemake.params["prop_missing_genotypes"]
        output_vcf = snakemake.output["vcf"]
        bedgraph_file = snakemake.output["bed"]

    except ImportError:
        parser = argparse.ArgumentParser(
            description="Modify VCF by introducing missing sites and genotypes."
        )
        parser.add_argument("-i", "--input_vcf", required=True, help="Input VCF file.")
        parser.add_argument(
            "-o", "--output_vcf", required=True, help="Output VCF file."
        )
        parser.add_argument(
            "-b", "--bedgraph_file", required=True, help="Output BEDGRAPH file."
        )
        parser.add_argument(
            "--prop_missing_sites",
            type=float,
            required=True,
            help="Proportion of missing sites.",
        )
        parser.add_argument(
            "--prop_missing_genotypes",
            type=float,
            required=True,
            help="Proportion of missing genotypes per site.",
        )

        args = parser.parse_args()
        input_vcf = args.input_vcf
        output_vcf = args.output_vcf
        prop_missing_genotypes = args.prop_missing_genotypes
        prop_missing_sites = args.prop_missing_sites
        bedgraph_file = args.bedgraph_file

    # Parse the input VCF
    header, body = parse_vcf(input_vcf)

    # Modify the VCF
    modified_body = modify_vcf(body, prop_missing_sites, prop_missing_genotypes)

    # Write the modified VCF
    write_vcf(header, modified_body, output_vcf)

    # Create the BEDGRAPH file
    create_bedgraph(modified_body, bedgraph_file)


if __name__ == "__main__":
    main()
