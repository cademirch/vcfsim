import argparse
import pysam


def add_invariant_sites(input_vcf, output_vcf):
    # Open the input VCF for reading
    vcf_in = pysam.VariantFile(input_vcf, "r")

    # Open the output VCF for writing, copying the header from the input VCF
    vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

    # Initialize variables
    last_position = None
    chrom = None

    # Iterate through the records in the VCF
    for record in vcf_in:
        # If this is the first record, initialize chromosome and last position
        if last_position is None:
            chrom = record.chrom
            last_position = record.pos
        else:
            # Insert invariant sites between last_position and the current record's position
            for pos in range(last_position + 1, record.pos):
                # Create a new record for the invariant site
                new_record = vcf_out.new_record()
                new_record.chrom = chrom
                new_record.pos = pos
                new_record.ref = "N"  # Placeholder for invariant reference
                new_record.alts = ["."]  # No alternate alleles
                new_record.qual = None
                new_record.filter.add("PASS")
                new_record.info.clear()  # Clear any INFO fields
                new_record.format["GT"] = ["0|0"] * len(record.samples)
                vcf_out.write(new_record)

            # Update last position to the current record's position
            last_position = record.pos

        # Write the current record to the output VCF
        vcf_out.write(record)

    # Close input and output VCF files
    vcf_in.close()
    vcf_out.close()


def main():
    try:
        from snakemake.script import snakemake

        inp = snakemake.input[0]
        out = snakemake.output[0]

    except ImportError:
        # Set up argument parser
        parser = argparse.ArgumentParser(
            description="Add invariant sites to a VCF file."
        )
        parser.add_argument("input_vcf", help="Path to the input VCF file")
        parser.add_argument(
            "output_vcf", help="Path to the output VCF file with invariant sites added"
        )

        # Parse arguments
        args = parser.parse_args()
        inp = args.input_vcf
        out = args.output_vcf
    # Call the function with provided arguments
    add_invariant_sites(inp, out)

    print(f"Invariant sites added. Output written to {args.output_vcf}")


if __name__ == "__main__":
    main()
