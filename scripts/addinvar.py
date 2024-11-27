import argparse
import pysam


def add_invariant_sites(variants_vcf, reference_fasta, output_vcf):
    # Open the variant VCF file and reference FASTA
    vcf_in = pysam.VariantFile(variants_vcf)
    ref_fasta = pysam.FastaFile(reference_fasta)

    # Open the output VCF file and copy the header
    header = vcf_in.header

    vcf_out = pysam.VariantFile(output_vcf, "w", header=header)

    for chrom in ref_fasta.references:
        chrom_len = ref_fasta.get_reference_length(chrom)
        last_pos = 0

        # Iterate over the variant VCF records and add invariant sites between variants
        for rec in vcf_in:
            # Add invariant sites between last variant and current variant
            if rec.pos > last_pos + 1:
                for pos in range(last_pos + 1, rec.pos):
                    ref_base = ref_fasta.fetch(
                        chrom, pos - 1, pos
                    )  # Fetch reference base (0-based)
                    invariant_record = create_invariant_record(
                        vcf_in.header, chrom, pos, ref_base, rec.samples
                    )
                    vcf_out.write(invariant_record)

            # Write the current variant record
            vcf_out.write(rec)
            last_pos = rec.pos

        # Add invariant sites after the last variant to the end of the chromosome
        if last_pos < chrom_len:
            for pos in range(last_pos + 1, chrom_len + 1):
                ref_base = ref_fasta.fetch(chrom, pos - 1, pos)
                invariant_record = create_invariant_record(
                    vcf_in.header, chrom, pos, ref_base, rec.samples
                )
                vcf_out.write(invariant_record)

    vcf_in.close()
    ref_fasta.close()
    vcf_out.close()


def create_invariant_record(header, chrom, pos, ref_base, samples):
    # Create a new VCF record for an invariant site
    record = header.new_record(
        contig=chrom, start=pos - 1, stop=pos, alleles=(ref_base, ".")
    )
    record.qual = None
    record.filter.add("PASS")

    # Set the genotype of invariant sites to 0|0 for all samples
    for sample in samples:
        record.samples[sample]["GT"] = (0, 0)  # 0|0 for invariant sites

    return record


if __name__ == "__main__":
    try:
        from snakemake.script import snakemake

        variants_vcf = snakemake.input["vcf"]
        reference_fasta = snakemake.input["ref"]
        output_vcf = snakemake.output["vcf"]

    except ImportError:
        parser = argparse.ArgumentParser(
            description="Add invariant sites to a variant-only VCF."
        )
        parser.add_argument("variants_vcf", help="Input variant-only VCF file")
        parser.add_argument("reference_fasta", help="Reference FASTA file")
        parser.add_argument("output_vcf", help="Output VCF file with invariant sites")

        args = parser.parse_args()
        variants_vcf = args.variants_vcf
        reference_fasta = args.reference_fasta
        output_vcf = args.output_vcf

    add_invariant_sites(variants_vcf, reference_fasta, output_vcf)
