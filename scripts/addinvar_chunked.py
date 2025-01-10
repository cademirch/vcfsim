import argparse
import pysam
import os
from multiprocessing import Pool
from tempfile import NamedTemporaryFile


def process_chunk(args):
    chrom, start, end, reference_fasta, variants_vcf, temp_output_path = args

    # Reopen files in the worker process
    ref_fasta = pysam.FastaFile(reference_fasta)
    chrom_seq = ref_fasta.fetch(chrom)

    vcf_in = pysam.VariantFile(variants_vcf)
    header = vcf_in.header

    # Add sample information to the header in the worker process
    with pysam.VariantFile(temp_output_path, "w", header=header) as temp_out:
        last_pos = start
        for rec in vcf_in.fetch(chrom, start - 1, end):
            if rec.pos > last_pos + 1:
                for pos in range(last_pos + 1, rec.pos):
                    if start <= pos <= end:
                        ref_base = chrom_seq[pos - 1]
                        invariant_record = create_invariant_record(
                            header, chrom, pos, ref_base, rec.samples
                        )
                        temp_out.write(invariant_record)

            temp_out.write(rec)
            last_pos = rec.pos

        # Add invariant sites after the last variant in this chunk
        if last_pos < end:
            for pos in range(last_pos + 1, end + 1):
                ref_base = chrom_seq[pos - 1]
                invariant_record = create_invariant_record(
                    header, chrom, pos, ref_base, rec.samples
                )
                temp_out.write(invariant_record)

    ref_fasta.close()
    vcf_in.close()


def add_invariant_sites_parallel(
    variants_vcf, reference_fasta, output_vcf, chunk_size=1000000
):
    # Open reference and VCF to get metadata
    ref_fasta = pysam.FastaFile(reference_fasta)
    vcf_in = pysam.VariantFile(variants_vcf)
    chrom = ref_fasta.references[0]
    chrom_len = ref_fasta.get_reference_length(chrom)
    ref_fasta.close()

    # Determine chunk boundaries
    chunks = [
        (chrom, start, min(start + chunk_size - 1, chrom_len))
        for start in range(1, chrom_len + 1, chunk_size)
    ]

    temp_files = []
    tasks = []

    # Create tasks for multiprocessing
    for start, end in [(start, end) for _, start, end in chunks]:
        temp_output = NamedTemporaryFile(delete=False, suffix=".vcf")
        temp_files.append(temp_output.name)

        tasks.append(
            (chrom, start, end, reference_fasta, variants_vcf, temp_output.name)
        )

    # Process chunks in parallel
    with Pool() as pool:
        pool.map(process_chunk, tasks)

    # Merge temporary files into the final output
    with pysam.VariantFile(output_vcf, "w", header=vcf_in.header) as vcf_out:
        for temp_file in temp_files:
            with pysam.VariantFile(temp_file) as temp_in:
                for rec in temp_in:
                    vcf_out.write(rec)
            os.remove(temp_file)

    # Close input files
    vcf_in.close()

def create_invariant_record(header, chrom, pos, ref_base, samples):
    # Create a new VCF record for an invariant site
    record = header.new_record(
        contig=chrom, start=pos - 1, stop=pos, alleles=(ref_base, ".")
    )
    record.filter.add("PASS")

    # Set genotypes to 0|0 for all samples
    for sample in samples:
        record.samples[sample]["GT"] = (0, 0)

    return record



if __name__ == "__main__":
    try:
        from snakemake.script import snakemake

        variants_vcf = snakemake.input["vcf"]
        reference_fasta = snakemake.input["ref"]
        output_vcf = snakemake.output["vcf"]
        chunk_size = snakemake.params["chunk_size"]

    except ImportError:
        parser = argparse.ArgumentParser(
            description="Add invariant sites to a variant-only VCF."
        )
        parser.add_argument("variants_vcf", help="Input variant-only VCF file")
        parser.add_argument("reference_fasta", help="Reference FASTA file")
        parser.add_argument("output_vcf", help="Output VCF file with invariant sites")
        parser.add_argument(
            "--chunk_size",
            type=int,
            default=1000000,
            help="Size of chunks to process in parallel (default: 1,000,000)",
        )

        args = parser.parse_args()
        variants_vcf = args.variants_vcf
        reference_fasta = args.reference_fasta
        output_vcf = args.output_vcf
        chunk_size = args.chunk_size

    add_invariant_sites_parallel(variants_vcf, reference_fasta, output_vcf, chunk_size)
