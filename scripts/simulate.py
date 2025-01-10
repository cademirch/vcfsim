import msprime
import tskit
import argparse
from pathlib import Path
import numpy as np
import concurrent.futures
import datetime
import pysam



def generate_windows(sequence_length, window_size):
    """
    Generate a list of window boundaries for a sequence of given length and window size.

    Args:
        sequence_length (int): The total length of the sequence.
        window_size (int): The size of each window.

    Returns:
        list: A list of increasing window boundaries starting with 0 and ending with sequence_length.
    """
    windows = list(range(0, sequence_length, window_size))
    if windows[-1] != sequence_length:
        windows.append(sequence_length)
    return windows


def write_vcf_for_sample(i, outdir, mts):
    """Helper function to write a VCF for a single sample."""
    with open(Path(outdir, f"sample_{i}.vcf"), "w") as f:
        mts.write_vcf(
            f,
            position_transform=lambda x: np.fmax(1, x),
            individuals=[i],
            contig_id="chr1",
        )


def simulate(
    seqlen, ne, mu, sample_size, seed, outdir, windows, max_workers, write_indv_vcfs
):
    # Simulate ancestry and mutations
    ts = msprime.sim_ancestry(
        sample_size,
        sequence_length=seqlen,
        random_seed=seed,
        population_size=ne,
        ploidy=2,
    )
    mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed, discrete_genome=True)
    reference_sequence = tskit.random_nucleotides(seqlen, seed=seed)
    num_mutations = mts.num_mutations

    # Write the reference sequence to a FASTA file
    ref_fasta_entry = "\n".join([">chr1", reference_sequence])
    Path(outdir, "ref.fa").write_text(ref_fasta_entry)

    if write_indv_vcfs:
        # Parallelize writing individual sample VCFs using ProcessPoolExecutor
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=max_workers
        ) as executor:
            futures = [
                executor.submit(write_vcf_for_sample, i, outdir, mts)
                for i in range(sample_size)
            ]
            for future in concurrent.futures.as_completed(futures):
                future.result()  # Ensures any exceptions are raised

    # Write the combined VCF for all samples
    all_vcf_path = Path(outdir, "all.vcf")
    with open(all_vcf_path, "w") as f:
        mts.write_vcf(f, position_transform=lambda x: np.fmax(1, x), contig_id="chr1")

    pysam.tabix_compress(
        all_vcf_path, all_vcf_path.with_suffix(all_vcf_path.suffix + ".gz")
    )
    pysam.tabix_index(
        str(all_vcf_path.with_suffix(all_vcf_path.suffix + ".gz")), preset="vcf"
    )

    # Calculate summary statistics
    pi = mts.diversity()
    theta = 4 * ne * mu
    Path(outdir, "sim_results.txt").write_text(
        f"{seqlen=}, {ne=}, {mu=}, {sample_size=}, {seed=}, {num_mutations=}, {pi=}, {theta=}"
    )

    # If windows are specified, calculate diversity per window and write to file
    if windows:
        pi_windows = mts.diversity(windows=windows)
        outfile = "msprime_pi_windows.tsv"
        start_positions = windows[:-1]
        end_positions = [pos - 1 for pos in windows[1:]]
        chrom = ["chr1"] * (len(windows) - 1)
        lines = "\n".join(
            [
                f"{chrom}\t{start}\t{end}\t{pi}"
                for chrom, start, end, pi in zip(
                    chrom, start_positions, end_positions, pi_windows
                )
            ]
        )

        with open(outfile, "w") as f:
            f.write("chrom\tstart\tend\tpi\n")
            f.writelines(lines)
def main():
    try:
        from snakemake.script import snakemake

        ne = snakemake.params["ne"]  # noqa: F821
        mu = snakemake.params["mu"]  # noqa: F821
        sample_size = snakemake.params["sample_size"]  # noqa: F821
        seed = int(snakemake.params["seed"])  # noqa: F821
        outdir = snakemake.params["outdir"]  # noqa: F821
        seqlen = snakemake.params["seqlen"]  # noqa: F821

        windows = snakemake.params.get("windows", None)
        write_indv_vcfs = snakemake.params.get("indv_vcfs", True)
        max_workers = snakemake.threads

    except ImportError:
        parser = argparse.ArgumentParser()

        parser.add_argument("--seqlen", type=int, help="Sequence length", required=True)
        parser.add_argument(
            "--ne", type=float, help="Effective population size", required=True
        )
        parser.add_argument("--mu", type=float, help="Mutation rate", required=True)
        parser.add_argument(
            "--sample-size", type=int, help="Sample size", required=True
        )
        parser.add_argument("--seed", type=int, help="Random seed", required=True)
        parser.add_argument(
            "--outdir", type=str, help="Output directory", required=True
        )
        parser.add_argument(
            "--windows", type=int, help="bp size of windows", required=False
        )
        parser.add_argument("--threads", type=int, help="num threads", required=True)
        parser.add_argument(
            "--indv-vcfs",
            type=bool,
            help="write indvidual vcfs",
            required=False,
            default=True,
        )

        args = parser.parse_args()

        outdir = args.outdir
        seqlen = args.seqlen
        ne = args.ne
        mu = args.mu
        sample_size = args.sample_size
        seed = args.seed
        windows = args.windows
        max_workers = args.threads
        write_indv_vcfs = args.indv_vcfs

    if windows is not None:
        windows = generate_windows(seqlen, windows)

    simulate(
        seqlen, ne, mu, sample_size, seed, outdir, windows, max_workers, write_indv_vcfs
    )


if __name__ == "__main__":
    main()
