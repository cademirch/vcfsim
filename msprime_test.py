import msprime
import tskit
import argparse
from pathlib import Path
import numpy as np


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


def simulate(seqlen, ne, mu, sample_size, seed, outdir, windows):

    ts = msprime.sim_ancestry(
        sample_size, sequence_length=seqlen, random_seed=seed, population_size=ne
    )
    mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)
    reference_sequence = tskit.random_nucleotides(seqlen, seed=seed)
    num_mutations = mts.num_mutations

    ref_fasta_entry = "\n".join([">chr1", reference_sequence])

    Path(outdir, "ref.fa").write_text(ref_fasta_entry)

    alignments = mts.alignments(reference_sequence=reference_sequence)

    with open(Path(outdir, "sim.vcf"), "w") as f:
        mts.write_vcf(f, position_transform=lambda x: np.fmax(1, x))

    for i, (hap1, hap2) in enumerate(zip(alignments, alignments)):
        fp = Path(outdir, f"sample_{i}.fa")
        fasta_entry = "\n".join([">hap1", hap1, ">hap2", hap2])
        fp.write_text(fasta_entry)
    Path(outdir, "sim_results.txt").write_text(
        f"{seqlen=}, {ne=}, {mu=}, {sample_size=}, {seed=}, {num_mutations=}"
    )
    if windows:
        pi_windows = mts.diversity(windows=windows)
        Path(outdir, "pi_windows.txt").write_text(
            "\n".join([str(windows), str(pi_windows)])
        )


def main():
    try:
        from snakemake.script import snakemake

        ne = snakemake.params["ne"]  # noqa: F821
        mu = snakemake.params["mu"]  # noqa: F821
        sample_size = snakemake.params["sample_size"]  # noqa: F821
        seed = snakemake.params["seed"]  # noqa: F821
        outdir = snakemake.params["outdir"]  # noqa: F821
        seqlen = snakemake.params["seqlen"]  # noqa: F821

        windows = snakemake.params.get("windows", None)

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

        args = parser.parse_args()

        outdir = args.outdir
        seqlen = args.seqlen
        ne = args.ne
        mu = args.mu
        sample_size = args.sample_size
        seed = args.seed
        windows = args.windows

    if windows is not None:
        windows = generate_windows(seqlen, windows)

    simulate(seqlen, ne, mu, sample_size, seed, outdir, windows)


if __name__ == "__main__":
    main()
