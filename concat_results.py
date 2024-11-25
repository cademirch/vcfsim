import pandas as pd
from pathlib import Path
import itertools


def read_truth_data(truth_file):
    """
    Reads a truth data file and creates a DataFrame with chrom, start, end, and pi columns.

    Parameters:
        truth_file (str): Path to the truth data file containing start-end and pi arrays.

    Returns:
        pd.DataFrame: DataFrame with columns 'chrom', 'start', 'end', and 'pi'.
    """
    with open(truth_file, "r") as f:
        # Read the start-end array
        lines = f.readlines()

    positions = eval(lines[0].strip())
    # Read the pi array
    pi_values = [l.strip().strip("[]").split() for l in lines[1:]]
    pi_values = [float(p) for p in itertools.chain.from_iterable(pi_values)]

    if len(positions) - 1 != len(pi_values):
        raise ValueError("Mismatch between positions and pi values in the truth file.")

    # Adjust positions to start at 1 and ensure inclusivity
    start_positions = positions[:-1]
    end_positions = [pos - 1 for pos in positions[1:]]  # Make ends inclusive

    # Build the DataFrame
    data = {
        "chrom": ["chr1"] * (len(positions) - 1),
        "start": start_positions,
        "end": end_positions,
        "pi": pi_values,
    }

    return pd.DataFrame(data)


results_dir = Path("./results")
clam_res = Path("clam_pi.tsv")
pixy_res = Path("pixy_pi.txt")
sim_res = Path("sim/pi_windows.txt")

clam_dfs = []
pixy_dfs = []
sim_dfs = []

for i in range(1, 400):
    r = Path(results_dir, str(i))
    if r.exists():
        clam_file = Path(r, clam_res)
        pixy_file = Path(r, pixy_res)
        sim_file = Path(r, sim_res)
        files = [clam_file, pixy_file, sim_file]
        if all([f.exists() for f in files]):
            clam_pd = pd.read_csv(clam_file, sep="\t")
            pixy_pd = pd.read_csv(pixy_file, sep="\t")
            sim_pd = read_truth_data(sim_file)
            clam_dfs.append(clam_pd)
            pixy_dfs.append(pixy_pd)
            sim_dfs.append(sim_pd)

clam_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(clam_dfs)]
pixy_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(pixy_dfs)]
sim_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(sim_dfs)]


all_clam = pd.concat(clam_dfs_with_index, ignore_index=True)
all_pixy = pd.concat(pixy_dfs_with_index, ignore_index=True)
all_sim = pd.concat(sim_dfs_with_index, ignore_index=True)

all_clam.to_csv("all_clam.csv")
all_sim.to_csv("all_sim.csv")
