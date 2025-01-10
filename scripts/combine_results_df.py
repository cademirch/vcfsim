import pandas as pd
from snakemake.script import snakemake

clam_dfs = [pd.read_csv(f, sep="\t") for f in snakemake.input.clam]
pixy_dfs = [pd.read_csv(f, sep="\t") for f in snakemake.input.pixy]
prime_dfs = [pd.read_csv(f, sep="\t") for f in snakemake.input.prime]

# Add source index and proportion columns
clam_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(clam_dfs)]
pixy_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(pixy_dfs)]
prime_dfs_with_index = [df.assign(source_idx=i) for i, df in enumerate(prime_dfs)]

# Concatenate and save to CSV
pd.concat(clam_dfs_with_index).to_csv(snakemake.output.clam, index=False)
pd.concat(pixy_dfs_with_index).to_csv(snakemake.output.pixy, index=False)
pd.concat(prime_dfs_with_index).to_csv(snakemake.output.prime, index=False)
