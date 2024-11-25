#!/bin/bash

# Activate your Snakemake conda environment

conda activate snakemake

# Function to run a single Snakemake command
run_snakemake() {
    local index=$1
    echo "Running Snakemake for index $index..."
    snakemake -d results/$index --use-conda --conda-prefix /home/cade/vcfsim/.snakemake/conda --cores 10 &> logs/$index.log
    if [ $? -eq 0 ]; then
        echo "Index $index completed successfully."
    else
        echo "Index $index failed."
    fi
}

# Export the function so it can be used in parallel
export -f run_snakemake

# Run the Snakemake commands in parallel
seq 1 1000 | parallel -j 10 run_snakemake

echo "All jobs are finished."
