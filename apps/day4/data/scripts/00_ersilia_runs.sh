#!/bin/bash

# Function to run the ersilia pipeline for a given EOS identifier
run_ersilia_pipeline() {
    eos_id=$1
    input_file=$2
    output_folder=$3

    # Fetch the model
    ersilia -v fetch "$eos_id"

    # Serve the model
    ersilia serve "$eos_id"

    # Run the model on the input file
    output_file="${output_folder}/${eos_id}_all_mols.csv"
    ersilia run -i "$input_file" -o "$output_file"

    # Close the model
    ersilia close
}

# Main script
input_file="../processed/all_mols.csv"
output_folder="../processed"

# Array of EOS identifiers
# eos_identifiers=("eos7kpb" "eos4zfy" "eos4rta" "eos7d58" "eos39dp")  # Add all the EOS identifiers you want
eos_identifiers=("eos39dp")

# Loop through each EOS identifier and run the pipeline
for eos_id in "${eos_identifiers[@]}"; do
    run_ersilia_pipeline "$eos_id" "$input_file" "$output_folder"
done
