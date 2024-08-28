#!/bin/bash

# Replace 'input.csv' with the path to your CSV file
csv_file="../atb_screening_data_top100.csv"

# Use awk to skip the header and print the second column
awk -F, 'NR>1 {print $2}' "$csv_file" | while read -r element; do
    # Pass each element to the Python script
    python caching_pipeline.py "$element"
done
