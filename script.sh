#!/bin/bash

# Initialize variables
num_iterations=10
temp_file="resultados/temp.log"
output_file="resultados/gcc_O0.log"

echo "gcc -O0" >> $output_file

# Loop over the command
for((i = 1; i <= num_iterations; i++)); do
    echo "Iteration $i"
    sudo perf stat -o $temp_file ./tiny_md
    cat $temp_file >> $output_file
done

# Cleanup temporary file
rm $temp_file
