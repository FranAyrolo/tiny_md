#!/bin/bash

# Valores que podemos usar de N = 256, 500, 864, 1372, 2048, 2916, 4000, 5324, 6912

# Initialize variables
num_iterations=5
temp_file="resultados/temp.log"
output_file="resultados/multicore_N1372.log"
# Output CSV file
csv_file="resultados/multicore_N1372.csv"

numThreads=(1 2 4 8)

# Loop over the command
for threads in "${numThreads[@]}"; do
    export OMP_NUM_THREADS=$threads
    make clean >> /dev/null
    make >> /dev/null
    for((i=1; i <= num_iterations; i++)); do
        echo "Num. Threads = $threads - Iteration $i"
        echo "Num. Threads $threads - Iteration $i" >> $output_file
        perf stat -o $temp_file ./tiny_md
        cat $temp_file >> $output_file
    done
done

# Cleanup temporary file
rm $temp_file


# Write CSV header
echo "Iteration number, Number of instructions, Ins per cycle, Percentage of misses, Seconds time elapsed" > "$csv_file"

# Extract data from text file
while IFS= read -r line
do
    # Extract iteration number
    if [[ $line =~ ^Iteration\ ([0-9]+) ]]; then
        iteration=${BASH_REMATCH[1]}
    fi

    # Extract number of instructions and ins per cycle
    if [[ $line =~ instructions ]]; then
        instructions=$(echo $line | awk '{print $1}' | tr -d ',')
        ins_per_cycle=$(echo $line | awk '{print $4}')
    fi

    # Extract percentage of misses
    if [[ $line =~ ([0-9]+\.[0-9]+)%\ of\ all\ branches ]]; then
        misses=${BASH_REMATCH[1]}
    fi

    # Extract secons time elapsed
    if [[ $line =~ ([0-9]+\.[0-9]+)\ seconds\ time\ elapsed ]]; then
        time_elapsed=${BASH_REMATCH[1]}
        
        # Write data to CSV file
        echo "$iteration,$instructions,$ins_per_cycle,$misses,$time_elapsed" >> "$csv_file"
    fi
    
done < $output_file
