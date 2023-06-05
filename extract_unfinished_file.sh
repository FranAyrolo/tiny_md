#!/bin/bash

# Valores que podemos usar de N = 256, 500, 864, 1372, 2048, 2916, 4000, 5324, 6912

# Initialize variables
input_file="resultados/input.log"
# Output CSV file
csv_file="resultados/output.csv"

# Write CSV header
echo "N value, thread number, Iteration number, Number of instructions, Ins per cycle, Percentage of misses, Seconds time elapsed" > "$csv_file"

# Extract data from text file
while IFS= read -r line
do
    # Extract N value
    if [[ $line =~ ^N=\ ([0-9]+) ]]; then
        n_value=${BASH_REMATCH[1]}
    fi
    # Extract thread number
    if [[ $line =~ Num.\ Threads\ ([0-9]+) ]]; then
        thread_num=${BASH_REMATCH[1]}
    fi

    # Extract iteration number
    if [[ $line =~ Iteration\ ([0-9]+) ]]; then
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
        echo "$n_value, $thread_num, $iteration,$instructions,$ins_per_cycle,$misses,$time_elapsed" >> "$csv_file"
    fi 
done < $input_file
