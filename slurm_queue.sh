#!/bin/bash

### Nombre de la tarea
#SBATCH --job-name=tiny_md_slurm_202335

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time 01:00:00
. /etc/profile

#export OMP_NUM_THREADS=4

#srun perf stat ./tiny_md > hello.log
srun ./extract_data.sh $1

