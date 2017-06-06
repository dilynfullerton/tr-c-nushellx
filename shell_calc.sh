#!/bin/bash
#PBS -l walltime=06:00:00
#PBS -l nodes=1:ppn=12
#PBS -l vmem=60gb
#PBS -m a
#PBS -j oe 

NTHREADS=12
export OMP_NUM_THREADS=${NTHREADS}

cd ${PBS_O_WORKDIR}

# Enter the job you want to be done on the cluster.

# helium
# python shell_calc.py 1 1 t 4 14 2

# lithium
python shell_calc.py 1 2 t 6 16 3

# n=z
# python shell_calc.py 1 1 t 4 16 2 8

# oxygen
# python shell_calc.py 2 1 t 16 28 8;
