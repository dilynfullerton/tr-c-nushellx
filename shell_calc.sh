#!/bin/bash
#PBS -l walltime=00:20:00
#PBS -l nodes=1:ppn=12
#PBS -l vmem=60gb
#PBS -m ae
#PBS -j oe 

NTHREADS=12
export OMP_NUM_THREADS=$NTHREADS

cd $PBS_O_WORKDIR

# Enter the job you want to be done on the cluster.
# helium
# python shell_calc.py -f 1 1 t 4 10 2
# oxygen
python shell_calc.py -f 2 1 t 16 24 8
