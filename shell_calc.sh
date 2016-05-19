#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=12
#PBS -l vmem=60gb
#PBS -m a
#PBS -M dilyn.fullerton@att.net
#PBS -j oe 

NTHREADS=12
export OMP_NUM_THREADS=${NTHREADS}

cd ${PBS_O_WORKDIR}

# Enter the job you want to be done on the cluster.
# helium
python shell_calc.py -v 1 1 t 4 10 2
# oxygen
# python shell_calc.py -v 2 1 t 16 24 8
