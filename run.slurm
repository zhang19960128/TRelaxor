#!/bin/bash
#SBATCH -n 96
#SBATCH -p all
#SBATCH -t 20:00:00
/workspace/jiahaoz/PACKAGE_INSTALL/INTEL/mpi/2021.2.0/bin/mpirun -n 10 ./sa.x control.PT > simulated.out
