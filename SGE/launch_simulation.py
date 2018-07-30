import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--n_proc", help="number of threads to use",required=True,type=int)
args = parser.parse_args()
n_proc = int(args.n_proc)

fout = open("launch_simulation.sh","w")

fout.write("""#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./out_sim
#$ -e ./err_sim
#$ -l arch=linux-x64    # Specify architecture, required
#$ -l mem_free=1G       # Memory usage, required.  Note that this is per slot
#$ -pe ompi_onehost {n_proc}            # Specify parallel environment and number of slots
#$ -R yes               # SGE host reservation, highly recommended
#$ -cwd                 # Current working directory
#$ -l h_rt=300:00:00
#$ -V

hostname
date

module load openmpi-1.8-x86_64
source /netapp/home/jaimefraser/plumed2/sourceme.sh
source /netapp/home/jaimefraser/gromacs-rdtscp/bin/GMXRC
export OMP_NUM_THREADS=1

mpirun -n {n_proc} gmx_mpi mdrun -plumed -multi 8 -ntomp 1 -cpi state

date
""".format(n_proc=n_proc))
fout.close()

import os
os.system("qsub launch_simulation.sh")
