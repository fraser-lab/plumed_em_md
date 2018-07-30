import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--map_dat", help="full path to the GMM dat file for the map",required=True)
parser.add_argument("--n_proc", help="number of threads to use",required=True)
args = parser.parse_args()
map_plumed_dat_file = args.map_dat
n_proc = int(args.n_proc)


fout = open("launch_prep_3.sh","w")

fout.write("""#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./out3
#$ -e ./err3
#$ -l arch=linux-x64    # Specify architecture, required
#$ -l mem_free=1G       # Memory usage, required.  Note that this is per slot
#$ -pe smp {n_proc}         # Specify parallel environment and number of slots
#$ -R yes               # SGE host reservation, highly recommended
#$ -cwd                 # Current working directory
#$ -l h_rt=1:00:00

hostname
date

module load openmpi-1.8-x86_64
source /netapp/home/jaimefraser/plumed2/sourceme.sh
source /netapp/home/jaimefraser/gromacs-rdtscp/bin/GMXRC
export OMP_NUM_THREADS=1

python /netapp/home/jaimefraser/plumed_em_md/prep_plumed_3.py --map_dat {map_dat} --n_proc {n_proc}

date
""".format(map_dat=map_dat,n_proc=n_proc))
fout.close()

import os
os.system("qsub launch_prep_3.sh")
