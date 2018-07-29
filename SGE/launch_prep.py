import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--pdb", help="pdb file name to run",required=True)
parser.add_argument("--n_proc", help="number of threads to use",required=True,type=int)

args = parser.parse_args()
pdb = args.pdb
n_proc = int(args.n_proc)

fout = open("launch_prep.sh","w")

fout.write("""#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./out
#$ -e ./err
#$ -l arch=linux-x64    # Specify architecture, required
#$ -l mem_free=1G       # Memory usage, required.  Note that this is per slot
#$ -pe ompi {n_proc}            # Specify parallel environment and number of slots
#$ -R yes               # SGE host reservation, highly recommended
#$ -cwd                 # Current working directory
#$ -l h_rt=80:00:00


hostname
date

module load openmpi-1.8-x86_64
source /netapp/home/jaimefraser/plumed2/sourceme.sh
source /netapp/home/jaimefraser/gromacs-2016.5-bin/bin/GMXRC
export OMP_NUM_THREADS=1

python /netapp/home/jaimefraser/plumed_em_md/prep_directory.py
python /netapp/home/jaimefraser/plumed_em_md/prep_plumed.py --input_pdb {pdb} --n_proc {n_proc}

date
""".format(pdb=pdb,n_proc=n_proc))
fout.close()

import os
os.system("qsub launch_prep.sh")
