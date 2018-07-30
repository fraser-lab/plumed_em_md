fout = open("launch_prep_2.sh","w")

fout.write("""#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./out2
#$ -e ./err2
#$ -l arch=linux-x64    # Specify architecture, required
#$ -l mem_free=1G       # Memory usage, required.  Note that this is per slot
#$ -R yes               # SGE host reservation, highly recommended
#$ -cwd                 # Current working directory
#$ -l h_rt=00:14:00

hostname
date

source /netapp/home/jaimefraser/plumed2/sourceme.sh
source /netapp/home/jaimefraser/gromacs-rdtscp/bin/GMXRC
export OMP_NUM_THREADS=1

python /netapp/home/jaimefraser/plumed_em_md/prep_plumed_2.py

date
""")
fout.close()

import os
os.system("qsub launch_prep_2.sh")
