import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--pdb", help="pdb file name to run",required=True)
args = parser.parse_args()
pdb = args.pdb

fout = open("launch_prep.sh","w")

fout.write("""#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./out
#$ -e ./err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=1G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=80:00:00

hostname
date

source /netapp/home/jaimefraser/plumed2/sourceme.sh
source /netapp/home/jaimefraser/gromacs-2016.5-bin/bin/GMXRC
/netapp/home/jaimefraser/plumed_em_md/prep_directory.py
/netapp/home/jaimefraser/plumed_em_md/prep_plumed.py --pdb

date
""".format(pdb=pdb))
fout.close()

import os
os.system("qsub launch_prep.sh")
