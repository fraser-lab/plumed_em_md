import argparse, sys
parser = argparse.ArgumentParser()
parser.add_argument("--input_map", help="full path to the input map for simulation",required=True)
parser.add_argument("--n_proc", help="number of threads to use",required=True)
parser.add_argument("--gmconvert_location", help="location of custom gmconvert",default="/home/jfraser/gmconvert/gmconvert")
parser.add_argument("--sbgrid_gmconvert_location", help="location of sbgrid gmconvert",default="/programs/x86_64-linux/gmconvert/20180516/bin/gmconvert")
args = parser.parse_args()

n_proc = int(args.n_proc)
input_map = args.input_map
gmconvert_location = args.gmconvert_location
sbgrid_gmconvert = args.sbgrid_gmconvert_location

fout = open("launch_prep.sh","w")

fout.write("""#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./out
#$ -e ./err
#$ -l arch=linux-x64    # Specify architecture, required
#$ -l mem_free=1G       # Memory usage, required.  Note that this is per slot
#$ -pe smp {n_proc}            # Specify parallel environment and number of slots
#$ -R yes               # SGE host reservation, highly recommended
#$ -cwd                 # Current working directory
#$ -l h_rt=80:00:00


hostname
date

python /netapp/home/jaimefraser/plumed_em_md/launch_gmm.py --n_proc {nproc} --input_map {input_map} --gmconvert_location /netapp/home/jaimefraser/gmconvert/gmconvert

date
""".format(pdb=pdb))
fout.close()

import os
os.system("qsub launch_prep.sh")
