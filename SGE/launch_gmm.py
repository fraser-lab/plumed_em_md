import argparse, sys
parser = argparse.ArgumentParser()
parser.add_argument("--input_map", help="full path to the input map for simulation",required=True)
parser.add_argument("--n_proc", help="number of threads to use",required=True)
parser.add_argument("--gmconvert_location", help="location of custom gmconvert",default="/netapp/home/jaimefraser/gmconvert/gmconvert")
parser.add_argument("--sbgrid_gmconvert_location", help="location of sbgrid gmconvert",default="/programs/x86_64-linux/gmconvert/20180516/bin/gmconvert")
args = parser.parse_args()


n_proc = int(args.n_proc)
input_map = args.input_map
gmconvert_location = args.gmconvert_location
#/home/jfraser/gmconvert/gmconvert
sbgrid_gmconvert = args.sbgrid_gmconvert_location
#identical on qb3 and orion

fout = open("launch_gmm.sh","w")

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

python /netapp/home/jaimefraser/plumed_em_md/generate_gmm.py --n_proc {n_proc} --input_map {input_map} --gmconvert_location {gmconvert_location}  --sbgrid_gmconvert_location {sbgrid_gmconvert}

date
""".format(n_proc=n_proc,input_map=input_map, gmconvert_location=gmconvert_location, sbgrid_gmconvert=sbgrid_gmconvert))
fout.close()

import os
os.system("qsub launch_gmm.sh")
