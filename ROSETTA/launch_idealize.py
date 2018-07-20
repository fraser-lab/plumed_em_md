import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--pdb", help="pdb file name to idealize",required=True)
args = parser.parse_args()
pdb = args.pdb

fout = open("launch_idealize.sh","w")

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

source /programs/sbgrid.shrc
idealize_jd2.linuxgccrelease -database /netapp/home/jaimefraser/database -in::path ./ -in::file::s {pdb}  -no_optH -out::path ./ -out::path::pdb ./ -chainbreaks  -overwrite
date
""".format(pdb=pdb))
fout.close()

import os
os.system("qsub launch_idealize.sh")
