from subprocess import check_output
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--xtc", help="trajectory xtc file",required=True)
args = parser.parse_args()
xtc = args.xtc


#GET NUMBER OF frames
os.system("gmx_mpi check -f {xtc} 2> temp.txt".format(xtc=xtc))
out = check_output(["grep Box temp.txt"],shell=True)
frames = int(out.split()[-1])
print(frames)
#divide into number of processors
#farm out RMSD calculations
#concatenate RMSD logs
#cluster
#pull out exemplars
#calcuate density in clusters by trajectory/replica
