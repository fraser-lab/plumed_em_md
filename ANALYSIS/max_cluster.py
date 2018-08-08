from subprocess import check_output
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--xtc", help="trajectory xtc file",required=True)
args = parser.parse_args()
xtc = args.xtc


#GET NUMBER OF frames
out = check_output(["gmx_mpi check -f {xtc}".format(xtc=xtc)],shell=True)

print(out)
#divide into number of processors
#farm out RMSD calculations
#concatenate RMSD logs
#cluster
#pull out exemplars
#calcuate density in clusters by trajectory/replica
