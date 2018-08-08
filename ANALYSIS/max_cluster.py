from subprocess import check_output
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--xtc", help="trajectory xtc file",required=True)
parser.add_argument("--nproc", help="number of cores",type=int,required=True)
args = parser.parse_args()
xtc = args.xtc
nproc = args.nproc

#GET NUMBER OF frames
os.system("gmx_mpi check -f {xtc} 2> temp.txt".format(xtc=xtc))
out = check_output(["grep Box temp.txt"],shell=True)
frames = int(out.split()[-1])

#divide into number of processors
#max number of calculations k1 = N * ( N - 1 ) / 2
max_calc = int((frames*(frames-1))/2)
#number of calcs per job
calcs_per_job = int(max_calc/nproc)
jobs = list(range(0,max_calc,calcs_per_job))

for i,calculation in enumerate(jobs):
    if i+1 > len(jobs):
        print(i,calculation,max_calc) 
    else:
        print(i,calculation,jobs[i+1])

#farm out RMSD calculations
#concatenate RMSD logs
#cluster
#pull out exemplars
#calcuate density in clusters by trajectory/replica
