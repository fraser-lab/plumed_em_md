from subprocess import check_output
import argparse
import os
from multiprocessing import Pool

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

if not os.path.isfile("output.log"): #so that we can skip redoing the RMSD calc
    #divide into number of processors
    #max number of calculations k1 = N * ( N - 1 ) / 2
    max_calc = int((frames*(frames-1))/2)
    #number of calcs per job
    calcs_per_job = int(max_calc/nproc)
    jobs = list(range(0,max_calc,calcs_per_job))
    print(jobs, max_calc, calcs_per_job,len(jobs))

    job_assign = []

    for i,calculation in enumerate(jobs):
        if i+1.1 > len(jobs):
            print(i,calculation,max_calc)
        else:
            print(i,calculation,jobs[i+1])


    #farm out RMSD calculations

    def calculate_rmsds(i):
        # os.system("python traj_rmsd.py structure.pdb reconstruct_small.xtc {i}.log 0 528 ")
        if i+1.1 > len(jobs):
            j = max_calc
        else:
            j = jobs[i+1]
        os.system("python ~/plumed_em_md/ANALYSIS/traj_rmsd.py structure.pdb {xtc} {i}_PRECLUST.log {start} {stop}".format(xtc=xtc,i=i,start=jobs[i],stop=j))
        return

    p = Pool(nproc)
    print(p.map(calculate_rmsds,range(0,len(jobs),1)))

    #concatenate RMSD logs

    os.system("cat *_PRECLUST.log > output.log")
else:
    pass

#cluster
rms = 2.0

os.system("~/max_cluster.exe output.log 1 {rms} {frames}".format(frames=frames,rms=rms))

#pull out exemplars
#calcuate density in clusters by trajectory/replica
#Select just backbones
#align them together based on NOT loop
#then do analysis
