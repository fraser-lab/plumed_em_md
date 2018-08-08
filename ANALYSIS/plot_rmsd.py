import MDAnalysis
import numpy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--pdb", help="starting conformation pdb, likely structure.pdb",required=True)
parser.add_argument("--xtc", help="trajectory xtc file",required=True)
args = parser.parse_args()
reference_pdb = args.pdb
xtc = args.xtc

u = MDAnalysis.Universe(reference_pdb)
u.load_new(xtc)
loop_ca = u.select_atoms('resid 36-48 and name CA')

# u.trajectory[0]   # rewind trajectory
print(len(u.trajectory))

import MDAnalysis.analysis.rms
R = MDAnalysis.analysis.rms.RMSD(u,select='resid 36-48 and name CA')
R.run()
rmsd = R.rmsd.T

replicas = {}
i = 0
last_time = 1
for j,time in enumerate(rmsd[1]):
    if last_time < time:
        # print(last_time,time,rmsd[2][j],i)
        if rmsd[2][j] < 10:
            # print (last_time,time,rmsd[2][j],i)
            replicas[i]["rmsds"].append(rmsd[2][j])
            replicas[i]["times"].append(rmsd[1][j])
        else:
            print (last_time,time,rmsd[2][j],i)
        last_time = time
    else:
        # print("RESET",last_time,time,rmsd[2][j],i)
        i += 1
        replicas[i] = {}
        replicas[i]["rmsds"] = [rmsd[2][j]]
        replicas[i]["times"] = [rmsd[1][j]]
        last_time = time

import matplotlib.pyplot as plt

for series in replicas:
    plt.plot(replicas[series]["times"],replicas[series]["rmsds"])
    print(numpy.mean(replicas[series]["rmsds"][1:]),numpy.std(replicas[series]["rmsds"][1:]))
plt.savefig("test.png")
