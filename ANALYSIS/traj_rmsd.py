import MDAnalysis
import MDAnalysis.analysis.rms
import sys
import math

# initialize MDAnalysis
# load PDB
u = MDAnalysis.Universe(sys.argv[1])
# load XTC/DCD/TRR
u.load_new(sys.argv[2])
# all-heavy selector - or any other selection
allheavy = u.select_atoms('resid 36-48 and name CA')
# number of frames
n = len(u.trajectory)

# open output log file
log=open(sys.argv[3], "w")

# first and last element of the 1D vector created below
# for everything include 0 as first
# and
# if N if the number of frames in the trajectory,
# then k1 = N * ( N - 1 ) / 2
# because last number will be excluded and python counts from 0

k0 = int(sys.argv[4])
k1 = int(sys.argv[5])

# cycle on pairs and write on a separate file
for k in range(k0, k1):
    # map one-dimensional index to two-dimensional indexes
    i = int(n - 2 - math.floor(math.sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5))
    j = int(k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2)
    # go to frame i
    u.trajectory[i]
    # copy positions of all-heavy atoms
    A = allheavy.positions.copy()
    # go to frame j
    u.trajectory[j]
    # copy positions of all-heavy atoms
    B = allheavy.positions.copy()
    # calculate RMSD with optimal alignment - put False to deactivate alignment
    rmsd = MDAnalysis.analysis.rms.rmsd(A, B, superposition=True)
    # print on file
    log.write("%6d %6d %6.3f\n" % (i, j, rmsd))

# close file
log.close()
