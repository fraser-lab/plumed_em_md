#!/bin/bash
#PBS -N rmsd_matrix 
#PBS -l nodes=1:ppn=1
#PBS -l walltime=48:00:00
#PBS -q l1
#PBS -o pbs.out
#PBS -e pbs.err
#PBS -m a
#PBS -V
#PBS -t 0-63

# go to right directory
cd $PBS_O_WORKDIR

# total number of entries in rmsd matrix
# if N if the number of frames in the trajectory,
# then nentries = N * ( N - 1 ) / 2
nentries=285234670
# number of entries per task (rounded up)
ntask=4456792
# id of the job (zero based!)
id=$PBS_ARRAYID
# range of entries for rmsd calculation
i0=$(( $id * $ntask ))
i1=$(( $i0 + $ntask ))
# check if over limit, and reset
if [ $i1 -gt $nentries ]; then i1=$nentries; fi

## DATA dir
DIR=/home/mb2006/SCRATCH/ClpP/CLUSTERING/DATA_FILES

# go
python traj_rmsd.py ${DIR}/ClpP_structure.pdb ${DIR}/traj_10_tstep.xtc RMSD.$id $i0 $i1
