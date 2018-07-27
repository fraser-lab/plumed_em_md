Clustering proceeds in steps:

1) you concatenate all the trajectories of all the replicas, using gromacs trjcat utility:

gmx_mpi trjcat -f traj*.xtc -o traj_all.xtc -cat 

It is better to work with the xtc files, which contain no water, just the protein.
If traj_all.xtc contains more than, let’s say, 30000 frames, you can write one every N frames:

gmx_mpi trjconv -f traj_all.xtc -o traj_all_small.xtc -skip N

1) then you need to calculate the RMSD between all pairs of frames in your trajectory.
    You need to modify the traj_rmsd.py python script at this line:

    allheavy = u.select_atoms('type C N O S’)

    to make a selection of the atoms you want for the RMSD calculation.
    The library (and syntax) used is MDAnalysis, which of course you need installed.

    You can calculate the RMSD matrix using multiple cores.
    Usually I do this on the cluster using an array job (see go_rmsd.sh).
    I think this file is sufficiently commented, so you can understand which parameters you need
    to change (total number of RMSD calculations, number of calculation per task, input files - pdb and trajectory).

2) once you have M RMSD sub-matrices (one per core used), you can do the actual clustering.
    You need to compile gromos_clustering.cpp, which implements the clustering method of GROMACS (gmx_mpi cluster), in a more efficient way.
    This code is serial, and an example of command line is in go_cluster.pbs. 

3) The output of gromos_clustering are two files: log.dat and trajectory.dat.
     log.dat contains a list of clusters, their population, and the frame number (zero based) of the cluster center
    (a representative structure of the cluster)
    trajectory.dat contains for each frame of the trajectory the assignment to a given cluster.

4) If you want to assess convergence of the simulation, you can post-process trajectory.dat.
    You need to remember the frame id of all the frames belonging to the first and second half of the simulation.
    At this point, you can recalculate the clusters populations using only the first or last half of the simulation.
    At convergence, they should be consistent within a few percent (even 10% is ok).
