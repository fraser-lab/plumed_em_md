import os

os.system("gmx_mpi trjcat -f traj_comp?.xtc -o traj_all.xtc -cat")
os.system("gmx_mpi trjconv -f traj_all.xtc -o AC_SYMM_reconstruct_all.xtc -pbc whole -s topol0.tpr")
