from multiprocessing import Pool
import os
trajectories = range(0,8)

from subprocess import check_output
out = check_output(["sed '3q;d' conf_box_oriented.gro"],shell=True)
ref = ",".join(out.split()[-3:])

def convert_to_gro(traj):
    print traj
    fout = open("scripted_fix_md_output%i.dat" %(traj),"w")
    fout.write(
"""
# include topology info: this is needed to identify atom types
MOLINFO STRUCTURE=structure.pdb

# define all heavy atoms using GROMACS index file
# which can be created with gmx_mpi make_ndx
protein-h: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
protein: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein
water: GROUP NDX_FILE=index.ndx NDX_GROUP=Water

# make protein whole: add reference position of first heavy atom (in nm)
WHOLEMOLECULES ADDREFERENCE ENTITY0=protein  REF0=%s

DUMPATOMS STRIDE=10 FILE=traj%i.gro ATOMS=protein
""" %(ref,traj))
    fout.close()
    os.system("plumed driver --plumed scripted_fix_md_output{traj}.dat --mf_trr traj{traj}.trr".format(traj=traj))
    os.system("echo 1 | gmx_mpi trjconv -f traj{traj}.gro -o traj{traj}.pdb".format(traj=traj))

# TO analyze water specifically FOLLOW RECIPE on page 18-19-20 of plumed_chapter.pdf

p = Pool(8)
print(p.map(convert_to_gro,range(0,1,1)))
