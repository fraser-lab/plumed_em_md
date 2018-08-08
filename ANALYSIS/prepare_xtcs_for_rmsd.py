import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--prefix", help="name for xtc",required=True)
args = parser.parse_args()
prefix = args.prefix


import os

os.system("gmx_mpi trjcat -f traj_comp?.xtc -o traj_all.xtc -cat")
os.system("printf \"14\" | gmx_mpi trjconv -f traj_all.xtc -o {prefix}.xtc -pbc whole -s topol0.tpr".format(prefix=prefix))
