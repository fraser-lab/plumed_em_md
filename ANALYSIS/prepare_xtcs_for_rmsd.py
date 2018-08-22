import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--prefix", help="name for xtc",required=True)
parser.add_argument("--replicas", help="nubmer of replicas",type=int, required=True)
parser.add_argument("--skip", help="number of frames to skip in writing xtc",type=int, required=True, default=1)
parser.add_argument("--end", help="end time",type=int, default=12000)
args = parser.parse_args()
replicas = args.replicas
prefix = args.prefix
skip = args.skip
end = args.end

import os

#create individual xtcs (for rmsf analysis)
for replica in range(0,replicas):
    print(replica)
    os.system("printf \"14\" | gmx_mpi trjconv -f traj_comp{replica}.xtc -o {prefix}_{replica}.xtc -pbc whole -s topol0.tpr -skip {skip} -e {end}".format(prefix=prefix,replica=replica,skip=skip, end=end))

# #create concatenated xtc
os.system("gmx_mpi trjcat -f {prefix}_?.xtc -o traj_all.xtc -cat".format(prefix=prefix))
os.system("printf \"14\" | gmx_mpi trjconv -f traj_all.xtc -o {prefix}.xtc -pbc whole -s topol0.tpr".format(prefix=prefix))
