map_plumed_dat_file = "AC-SYMM-GM-PLUMED.dat"
n_proc = 48

import os

from subprocess import check_output
out = check_output(["sed '3q;d' conf_box_oriented.gro"],shell=True)
ref = ",".join(out.split()[-3:])

from multiprocessing import Pool

def create_topols(i):
    os.system("gmx_mpi grompp -f nvt_2016.mdp -c conf_{i}.gro -o topol{i}.tpr".format(i=i))

p = Pool(n_proc)
print(p.map(create_topols,range(0,8,1)))

fout = open("plumed.dat","w")
fout.write("""
# RESTART
# include topology info: this is needed to identify atom types
MOLINFO STRUCTURE=structure.pdb

# define all heavy atoms using GROMACS index file
# which can be created with gmx_mpi make_ndx
protein-h: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
protein: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein
# water: GROUP NDX_FILE=index.ndx NDX_GROUP=Water
# make protein whole: add reference position of first heavy atom (in nm)
WHOLEMOLECULES ADDREFERENCE ENTITY0=protein REF0={ref}

# create EMMI score
EMMI ...
LABEL=gmm NOPBC TEMP=300.0 NL_STRIDE=50 NL_CUTOFF=0.01
ATOMS=protein-h GMM_FILE={map_dat}
# resolution is not used - no Bfactor fitting
# SIGMA_MIN can go down to 0.02, must be > NL_CUTOFF,
# if crashes, raise back to 0.05,
# the larger the number the softer the contribution of the EM data
SIGMA_MIN=0.05 RESOLUTION=0.1 NOISETYPE=MARGINAL
# this is needed only with GAUSS noise model
#SIGMA0=0.5 DSIGMA=0.01 WRITE_STRIDE=1000
# this is used to determine scaling factor between predicted and exp maps
#REGRESSION=1000 REG_SCALE_MIN=0.4 REG_SCALE_MAX=10.0 REG_DSCALE=0.001
#SCALE=0.955443
#If there are missing atoms, start a simulation with regression on, get the scale and hard code it by turning on the scale line
...

# translate into bias
emr: BIASVALUE ARG=gmm.scoreb STRIDE=2

# print useful info to file
PRINT ARG=gmm.* FILE=COLVAR STRIDE=500
""".format(ref=ref,map_dat=map_plumed_dat_file))
fout.close()

print("""
You should be set to run:

mpirun -n 24 gmx_mpi mdrun -plumed -multi 8 -ntomp 1

and can analyze on the fly or at end with generate_trajectories_as_pdbs.py
""")
