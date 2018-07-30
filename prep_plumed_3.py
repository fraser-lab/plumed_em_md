import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--map_dat", help="full path to the GMM dat file for the map",required=True)
parser.add_argument("--n_proc", help="number of threads to use",required=True)
args = parser.parse_args()
map_plumed_dat_file = args.map_dat
n_proc = int(args.n_proc)

import os

try:
    from subprocess import STDOUT, check_output, CalledProcessError
except ImportError:  # pragma: no cover
    # python 2.6 doesn't include check_output
    # monkey patch it in!
    # from: https://stackoverflow.com/questions/4814970/subprocess-check-output-doesnt-seem-to-exist-python-2-6-5
    import subprocess
    STDOUT = subprocess.STDOUT

    def check_output(*popenargs, **kwargs):
        if 'stdout' in kwargs:  # pragma: no cover
            raise ValueError('stdout argument not allowed, '
                             'it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE,
                                   *popenargs, **kwargs)
        output, _ = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd,
                                                output=output)
        return output
    subprocess.check_output = check_output

    # overwrite CalledProcessError due to `output`
    # keyword not being available (in 2.6)
    class CalledProcessError(Exception):

        def __init__(self, returncode, cmd, output=None):
            self.returncode = returncode
            self.cmd = cmd
            self.output = output

        def __str__(self):
            return "Command '%s' returned non-zero exit status %d" % (
                self.cmd, self.returncode)
    subprocess.CalledProcessError = CalledProcessError

out = check_output(["sed '3q;d' conf_box_oriented.gro"],shell=True)
ref = ",".join(out.split()[-3:])

#create protein-no-negative group in index.ndx
conf_in = open("conf_box_oriented.gro")
negative_charges = ["OD1","OD2","OE1","OE2"]
negative_residues = ["ASP","GLU"]
negative_atoms = []
for line in conf_in:
    try:
        if line.split()[1] in negative_charges and line.split()[0][-3:] in negative_residues:
            negative_atoms.append(line.split()[2])
    except:
        pass

os.system("cp index.ndx index.ndx_BACKUP")
index_file = open("index.ndx","r+")
in_protein = 0
protein_atoms = []
for line in index_file:
    if "[ C-alpha ]" in line:
        in_protein = 0
    if in_protein:
        protein_atoms.extend(line.split())
    if "[ Protein-H ]" in line:
        in_protein = 1

for negative_atom in negative_atoms:
    protein_atoms.remove(negative_atom)

index_file.write("[ Protein-no-negative-no-h ]\n")
i = 1
newline = ""
for atom in protein_atoms:
    newline = newline + atom + " "
    i +=1
    if i % 15 == 0:
        index_file.write(newline)
        index_file.write("\n")
        newline = ""
index_file.write(newline)
index_file.close()

from multiprocessing import Pool

def create_topols(i):
    os.system("gmx_mpi grompp -f nvt_2016.mdp -c conf_{i}.gro -o topol{i}.tpr".format(i=i))

p = Pool(n_proc)
print(p.map(create_topols,range(0,8,1)))

fout = open("plumed.dat","w")
fout.write("""# RESTART
# include topology info: this is needed to identify atom types
MOLINFO STRUCTURE=structure.pdb

# define all heavy atoms using GROMACS index file
# which can be created with gmx_mpi make_ndx
protein-h: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
protein: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein
protein-no-negative: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-no-negative-no-h
# water: GROUP NDX_FILE=index.ndx NDX_GROUP=Water
# make protein whole: add reference position of first heavy atom (in nm)
WHOLEMOLECULES ADDREFERENCE ENTITY0=protein REF0={ref}

# create EMMI score
EMMI ...
LABEL=gmm NOPBC TEMP=300.0 NL_STRIDE=50 NL_CUTOFF=0.01
ATOMS=protein-no-negative GMM_FILE={map_dat}
#ATOMS selects what is refined
# no-negative-no-h or no-h are reasonable
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

mpirun -n 8 gmx_mpi mdrun -plumed -multi 8 -ntomp 1

if you need to restart:
Edit plumed.dat
Uncomment out #RESTART
Then on command line try:
mpirun -n 8 gmx_mpi mdrun -plumed -multi 8 -ntomp 1 -cpi state


you can analyze on the fly or at end with:
generate_trajectories_as_pdbs.py
""")
