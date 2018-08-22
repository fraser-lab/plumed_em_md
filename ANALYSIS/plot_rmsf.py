import MDAnalysis
import numpy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--pdb", help="starting conformation pdb, usually structure.pdb",required=True)
parser.add_argument("--prefix", help="trajectory xtc prefix",required=True)
parser.add_argument("--replicas", help="nubmer of replicas",type=int, required=True)
parser.add_argument("--output", help="output file prefix",required=True)
args = parser.parse_args()
replicas = args.replicas
reference_pdb = args.pdb
prefix = args.prefix
output = args.output

for replica in range(0,replicas):
    xtc = "{prefix}_{replica}.xtc".format(prefix=prefix,replica=replica)
    print(xtc)
    u = MDAnalysis.Universe(reference_pdb)
    u.load_new(xtc)
    loop_ca = u.select_atoms('resid 34-50 and name CA')

    from MDAnalysis.analysis.rms import RMSF

    rmsfer = RMSF(loop_ca, verbose=True).run()
    print(rmsfer.rmsf)
    import matplotlib.pyplot as plt

    plt.plot(loop_ca.resnums, rmsfer.rmsf)

plt.ylim(0,6)
plt.savefig("{output}_rmsf.png".format(output=output))
