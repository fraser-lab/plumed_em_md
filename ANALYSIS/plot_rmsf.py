import MDAnalysis
import numpy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--pdb", help="starting conformation pdb",required=True)
parser.add_argument("--xtc", help="trajectory xtc file",required=True)
args = parser.parse_args()
reference_pdb = args.pdb
xtc = args.xtc

u = MDAnalysis.Universe(reference_pdb)
u.load_new(xtc)
loop_ca = u.select_atoms('resid 36-48 and name CA')

from MDAnalysis.analysis.rms import RMSF

rmsfer = RMSF(loop_ca, verbose=True).run()

import matplotlib.pyplot as plt

plt.plot(calphas.resnums, rmsfer.rmsf)
