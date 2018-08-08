import MDAnalysis
import MDAnalysis.analysis.encore as encore

reference_pdb = "structure.pdb"
xtc1 = "AC_NCS_0.xtc"
xtc2 = "AC_NCS_1.xtc"

ens1 = MDAnalysis.Universe(reference_pdb)
ens2 = MDAnalysis.Universe(reference_pdb)
ens1.load_new(xtc1)
ens2.load_new(xtc2)
    # loop_ca = u.select_atoms('resid 34-50 and name CA')

coordinates, details = encore.reduce_dimensionality([ens1,ens2],selection='resid 34-50 and name CA',ncores=48)
plt.scatter(coordinates[0], coordinates[1],
            color=[["red", "blue"][m-1] for m
            in details["ensemble_membership"]])
plt.savefig("dim_red.png")
