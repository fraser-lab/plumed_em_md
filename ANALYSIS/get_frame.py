import MDAnalysis, argparse

parser = argparse.ArgumentParser()
# parser.add_argument("--pdb", help="starting conformation pdb, likely structure.pdb",required=True)
# parser.add_argument("--xtc", help="trajectory xtc file",required=True)
# args = parser.parse_args()
# reference_pdb = args.pdb
# xtc = args.xtc

xtcs = ["AC_NCS_med_0.xtc","AC_NCS_med_1.xtc","AC_NCS_med_2.xtc","AC_NCS_med_3.xtc","AC_NCS_med_4.xtc","AC_NCS_med_5.xtc","AC_NCS_med_6.xtc","AC_NCS_med_7.xtc","dAC_NCS_med_0.xtc","dAC_NCS_med_1.xtc","dAC_NCS_med_2.xtc","dAC_NCS_med_3.xtc","dAC_NCS_med_4.xtc","dAC_NCS_med_5.xtc","dAC_NCS_med_6.xtc","dAC_NCS_med_7.xtc"]

pdbs = ["AC_NCS_structure.pdb","dAC_NCS_structure.pdb"]

centroids = [9642,1202,11850,3642,4903,15720,9557,17417,10478,35,2410,3964]

for i,centroid in enumerate(centroids):
    xtc_index = divmod(centroid,1201)[0]
    frame_index = divmod(centroid,1201)[1]
    if "dAC" in xtcs[xtc_index]:
        u = MDAnalysis.Universe(pdbs[1])
    else:
        u = MDAnalysis.Universe(pdbs[0])
    u.load_new(xtcs[xtc_index])
    protein = u.select_atoms("resid 20-60")
    for ts in u.trajectory:
        # if u.trajectory.time == 0:
        #     print(ts.frame, u.trajectory.time)
        if ts.frame == frame_index:
            with MDAnalysis.Writer("{centroid}_{xtc_index}_{frame}.pdb".format(xtc_index=xtc_index,centroid=i,frame=ts.frame)) as W:
                W.write(protein)
