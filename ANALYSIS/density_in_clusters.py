xtcs = ["AC_NCS_med_0.xtc","AC_NCS_med_1.xtc","AC_NCS_med_2.xtc","AC_NCS_med_3.xtc","AC_NCS_med_4.xtc","AC_NCS_med_5.xtc","AC_NCS_med_6.xtc","AC_NCS_med_7.xtc","dAC_NCS_med_0.xtc","dAC_NCS_med_1.xtc","dAC_NCS_med_2.xtc","dAC_NCS_med_3.xtc","dAC_NCS_med_4.xtc","dAC_NCS_med_5.xtc","dAC_NCS_med_6.xtc","dAC_NCS_med_7.xtc"]

pdbs = ["AC_NCS_structure.pdb","dAC_NCS_structure.pdb"]

#COUNT 

f = open("trajectory.dat")
for line in f:
    frame = line.split()[0]
    cluster = line.split()[1]
    xtc_index = divmod(frame,1201)[0]
    frame_index = divmod(frame,1201)[1]


# import MDAnalysis
# u = MDAnalysis.Universe("structure_test.pdb")
# u.load_new("both_med_bb.xtc")
#     # protein = u.select_atoms("resid 20-60")
# for ts in u.trajectory:
#     if u.trajectory.time == 0:
#         print(ts.frame, u.trajectory.time)
