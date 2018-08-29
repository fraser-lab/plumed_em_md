xtcs = ["AC_NCS_med_0.xtc","AC_NCS_med_1.xtc","AC_NCS_med_2.xtc","AC_NCS_med_3.xtc","AC_NCS_med_4.xtc","AC_NCS_med_5.xtc","AC_NCS_med_6.xtc","AC_NCS_med_7.xtc","dAC_NCS_med_0.xtc","dAC_NCS_med_1.xtc","dAC_NCS_med_2.xtc","dAC_NCS_med_3.xtc","dAC_NCS_med_4.xtc","dAC_NCS_med_5.xtc","dAC_NCS_med_6.xtc","dAC_NCS_med_7.xtc"]

pdbs = ["AC_NCS_structure.pdb","dAC_NCS_structure.pdb"]

import numpy
# cluster_count = {}
cluster_counter = numpy.zeros((16,26))
ac_dc_counter = numpy.zeros((26,2))


f = open("trajectory.dat")
for line in f:
    frame = int(line.split()[0])
    cluster = int(line.split()[1])
    xtc_index = divmod(frame,1201)[0]
    frame_index = divmod(frame,1201)[1]
    cluster_counter[xtc_index][cluster] +=1
    if "dAC" in xtcs[xtc_index]:
        ac_dc_counter[cluster][1] +=1
    else:
        ac_dc_counter[cluster][0] +=1
print(cluster_counter)
print(ac_dc_counter)


import matplotlib
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
im = ax.imshow(ac_dc_counter)
ax.set_xticks(numpy.arange(2))
ax.set_xticklabels(("A","D"))

for i in range(26):
    for j in range(2):
        print(ac_dc_counter[i,j], end='\t')
    print("")
        # text = ax.text(j, i, ac_dc_counter[i, j],
        #                ha="center", va="center", color="w")

plt.savefig("ac_dc.png")
