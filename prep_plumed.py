input_pdb = "/home/jfraser/working_with_max/07012018_AC_SYMM/AC_SYMM.pdb"

import os
os.system("printf \"9\n1\" | gmx_mpi pdb2gmx -f {pdb} -ignh".format(pdb=input_pdb))
os.system("gmx_mpi editconf -f conf.gro -bt dodecahedron -d 1.0 -o conf_box.gro")
os.system("printf \"0\n0\" | gmx_mpi trjconv -f conf_box.gro -s conf.gro -fit translation -o conf_box_oriented.gro")
os.system("gmx_mpi solvate -cp conf_box_oriented.gro -cs spc216.gro -p topol.top -o conf_water.gro")
os.system("gmx_mpi grompp -f em_2016.mdp -c conf_water.gro")
os.system("gmx_mpi mdrun -c conf_emin.gro")
os.system("gmx_mpi grompp -f npt_2016.mdp -c conf_emin.gro")

#this one takes some time
os.system("gmx_mpi mdrun -c conf_npt.gro -nsteps 250000 -ntomp 16")
#this one takes some time

os.system("gmx_mpi grompp -f nvt_2016_EQUIL.mdp -c conf_npt.gro")

#consider editing posre.itp to less than 1000 for the loop residues
# 5 or 10 likely generates a displacement of >1A
# 200 likely generates displacement of ~1A
# 2/(displacement in nm^2) = number to input
# this will generate greater diversity in the nvt simulation

#this one takes some time
os.system("gmx_mpi mdrun -c conf_nvt.gro -nsteps 1000000 -ntomp 16")
#this one takes some time

for i in range(0,8,1):
    os.system("echo 0 | gmx_mpi trjconv -f traj.trr -s topol.tpr -dump {f} -o conf_{i}.gro".format(i=i,f=i*140))

os.system("echo 0 | gmx_mpi trjconv -f conf_box_oriented.gro -s conf_box_oriented.gro -o structure.pdb")

print("""

CHECK that structure.pdb fits in your original map!

when ready proced to prep_plumed2.py
""")
