import os
from subprocess import check_output
out = check_output(["sed '3q;d' conf_box_oriented.gro"],shell=True)
ref = ",".join(out.split()[-3:])

os.system("echo \"q\" | gmx_mpi make_ndx -f topol.tpr")

for i in range(0,8,1):
    fout = open("reconstruct_{i}.dat".format(i=i),"w")
    fout.write("""
    # include topology info: this is needed to identify atom types
    MOLINFO STRUCTURE=structure.pdb

    # define all heavy atoms using GROMACS index file
    # which can be created with gmx_mpi make_ndx
    protein-h: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein-H
    protein: GROUP NDX_FILE=index.ndx NDX_GROUP=Protein
    # make protein whole: add reference position of first heavy atom (in nm)
    WHOLEMOLECULES ADDREFERENCE ENTITY0=protein REF0={ref}

    DUMPATOMS STRIDE=1 FILE=conf_{i}_recon.gro ATOMS=protein
    """.format(ref=ref,i=i))
    fout.close()

    os.system("plumed driver --plumed reconstruct_{i}.dat --igro conf_{i}.gro".format(i=i))
    os.system("echo 0 | gmx_mpi trjconv -f conf_{i}_recon.gro -s conf_{i}_recon.gro -o structure{i}.pdb".format(i=i))

print("""

    CHECK that structure{0-7}.pdb fits in your original map!

    when ready proced to prep_plumed3.py
""")
