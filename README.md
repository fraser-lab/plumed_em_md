# plumed_em_md

requires:
- PLUMED (https://github.com/plumed/plumed2.git)
- GROMACS 2016.5 (http://ftp.gromacs.org/pub/gromacs/gromacs-2016.5.tar.gz)
- gmmconvert (from Max)
- gmmconvert (from sbgrid)
- CHARMM forcefield (for ALY)


Scripts to make running a PLUMED EM MD simulation (https://doi.org/10.1016/j.bpj.2018.02.028) much easier.

Given:
- Input PDB (cleaned of heteroatoms)
- Input Map (MRC format)

prep_directory.py - copies instruction files into working directory:
- em_2016.mdp
- npt_2016.mdp
- nvt_2016_EQUIL.mdp
- nvt_2016.mdp

generate_gmm.py - prepares map in gmm format
convert_GMM2PLUMED.sh - generates input gmm.dat file for plumed (from Max)

Prep_plumed.py - equilibrates
Prep_plumed_2.py - sets up simulation
Prep_plumed_3.py - outputs command to run simulation

generate_trajectories_as_pdbs.py - outputs md simulation as gro or pdb files
