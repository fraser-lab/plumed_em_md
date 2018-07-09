# plumed_em_md

requires:
- PLUMED (https://github.com/plumed/plumed2.git)
- GROMACS 2016.5 (http://ftp.gromacs.org/pub/gromacs/gromacs-2016.5.tar.gz)
- gmmconvert (from Max)
- gmmconvert (from sbgrid)


Scripts to make running a PLUMED EM MD simulation (https://doi.org/10.1016/j.bpj.2018.02.028) much easier.

Copy into working directory:

Input pdb
Input Map
em_2016.mdp
npt_2016.mdp
nvt_2016_EQUIL.mdp
nvt_2016.mdp

Prepare map:
Generate_gmm.py (set input map, gmconvert locations)
convert_GMM2PLUMED.sh (generates input gmm.dat file)

Prepare md simulation:
Prep_plumed.py (set input pdb)
Prep_plumed_2.py
Prep_plumed_3.py (set input gmm.dat file)

Analyze md simulation:
generate_trajectories_as_pdbs.py
