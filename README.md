# plumed_em_md

Scripts to make running a PLUMED EM MD simulation (https://doi.org/10.1016/j.bpj.2018.02.028) much easier.

Requires:
- python 2.6.5 or higher
- PLUMED (https://github.com/plumed/plumed2.git)
- GROMACS 2016.5 (http://ftp.gromacs.org/pub/gromacs/gromacs-2016.5.tar.gz)
- gmmconvert (from Max, in this repo)
- gmmconvert (from sbgrid, for correlations)
- CHARMM forcefield (for ALY from Max, based on http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs)

Install:
- confirm gcc and mpicc match (on qb3 scl enable devtoolset-7 bash, module load openmpi-1.8-x86_64)
- compile plumed
  - git clone https://github.com/plumed/plumed2.git
  - cd plumed2
  - git checkout isdb
  - ./configure --disable-python
  - make -j 12
  compile GROMACS
  - wget http://ftp.gromacs.org/pub/gromacs/gromacs-2016.5.tar.gz
  - tar -xvf gromacs-2016.5.tar.gz
  - cd /netapp/home/jaimefraser/gromacs-2016.5
  - plumed-patch -p --shared
  - Choose gromacs-2016.5 (option 2)
  - mkdir build
  - cd build
  - mkdir /netapp/home/jaimefraser/gromacs-2016.5-bin/
  - cmake ../ -DBUILD_SHARED_LIBS=ON -DGMX_OPENMP=OFF -DGMX_THREAD_MPI=OFF -DGMX_GPU=OFF -DCMAKE_INSTALL_PREFIX=/netapp/home/jaimefraser/gromacs-2016.5-bin/ -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc -DGMX_MPI=ON -DGMX_USE_RDTSCP=off
    - Need the DGMX_USE_RDTSCP flag off b/c iqint is different subtly than some other nodes on qb3 cluster! if using in an environment with only one CPU type, delete that flag.
  - make -j 12
  - make install

Edit bash_profile or launch scripts to include:
- source ~/plumed2/sourceme.sh
- source ~/gromacs-2016.5-bin/bin/GMXRC

Prep maps:
- generate_gmm.py - prepares map in gmm format
- convert_GMM2PLUMED.sh - generates input gmm.dat file for plumed (from Max)

Prep simulation:
- prep_directory.py - copies instruction files into working directory:
- prep_plumed.py - equilibrates
- prep_plumed_2.py - sets up simulation
- prep_plumed_3.py - outputs command to run simulation
  - One significant change relative to published method, we remove negative scatterers from contributing to the density calculation.
- simulation can then be run on command line

On SGE:
- launch_gmm.py - generates GMM
- launch_prep.py - equilibrates
- lauch_prep2.py - sets up simulation
- launch_prep3.py - outputs commands to run simulation
- launch_simulation.py - runs simulation

Analysis:
- still a work in progress
- generate_trajectories_as_pdbs.py - outputs md simulation as gro or pdb files
