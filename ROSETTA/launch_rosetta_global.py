import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--pdb", help="pdb file name output from idealize",required=True)
parser.add_argument("--map", help="map file name ",required=True)
parser.add_argument("--rms", help="rms to sample ",default=2.0)
parser.add_argument("--resolution", help="reported resolution ",required=True, type=float)
parser.add_argument("--tasks", help="number of jobs to submit",required=True,type=int)
parser.add_argument("--nstruct", help="number of structures to generate per job",required=True,type=int)

args = parser.parse_args()
pdb_ideal = args.pdb
map = args.map
rms = args.rms
resolution = args.resolution
nstruct = args.nstruct
tasks = args.tasks

### FROM GABE LANDER:
# determine electron density weight
# according to Rosetta:
# this value should be 10 for 3A, 20-25 for 3.5-4A, and 50-60 for 4.5-5A
dens = 5
if resolution > 2.8:
    dens = int((resolution-2.58)/0.0411)
print "\n\n\ndensity weight set to {dens}\n\n\n".format(dens=dens)

fout = open("launch_rosetta_refine.sh","w")

fout.write("""#!/bin/bash
#
#$ -S /bin/bash
#$ -o ./out
#$ -e ./err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=1G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=80:00:00
#$ -t 1-{tasks}

hostname
date

source /programs/sbgrid.shrc
rosetta_scripts.linuxgccrelease -database /netapp/home/jaimefraser/database -in::file::s {pdb} -edensity::mapfile {map} -parser::protocol new_multi_global.xml   -edensity::mapreso {resolution} -default_max_cycles 200 -edensity::cryoem_scatterers -out::suffix $SGE_TASK_ID -crystal_refine -beta -parser::script_vars denswt={dens} rms={rms} -nstruct {nstruct}

date
""".format(pdb=pdb_ideal,map=map,dens=dens,rms=rms,resolution=resolution,tasks=tasks,nstruct=nstruct))

fout.close()

fxml = open("new_multi_global.xml","w")
fxml.write("""<ROSETTASCRIPTS>
   <SCOREFXNS>
      <ScoreFunction name="cen" weights="score4_smooth_cart">
         <Reweight scoretype="elec_dens_fast" weight="20" />
      </ScoreFunction>
      <ScoreFunction name="dens_soft" weights="beta_soft">
         <Reweight scoretype="cart_bonded" weight="0.5" />
         <Reweight scoretype="pro_close" weight="0.0" />
         <Reweight scoretype="elec_dens_fast" weight="%%denswt%%" />
      </ScoreFunction>
      <ScoreFunction name="dens" weights="beta_cart">
         <Reweight scoretype="elec_dens_fast" weight="%%denswt%%" />
         <Set scale_sc_dens_byres="R:0.76,K:0.76,E:0.76,D:0.76,M:0.76, C:0.81,Q:0.81,H:0.81,N:0.81,T:0.81,S:0.81,Y:0.88,W:0.88, A:0.88,F:0.88,P:0.88,I:0.88,L:0.88,V:0.88" />
      </ScoreFunction>
   </SCOREFXNS>
   <MOVERS>
      <SetupForDensityScoring name="setupdens" />
      <SwitchResidueTypeSetMover name="tocen" set="centroid" />
      <MinMover name="cenmin" scorefxn="cen" type="lbfgs_armijo_nonmonotone" max_iter="200" tolerance="0.00001" bb="1" chi="1" jump="ALL" />
      <CartesianSampler name="cen5_50" automode_scorecut="-0.5" scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" strategy="auto" fragbias="density" rms="%%rms%%" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4" fraglens="7" nfrags="25" />
      <CartesianSampler name="cen5_60" automode_scorecut="-0.3" scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" strategy="auto" fragbias="density" rms="%%rms%%" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4" fraglens="7" nfrags="25" />
      <CartesianSampler name="cen5_70" automode_scorecut="-0.1" scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" strategy="auto" fragbias="density" rms="%%rms%%" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4" fraglens="7" nfrags="25" />
      <CartesianSampler name="cen5_80" automode_scorecut="0.0" scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" strategy="auto" fragbias="density" rms="%%rms%%" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4" fraglens="7" nfrags="25" />
      <BfactorFitting name="fit_bs" max_iter="50" wt_adp="0.0005" init="1" exact="1" />
      <FastRelax name="relaxcart" scorefxn="dens" repeats="1" cartesian="1" />
   </MOVERS>
   <PROTOCOLS>
      <Add mover="setupdens" />
      <Add mover="tocen" />
      <Add mover="cenmin" />
      <Add mover="relaxcart" />
      <Add mover="cen5_50" />
      <Add mover="relaxcart" />
      <Add mover="cen5_60" />
      <Add mover="relaxcart" />
      <Add mover="cen5_70" />
      <Add mover="relaxcart" />
      <Add mover="cen5_80" />
      <Add mover="relaxcart" />
      <Add mover="relaxcart" />
   </PROTOCOLS>
   <OUTPUT scorefxn="dens" />
</ROSETTASCRIPTS>
""")
fxml.close()

import os
os.system("qsub launch_rosetta_refine.sh")
