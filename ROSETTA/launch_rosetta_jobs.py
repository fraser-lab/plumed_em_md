import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--pdb", help="pdb file name output from idealize",required=True)
parser.add_argument("--map", help="map file name ",required=True)
args = parser.parse_args()
pdb_ideal = args.pdb
map = args.map

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

hostname
date

source /programs/sbgrid.shrc
rosetta_scripts.linuxgccrelease -database /netapp/home/jaimefraser/database -in::file::s {pdb} -edensity::mapfile {map} -parser::protocol new_multi_ray.xml   -edensity::mapreso 3.5 -default_max_cycles 200 -edensity::cryoem_scatterers -out::suffix _asymm -crystal_refine -beta

date
""".format(pdb=pdb_ideal,map=map)) #SUFFIX can be $SGE_TASK_ID
fout.close()

fxml = open("new_multi_ray.xml","w")
fxml.write("""<ROSETTASCRIPTS>
     <SCOREFXNS>
         <ScoreFunction name="cen" weights="score4_smooth_cart">
             <Reweight scoretype="elec_dens_fast" weight="15"/>
         </ScoreFunction>
         <ScoreFunction name="dens_soft" weights="beta_soft">
             <Reweight scoretype="cart_bonded" weight="0.5"/>
             <Reweight scoretype="pro_close" weight="0.0"/>
             <Reweight scoretype="elec_dens_fast" weight="25"/>
         </ScoreFunction>
             <ScoreFunction name="dens" weights="beta_cart">
                 <Reweight scoretype="elec_dens_fast" weight="20"/>
                 <Set scale_sc_dens_byres="R:0.76,K:0.76,E:0.76,D:0.76,M:0.76, C:0.81,Q:0.81,H:0.81,N:0.81,T:0.81,S:0.81,Y:0.88,W:0.88, A:0.88,F:0.88,P:0.88,I:0.88,L:0.88,V:0.88"/>
               </ScoreFunction>
             <ScoreFunction name="dens_tor" weights="beta">
                 <Reweight scoretype="elec_dens_fast" weight="20"/>
                 <Set scale_sc_dens_byres="R:0.76,K:0.76,E:0.76,D:0.76,M:0.76, C:0.81,Q:0.81,H:0.81,N:0.81,T:0.81,S:0.81,Y:0.88,W:0.88, A:0.88,F:0.88,P:0.88,I:0.88,L:0.88,V:0.88"/>
         </ScoreFunction>
     </SCOREFXNS>
     <MOVERS>
         <SetupForDensityScoring name="setupdens"/>
         <CartesianSampler name="longfrag" strategy="user" residues="21A-63A" rsd_wdw_to_refine="5"
             scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" fragbias="density"
             rms="4" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4"
             fraglens="11" nfrags="25"/>

         <CartesianSampler name="shortfrag" strategy="user" residues="21A-63A" rsd_wdw_to_refine="5"
             scorefxn="cen" mcscorefxn="cen" fascorefxn="dens_soft" fragbias="density"
             rms="4" ncycles="200" fullatom="0" bbmove="1" nminsteps="25" temp="4"
             fraglens="5" nfrags="25"/>


         <BfactorFitting name="fit_bs" max_iter="50" wt_adp="0.0005" init="1" exact="1"/>
         <FastRelax name="relaxcart" scorefxn="dens" repeats="1" cartesian="1"/>
         <FastRelax name="relaxtors" scorefxn="dens_tor" repeats="3" />
     </MOVERS>
     <PROTOCOLS>
         <Add mover="setupdens"/>
         <Add mover="longfrag"/>
         <Add mover="shortfrag"/>
         <Add mover="relaxtors"/>
         <Add mover="relaxcart"/>
     </PROTOCOLS>
     <OUTPUT scorefxn="dens"/>
</ROSETTASCRIPTS>
""")
fxml.close()

import os
os.system("qsub launch_rosetta_refine.sh")
