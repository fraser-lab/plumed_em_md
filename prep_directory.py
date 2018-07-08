import os
plumed_working_dir = "/home/jfraser/plumed_em_md"

os.system("cp {dir}/*.mdp .".format(dir=plumed_working_dir))
