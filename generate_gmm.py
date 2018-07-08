from multiprocessing import Pool
import os
gmconvert_location = "/home/jfraser/gmconvert/gmconvert"
sbgrid_gmconvert = "/programs/x86_64-linux/gmconvert/20180516/bin/gmconvert"

import argparse, sys
parser = argparse.ArgumentParser()
parser.add_argument("--input_map", help="full path to the input map for simulation",required=True)
parser.add_argument("--n_proc", help="number of threads to use",required=True)
args = parser.parse_args()

n_proc = int(args.n_proc)
input_map = args.input_map

os.system("{gmconvert} -imap {map} -oimap output-0.0.mrc -zth 0.0".format(map=input_map, gmconvert=gmconvert_location))
os.mkdir("ITER_0")
os.mkdir("ITER_1")
os.chdir("ITER_0")
os.system("{gmconvert} V2G -imap ../output-0.0.mrc -ogmm 0.gmm -ng 150 -zth 0.0".format(gmconvert=gmconvert_location))
os.chdir("../")

os.chdir("ITER_1")
for i in range(0,150,1):
    os.mkdir("map_0_{i}".format(i=i))
    os.chdir("map_0_{i}".format(i=i))
    os.system("echo \"../../ITER_0/0.gmm {i}\" > GMM.list".format(i=i))
    os.chdir("../")

def create_gmms(i):
    os.chdir("map_0_{i}".format(i=i))
    os.system("{gmconvert} -imap ../../{map} -gmml GMM.list -ng 100 -ogmm {i}_0.gmm -zth 0.0".format(map=input_map,i=i,gmconvert=gmconvert_location))
    os.chdir("../")
    return

p = Pool(n_proc)
print(p.map(create_gmms,range(0,150,1)))

os.system("cat map_*/*.gmm > 1.gmm")
os.system("{gmconvert} VcmpG -igmm 1.gmm -imap {map} -zth 0.0 -omap iter_1.mrc > iter_1.log".format(map=input_map, gmconvert=sbgrid_gmconvert))
os.system("grep CCgrid iter_1.log")
