from multiprocessing import Pool
import os

input_map = "/home/jfraser/working_with_max/07012018_AC_SYMM/AC_SYMM_MASKED_10_ORIGIN_0.mrc"
n_proc = 48
gmconvert_location = "/home/jfraser/gmconvert/gmconvert"
sbgrid_gmconvert = "/programs/x86_64-linux/system/sbgrid_bin/gmconvert"

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
    os.system("{gmconvert} -imap {map} -gmml GMM.list -ng 100 -ogmm {i}_0.gmm -zth 0.0".format(map=input_map,i=i,gmconvert=gmconvert_location))
    os.chdir("../")
    return

p = Pool(n_proc)
print(p.map(create_gmms,range(0,150,1)))

os.system("cat map_*/*.gmm > 1.gmm")
os.system("{gmconvert} VcmpG -igmm 1.gmm -imap {map} -zth 0.0 -omap iter_1.mrc > iter_1.log".format(map=input_map, gmconvert=sbgrid_gmconvert))
os.system("grep CCgrid iter_1.log")
