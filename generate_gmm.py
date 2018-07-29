from multiprocessing import Pool
import os


import argparse, sys
parser = argparse.ArgumentParser()
parser.add_argument("--input_map", help="full path to the input map for simulation",required=True)
parser.add_argument("--n_proc", help="number of threads to use",required=True)
parser.add_argument("--gmconvert_location", help="location of custom gmconvert",default="/home/jfraser/gmconvert/gmconvert")
parser.add_argument("--sbgrid_gmconvert_location", help="location of sbgrid gmconvert",default="/programs/x86_64-linux/gmconvert/20180516/bin/gmconvert")
args = parser.parse_args()

n_proc = int(args.n_proc)
input_map = args.input_map
gmconvert_location = args.gmconvert_location
sbgrid_gmconvert = args.sbgrid_gmconvert_location

# os.system("{gmconvert} -imap {map} -oimap output-0.0.mrc -zth 0.0".format(map=input_map, gmconvert=gmconvert_location))
# os.mkdir("ITER_0")
# os.mkdir("ITER_1")
# os.chdir("ITER_0")
# os.system("{gmconvert} V2G -imap ../output-0.0.mrc -ogmm 0.gmm -ng 150 -zth 0.0".format(gmconvert=gmconvert_location))
# os.chdir("../")

os.chdir("ITER_1")
# for i in range(0,150,1):
#     os.mkdir("map_0_{i}".format(i=i))
#     os.chdir("map_0_{i}".format(i=i))
#     os.system("echo \"../../ITER_0/0.gmm {i}\" > GMM.list".format(i=i))
#     os.chdir("../")
#
# def create_gmms(i):
#     os.chdir("map_0_{i}".format(i=i))
#     os.system("{gmconvert} -imap ../../{map} -gmml GMM.list -ng 100 -ogmm {i}_0.gmm -zth 0.0".format(map=input_map,i=i,gmconvert=gmconvert_location))
#     os.chdir("../")
#     return
#
# p = Pool(n_proc)
# print(p.map(create_gmms,range(0,150,1)))
#
# os.system("cat map_*/*.gmm > 1.gmm")
# os.system("{gmconvert} VcmpG -igmm 1.gmm -imap ../{map} -zth 0.0 -omap iter_1.mrc > iter_1.log".format(map=input_map, gmconvert=sbgrid_gmconvert))

from subprocess import check_output

try:
    from subprocess import STDOUT, check_output, CalledProcessError
except ImportError:  # pragma: no cover
    # python 2.6 doesn't include check_output
    # monkey patch it in!
    import subprocess
    STDOUT = subprocess.STDOUT

    def check_output(*popenargs, **kwargs):
        if 'stdout' in kwargs:  # pragma: no cover
            raise ValueError('stdout argument not allowed, '
                             'it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE,
                                   *popenargs, **kwargs)
        output, _ = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            raise subprocess.CalledProcessError(retcode, cmd,
                                                output=output)
        return output
    subprocess.check_output = check_output

    # overwrite CalledProcessError due to `output`
    # keyword not being available (in 2.6)
    class CalledProcessError(Exception):

        def __init__(self, returncode, cmd, output=None):
            self.returncode = returncode
            self.cmd = cmd
            self.output = output

        def __str__(self):
            return "Command '%s' returned non-zero exit status %d" % (
                self.cmd, self.returncode)
    subprocess.CalledProcessError = CalledProcessError


out = check_output(["grep CCgrid iter_1.log"],shell=True)

print("""

CORRELATION BETWEEN GMM and ORIGINAL MAP:

%s

run convert_GMM2PLUMED.sh to generate final dat file for PLUMED
""" %out)
