import os

try:
    from subprocess import STDOUT, check_output, CalledProcessError
except ImportError:  # pragma: no cover
    # python 2.6 doesn't include check_output
    # monkey patch it in!
    # from: https://stackoverflow.com/questions/4814970/subprocess-check-output-doesnt-seem-to-exist-python-2-6-5
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
