#!/usr/bin/env python
##
## <InGmm2atoms.py>
##


import sys
import os
import random
from datetime import datetime
import math
import gmm

LastModDate = "Oct 26, 2014"

def read_option(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")

  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


def read_list_file(ifname,list):
  print "#read_list_file('%s')"%(ifname)
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open filename '%s'" %(ifname)
    return(0)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('#')==0):
     #print line
     fields = line.split()
     list.append(fields[0])
  f.close()



def write_CA_CB_pdb_file(ofname,Dcacb):
#ATOM    124  CA  ALA A  15      22.093  16.971  17.520  1.00 12.79           C  
#ATOM    127  CB  ALA A  15      22.039  16.565  19.066  1.00 12.85           C  
  print "#write_CA_CB_pdb_file() --> '%s'"%(ofname)
  of = open(ofname,'w')
  of.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 12.79           C\n")
  of.write("ATOM      2  CB  ALA A   1    %8.3f   0.000   0.000  1.00 12.79           C\n"%(Dcacb))
  of.close()


###############
#### MAIN #####
###############

OPT = {}

OPT['odpdb'] = 'PDB'
OPT['odgmm'] = 'GMM'
OPT['dmin'] = '0.0'
OPT['dmax'] = '10.0'
OPT['ndiv'] = '10'
OPT['A'] = 'F'
OPT['of'] = 'summary.out'
OPT['emalg'] = 'G'
OPT['I']     = 'O'
OPT['varatm'] = 'A'
OPT['rr2var'] = '0.2'
OPT['resoatm'] = '0.0'

if (len(sys.argv)<3):
  print "InGmm2atoms.py <options>"
  print " for testing GMMinput-EM algorithm by testing by input two atoms."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -odpdb  : output dir for generated PDB for two atoms. [%s]"%(OPT['odpdb'])
  print " -odgmm  : output dir for generated GMM for two gdfs.  [%s]"%(OPT['odgmm'])
  print " -dmin   : minimum distance  [%s]"%(OPT['dmin'])
  print " -dmax   : maximum distance  [%s]"%(OPT['dmax'])
  print " -ndiv   : division          [%s]"%(OPT['ndiv'])
  print " -of     : output file for summary [%s]"%(OPT['of'])
  print " -A      : Action (T or F) [%s]"%(OPT['A'])
  print "<options for gmconvert>"
  print " -emalg  : type for EM algorithm. 'P'oint_observed_EM,'G'mm_observe_EM, 'O':one-to-one_atom/voxel [%s]"%(OPT['emalg'])
  print " -I      : Initialization of GMM. 'K'-means, 'R'andom 'O':one-to-one_atom/voxel [%s]"%(OPT['I'])
  print " -varatm : Variance type for atom (for -emalg G). 'A': var = rr2var * Rvdw*Rvdw for each atom, 'R': var = (resoatm/2.0)^2.[%s]"%(OPT['varatm'])
  print " -rr2var : Constant for variance = Const * Rvdw*Rvdw for -emalg G -varatm A. Default is 1/5 [%s]."%(OPT['rr2var'])
  print " -resoatm: Resolution for atom for -emalg G -varatm R.  [%s]."%(OPT['resoatm'])


  sys.exit(1)

read_option(sys.argv,OPT)
ndiv = int(OPT['ndiv'])
dmin = float(OPT['dmin'])
dmax = float(OPT['dmax'])

DATA = []
for i in range(0,ndiv+1):
  a = {}
  dis = dmin + i*(dmax-dmin)/ndiv
  print "#[%d] dis %f"%(i,dis)
  opdbfile = OPT['odpdb'] + '/' + "%d.pdb"%(i)
  ogmmfile = OPT['odgmm'] + '/' + "%d.gmm"%(i)
  command = "gmconvert -ipdb %s -ogmm %s -emalg %s -I %s -varatm %s -rr2var %s -resoatm %s -ng 2 -ccatm T -delzw F -delid F"%(opdbfile,ogmmfile,OPT['emalg'],OPT['I'],OPT['varatm'],OPT['rr2var'],OPT['resoatm'])
  print "#command '%s'"%(command)
  a['Datom'] = dis
  a['command'] = command 
  if (OPT['A']=='T'): 
    write_CA_CB_pdb_file(opdbfile,dis)
    fg = os.popen(command)
    for line in fg:
# NGAUSS 2 CorrCoeff 1.000000 logLike -8.827242e+00 Sigma 5.000000
      if (line.startswith("NGAUSS")):
        line = line.rstrip('\n')
        fields = line.split()
        a['CorrCoeff'] = fields[3] 
        a['logLike']   = fields[5] 
    fg.close()
    g = gmm.GAUSS_MOLECULE()
    g.read(ogmmfile)
    Dgauss = 0.0
    for k in range(3):
      dx = g.gauss[0].M[k] - g.gauss[1].M[k]
      Dgauss += dx*dx
    Dgauss = math.sqrt(Dgauss)
    a['Dgauss'] = Dgauss

  DATA.append(a)


### [2] Output summary file ###
if (OPT['A']=='T'): 
  print "#write_summary() --> '%s'"%(OPT['of'])
  of = open(OPT['of'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#DATE    %s\n"%(OPT['START_DATE']))
  of.write("#EXAMPLE_FOR_GMCONVERT_COMMAND %s\n"%(DATA[0]['command']))
  of.write("#[num] [Datom] [Dgauss] [CorrCoeff] [logLike]\n")
  for i in range(0,ndiv+1):
    a = DATA[i] 
    of.write("%d %f %f %s %s\n"%(i,a['Datom'],a['Dgauss'],a['CorrCoeff'],a['logLike']))
  of.close()



