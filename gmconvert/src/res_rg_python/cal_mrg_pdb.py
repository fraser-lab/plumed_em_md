#!/usr/bin/python
#
# <cal_mrg_pdb.py>
#

import sys
import os
import math
from datetime import datetime
#sys.path.append('/home/takawaba/work/gmconvert/src/mrg_python')
import pdb 
import radius

LastModDate = "Oct 18, 2014"

def read_option(argv,opt_dic):
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]



############
### MAIN ###
############

OPT = {}
OPT['ipdb'] = ''

if (len(sys.argv)<2):
  print "cal_mrg_pdb.py <options>"
  print " for calculating molecular weight and radius of gyration."
  print " coded by T.Kawabata.  LastModified:%s"%(LastModDate)
  print "<options>"
  print " -ipdb  : input PDB file [%s]"%(OPT['ipdb'])
  sys.exit(1)


read_option(sys.argv,OPT)
RADIUS = {}
radius.set_vdw_radius_dic(RADIUS) 
P = pdb.Molecule()
P.read(OPT['ipdb'])
radius.set_radius_to_pdb(RADIUS,P)
sum_weight = radius.set_weight_to_pdb(P)

for a in range(P.Natom):
  print "#'%s' '%s' radius %.2f weight %.2f"%(P.resi[a],P.atom[a],P.radius[a],P.weight[a])

print "#sum_weight %f"%(sum_weight)


### calculate covariance matrix ###
CovM = [[0.0 for i in range(3)] for j in range(3)]
M = [0.0 for i in range(3)]
for a in range(P.Natom):
  w = P.weight[a]/sum_weight
  for i in range(3):
    M[i] += w * P.Pos[a][i] 
print "#M %f %f %f"%(M[0],M[1],M[2])
for a in range(P.Natom):
  w = P.weight[a]/sum_weight
  variance = P.radius[a]*P.radius[a]/5.0

  for i in range(3):
    for j in range(3):
      if (i==j):
        CovM[i][j] += w * (variance + (P.Pos[a][i]-M[i])*(P.Pos[a][j]-M[j]))
      else: 
        CovM[i][j] += w * (P.Pos[a][i]-M[i])*(P.Pos[a][j]-M[j])
  
print "#CovM %f %f %f"%(CovM[0][0],CovM[0][1],CovM[0][2])
print "#CovM %f %f %f"%(CovM[1][0],CovM[1][1],CovM[1][2])
print "#CovM %f %f %f"%(CovM[2][0],CovM[2][1],CovM[2][2])

Rg_pnt = P.radius_of_gylation_gcen()

Rg = math.sqrt(CovM[0][0]+CovM[1][1]+CovM[2][2])

print "#Rg %f Rg_pnt %f"%(Rg,Rg_pnt)
print "%-14s sum_weight %10.5f Rg %f"%(OPT['ipdb'],sum_weight,Rg)
