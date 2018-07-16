#!/usr/bin/env python

import sys
import os
from datetime import datetime

LastModDate = "Oct 21, 2013"

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


def read_pdb_xyz_coodinates(ifname,atomdat):
  print "#read_pdb_xyz_coodinates('%s')"%(ifname)
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open '%s'"%(ifname)
    sys.exit(1)
  f = open(ifname)
  for line in f:
    if (line.startswith('ATOM')) or (line.startswith('HETATM')):
      line = line.rstrip('\n')
      #print line
      anum  = line[6:11]
      atom  = line[12:16]
      resi  = line[17:20]
      chain = line[21:22]
      rnum  = line[22:27]
      index = atom + ':' + resi + ':' + rnum
      #print index
      atomdat[index] = {}
      atomdat[index]['x'] = line[31:38]
      atomdat[index]['y'] = line[38:46]
      atomdat[index]['z'] = line[46:54]
  f.close()

##############
#### MAIN ####
##############

OPT = {}

OPT['impdb'] = ''
OPT['irpdb'] = ''
OPT['ompdb'] = ''

if (len(sys.argv)<2):
  print "repxyz_mpdb.py <options>"
  print " for replacing xyz of membership_pdb file."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -impdb: input  membership pdb file [%s]"%(OPT['impdb'])
  print " -irpdb: input  pdb file containing replacing XYZ coordinates[%s]"%(OPT['irpdb'])
  print " -ompdb: output membership pdb file [%s]"%(OPT['ompdb'])
  sys.exit(1)

read_option(sys.argv,OPT)

### [1] read OPT['irpdb'] and store as repxyz ###
repxyz = {}
read_pdb_xyz_coodinates(OPT['irpdb'],repxyz)


### [2] read OPT['impdb'] and write to OPT['ompdb']  ###
if not os.access(OPT['impdb'],os.R_OK):
  print "#ERROR:Can't open '%s'"%(OPT['impdb'])
  sys.exit(1)

fi = open(OPT['impdb']) 
print "#open membership_pdb_file('%s')"%(OPT['impdb'])

if (OPT['ompdb']==''):
  fo = sys.stdout
else:
  fo = open(OPT['ompdb'],'w')
  print "#write_replace_xyz_memberships_pdb_file() --> '%s'"%(OPT['ompdb'])

fo.write("REMARK   COMMAND %s\n"%(OPT['COMMAND']))
fo.write("REMARK   DATE    %s\n"%(OPT['START_DATE']))

for line in fi:
  if (line.startswith('ATOM')) or (line.startswith('HETATM')):
    line = line.rstrip('\n')

#REMARK                                   MEMBERSHIP VALUE FOR GDF: [  1] [  2]
#          1         2         3         4         5         6         7 
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#ATOM      1  N   ALA A   2      63.130 -67.331 -18.137  1.65  1.00 1.000 0.000
#ATOM      2  CA  ALA A   2      63.914 -67.078 -16.894  1.87  1.00 1.000 0.000
#ATOM      3  C   ALA A   2      63.317 -67.813 -15.708  1.76  1.00 1.000 0.000
#ATOM      4  O   ALA A   2      62.179 -68.283 -15.768  1.40  1.00 1.000 0.000
#ATOM      5  CB  ALA A   2      63.950 -65.601 -16.605  1.87  1.00 1.000 0.000
    atom  = line[12:16]
    resi  = line[17:20]
    rnum  = line[22:27]
    index = atom + ':' + resi + ':' + rnum
    if (index in repxyz):
      head = line[0:30] 
      tail = line[54:] 
      xyzstr = "%8s%8s%8s"%(repxyz[index]['x'],repxyz[index]['y'],repxyz[index]['z'])
      fo.write("%s%s%s\n"%(head,xyzstr,tail))
  else:
    fo.write("%s"%(line))

fi.close()

if (OPT['ompdb']!=''):
  fo.close()



