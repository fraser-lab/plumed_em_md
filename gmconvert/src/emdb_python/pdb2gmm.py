#!/usr/bin/env python
##
## <pdb2gmm.py>
##


import sys
import os
import random
from datetime import datetime
from xml.etree.ElementTree import *
import math

LastModDate = "July 5, 2015"

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

###############
#### MAIN #####
###############

OPT = {}

OPT['ilist'] = ''
OPT['A'] = 'F'
OPT['div'] = '0/1'
OPT['ng'] = '10'
OPT['odir'] = 'tmpout/'
OPT['ipdbdir'] = '/DB/PDBv3'
OPT['emalg'] = 'G'
OPT['delzw'] = 'T'
OPT['delid'] = 'T'
OPT['atmrw'] = 'C'
OPT['subout'] = 'T'
OPT['update'] = 'F'
OPT['odirup'] = ''

if (len(sys.argv)<3):
  print "pdb2gmm.py <options>"
  print " for making GMMs from PDB atomic models."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -ilist  : input file for molecular list [%s]"%(OPT['ilist'])
  print " -ipdbdir : input PDB file directory [%s]"%(OPT['ipdbdir'])
  print " -odir   : output GMM directory [%s]"%(OPT['odir'])
  print " -subout : output with subdirectory ('T' or 'F') [%s]"%(OPT['subout'])
  print " -ng     : Ngauss [%s]"%(OPT['ng'])
  print " -atmrw : Model for radius and weight. 'A':atom model, 'R':residue model, 'C':decide from content. [%s]"%(OPT['atmrw'])  
  print " -emalg  : type for EM algorithm. 'P'oint_observed_EM,'G'mm_observe_EM [%s]"%(OPT['emalg'])
  print " -delzw : Delete Zero-weight gdfs from the GMM. ('T' or 'F') [%s]"%(OPT['delzw'])
  print " -delid : Delete identical gdfs in the GMM. ('T' or 'F') [%s]"%(OPT['delid'])
  print " -A      : Action ('T' or 'F') [%s]"%(OPT['A'])
  print " -div    : Job division (bunshi)/(bunbo)[%s]"%(OPT['div'])
  print " -update : make gmm only updating files (T or F) [%s]"%(OPT['update'])
  print " -odirup : output dir only for update GMMs [%s]"%(OPT['odirup'])
  sys.exit(1)

read_option(sys.argv,OPT)
[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)

read_option(sys.argv,OPT)


### (1) Get pdbid_list[] from OPT['ilist'] or the directory OPT['ipdbdir'] ###

idlist_orig = []
pdbid_list = []

if (OPT['ilist'] != ''):
  read_list_file(OPT['ilist'],idlist_orig)
  #>> examples of idlist[] <<
  #EMD-1895 30s
  #PDB-1p6g 30s
  #PDB-1p87 30s
  #EMD-1963 c2
  #PDB-3iyf c2
  for x in (idlist_orig):
    print x
    if (x.startswith('PDB-')):
      (head,pdbid) = x.split('-')
      pdbid_list.append(pdbid)
   # else:
   #   pdbid_list.append(x)
else:
  if (os.path.isdir(OPT['ipdbdir'])):
    dirlist = os.listdir(OPT['ipdbdir'])
    for subdir in (dirlist):
      dir = OPT['ipdbdir'] + '/' + subdir
      if (os.path.isfile(dir)):
        pass
      elif (os.path.isdir(dir)):
        filelist = os.listdir(dir)
        for file in (filelist):
          filefull = OPT['ipdbdir'] + '/' + subdir + '/' + file 
          if (os.path.isfile(filefull)):
            if (file.startswith('pdb')) and (file.endswith('.ent')):
              # ex) pdb4hhb.ent
              pdbid = file[3:7]
              #print pdbid  
              pdbid_list.append(pdbid)
print "#len(pdbid_list):%d"%(len(pdbid_list))


### (2) Make target_pdbid_list[] from pdbid_list ###
if (OPT['update'] == 'T'):
  if (os.path.isdir(OPT['odirup'])==0):
    print "#ERROR:update output directory (-odirup) '%s' does not exist."%(OPT['odirup'])
    sys.exit(1)

  ogmmdir_up = OPT['odirup'] + '/' + OPT['odir']
  if (os.path.isdir(ogmmdir_up)==0):
    if (OPT['A']=='T'):
      os.system("mkdir %s"%(ogmmdir_up)) 
    else:
      print "#mkdir %s"%(ogmmdir_up) 



target_pdbid_list = []
for pdbid in (pdbid_list):
  ipdbfile = OPT['ipdbdir'] + '/' + pdbid[1:3] + '/' + 'pdb' + pdbid + '.ent'
  if (OPT['subout']=='T'):
    ogmmdir = OPT['odir'] + '/' + pdbid[1:3]
    ogmmfile = ogmmdir + '/' + pdbid + '.gmm' 
  else:
    ogmmdir = OPT['odir']
    ogmmfile = ogmmdir + '/' + 'PDB-' + pdbid + '.gmm'

  #print pdbid,ipdbfile
  #print command
  accept_pdbid = 0
  if (OPT['update'] != 'T'):
    accept_pdbid = 1
  elif (OPT['update'] == 'T'):
    if (os.path.isfile(ogmmfile)==0):
      accept_pdbid = 1

  if (accept_pdbid == 1):
    target_pdbid_list.append(pdbid)
   
    if (os.path.isdir(ogmmdir)==0):
      if (OPT['A']=='T'):
        os.system("mkdir %s"%(ogmmdir)) 
      else:
        print "#mkdir %s"%(ogmmdir) 

    if (OPT['update']=='T'):
      if (OPT['subout']=='T'):
        ogmmdir_up = OPT['odirup'] + '/' + OPT['odir'] + '/' + pdbid[1:3]
      else:
        ogmmdir_up = OPT['odirup'] + '/' + OPT['odir'] 
      if (os.path.isdir(ogmmdir_up)==0):
        if (OPT['A']=='T'):
          os.system("mkdir %s"%(ogmmdir_up)) 
        else:
          print "#mkdir %s"%(ogmmdir_up) 



### (3) Assign Nstart,Nend by OPT['div']  ###
(bunshi,bunbo) = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)

Npdbid = len(target_pdbid_list);
Nstart   = bunshi*int(Npdbid/bunbo);
Nend     = (bunshi+1)*int(Npdbid/bunbo);
if (bunshi>=(bunbo-1)):
  Nend = Npdbid
print "#Npdbid %d bunshi/bunbo %d/%d start %d end %d"%(Npdbid,bunshi,bunbo,Nstart,Nend)

### (3) do commands ###
for i in range(Nstart,Nend):
  pdbid = target_pdbid_list[i]
  ipdbfile = OPT['ipdbdir'] + '/' + pdbid[1:3] + '/' + 'pdb' + pdbid + '.ent'
  if (OPT['subout']=='T'):
    ogmmdir = OPT['odir'] + '/' + pdbid[1:3]
    ogmmfile = ogmmdir + '/' + pdbid + '.gmm' 
  else:
    ogmmdir = OPT['odir']
    ogmmfile = ogmmdir + '/' + 'PDB-' + pdbid + '.gmm'
  
  command = "gmconvert -ipdb %s -ng %s -atmrw %s -emalg %s -delzw %s -delid %s -stdlog F -ogmm %s"%(ipdbfile,OPT['ng'],OPT['atmrw'],OPT['emalg'],OPT['delzw'],OPT['delid'],ogmmfile)
  print "#%s"%(command)
  if (OPT['A']=='T'):
    os.system(command)

  if (OPT['update']=='T'):
    if (OPT['subout']=='T'):
      ogmmdir_up  = OPT['odirup'] + '/' + OPT['odir'] + '/' + pdbid[1:3]
      ogmmfile_up = ogmmdir_up + '/' + pdbid + '.gmm' 
    else:
      ogmmdir_up  = OPT['odirup'] + '/' + OPT['odir'] 
      ogmmfile_up = ogmmdir_up + '/' + 'PDB-' + pdbid + '.gmm'
    cpcommand = "cp %s %s"%(ogmmfile,ogmmfile_up)
    print "#%s"%(cpcommand)
    if (OPT['A']=='T'):
      os.system(cpcommand)
