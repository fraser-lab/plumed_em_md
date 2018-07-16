#!/usr/bin/env python
##
## <biounit2gmm.py>
##


import sys
import os
import random
from datetime import datetime
from xml.etree.ElementTree import *
import math

LastModDate = "Sep 21, 2014"

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
OPT['ipdbdir'] = '/home/takawaba//DB/PDBbiounit'
OPT['emalg'] = 'G'
OPT['delzw'] = 'T'
OPT['delid'] = 'T'
OPT['subout'] = 'T'
OPT['update'] = 'F'
OPT['odirup'] = ''
OPT['atmsel'] = 'A'
OPT['atmrw']  = 'C'
OPT['maxatm'] = '100000'

if (len(sys.argv)<3):
  print "biounit2gmm.py <options>"
  print " for making GMMs from PDB atomic models of Biological Unit."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -ilist  : input file for molecular list [%s]"%(OPT['ilist'])
  print " -ipdbdir : input PDB file directory for biological units [%s]"%(OPT['ipdbdir'])
  print " -odir   : output GMM directory [%s]"%(OPT['odir'])
  print " -subout : output with subdirectory ('T' or 'F') [%s]"%(OPT['subout'])
  print " -ng     : Ngauss [%s]"%(OPT['ng'])
#  print " -atmres : Model for radius and weight. 'A':atom model, 'R':residue model, 'C':decide from content. [%s]"%(OPT['atmres'])  
  print " -atmsel : Atom selection 'A'll atom except hydrogen, 'R'esidue-based (only ' CA ' and ' P ') [%s]"%(OPT['atmsel'])
  print " -atmrw  : Model for radius and weight. 'A':atom model, 'R':residue model, 'C':decide from content. [%s]"%(OPT['atmrw'])
  print " -maxatm : maximum allowed number of atoms for '-atmsel A'. If over '-minatmA', then change '-atmsel R'.[%s]"%(OPT['maxatm'])
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


### (1) Get pdb_asmbl_id_list[] from OPT['ilist'] or the directory OPT['ipdbdir'] ###

idlist_orig = []
pdb_asmbl_id_list = []

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
      (head,pdb_id) = x.split('-')
      pdb_asmbl_id_list.append(pdb_id)
   # else:
   #   pdb_asmbl_id_list.append(x)
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
            if (file.endswith('.gz')) and (len(file)>=12):
            # ex)
            # 2hhg.pdb1.gz
            # 2hhg.pdb2.gz
            # 4njn.pdb1.gz
            # 4njn.pdb2.gz
            # 4njn.pdb3.gz
            # 4njn.pdb4.gz
            # 2zzs.pdb1.gz
            # 2zzs.pdb2.gz
            # :
            # 2zzs.pdb10.gz
            # 2zzs.pdb11.gz
            # :
            # 2zzs.pdb28.gz
            # 2zzs.pdb29.gz
              (pdb_id,tail,gz) = file.split('.')
              #pdb_id       = file[0:4]
              asmbl_id = tail[3:]
              pdb_asmbl_id = pdb_id + '-' + asmbl_id
              #print "%s '%s' '%s'"%(file,pdb_id,asmbl_id)
              #print pdb_id  
              pdb_asmbl_id_list.append(pdb_asmbl_id)
print "#len(pdb_asmbl_id_list):%d"%(len(pdb_asmbl_id_list))


### (2) Make target_pdb_asmbl_id_list[] from pdb_asmbl_id_list ###
if (OPT['update'] == 'T'):
  ogmmdir_up = OPT['odirup'] + '/' + OPT['odir']
  if (os.path.isdir(ogmmdir_up)==0):
    if (OPT['A']=='T'):
      os.system("mkdir %s"%(ogmmdir_up)) 
    else:
      print "#mkdir %s"%(ogmmdir_up) 



target_pdb_asmbl_id_list = []
for pdb_asmbl_id in (pdb_asmbl_id_list):
  (pdb_id,asmbl_id) = pdb_asmbl_id.split('-') 
  ipdbfile = OPT['ipdbdir'] + '/' + pdb_id[1:3] + '/' + pdb_id + '.pdb' + asmbl_id + '.gz'
  if (OPT['subout']=='T'):
    ogmmdir = OPT['odir'] + '/' + pdb_id[1:3]
    ogmmfile = ogmmdir + '/' + pdb_id + '-' + asmbl_id + '.gmm' 
  else:
    ogmmdir = OPT['odir']
    ogmmfile = ogmmdir + '/' + 'PDB-' + pdb_id + '-' + asmbl_id + '.gmm'

  #print pdb_id,ipdbfile
  #print command
  accept_pdb_id = 0
  if (OPT['update'] != 'T'):
    accept_pdb_id = 1
  elif (OPT['update'] == 'T'):
    if (os.path.isfile(ogmmfile)==0):
      accept_pdb_id = 1

  if (accept_pdb_id == 1):
    target_pdb_asmbl_id_list.append(pdb_asmbl_id)
   
    if (os.path.isdir(ogmmdir)==0):
      if (OPT['A']=='T'):
        os.system("mkdir %s"%(ogmmdir)) 
      else:
        print "#mkdir %s"%(ogmmdir) 

    if (OPT['update']=='T'):
      if (OPT['subout']=='T'):
        ogmmdir_up = OPT['odirup'] + '/' + OPT['odir'] + '/' + pdb_id[1:3]
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

Npdb_id = len(target_pdb_asmbl_id_list);
Nstart   = bunshi*int(Npdb_id/bunbo);
Nend     = (bunshi+1)*int(Npdb_id/bunbo);
if (bunshi>=(bunbo-1)):
  Nend = Npdb_id
print "#Npdb_id %d bunshi/bunbo %d/%d start %d end %d"%(Npdb_id,bunshi,bunbo,Nstart,Nend)

### (3) do commands ###
for i in range(Nstart,Nend):
  pdb_asmbl_id = target_pdb_asmbl_id_list[i]
  (pdb_id,asmbl_id) = pdb_asmbl_id.split('-') 
  ipdbfile = OPT['ipdbdir'] + '/' + pdb_id[1:3] + '/' + pdb_id + '.pdb' + asmbl_id + '.gz'
  if (OPT['subout']=='T'):
    ogmmdir = OPT['odir'] + '/' + pdb_id[1:3]
    ogmmfile = ogmmdir + '/' + pdb_id + '-' + asmbl_id + '.gmm' 
  else:
    ogmmdir = OPT['odir']
    ogmmfile = ogmmdir + '/' + 'PDB-' + pdb_id + '-' + asmbl_id + '.gmm'
  
  command = "gmconvert -ipdb %s"%(ipdbfile)
  command += " -model M -ng %s -atmsel %s -atmrw %s -emalg %s -delzw %s -delid %s -stdlog F"%(OPT['ng'],OPT['atmsel'],OPT['atmrw'],OPT['emalg'],OPT['delzw'],OPT['delid'])
  if (int(OPT['maxatm'])>0):
    command += " -maxatm %s"%(OPT['maxatm'])
  command += " -ogmm %s"%(ogmmfile)

  print "#%s"%(command)
  if (OPT['A']=='T'):
    os.system(command)

  if (OPT['update']=='T'):
    if (OPT['subout']=='T'):
      ogmmdir_up  = OPT['odirup'] + '/' + OPT['odir'] + '/' + pdb_id[1:3]
      ogmmfile_up = ogmmdir_up + '/' + pdb_id + '.gmm' 
    else:
      ogmmdir_up  = OPT['odirup'] + '/' + OPT['odir'] 
      ogmmfile_up = ogmmdir_up + '/' + 'PDB-' + pdb_id + '.gmm'
    cpcommand = "cp %s %s"%(ogmmfile,ogmmfile_up)
    print "#%s"%(cpcommand)
    if (OPT['A']=='T'):
      os.system(cpcommand)
