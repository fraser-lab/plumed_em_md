#!/usr/bin/env python
##
## <mmCIF2gmm.py>
##


import sys
import os
import random
from datetime import datetime
import math

import mmCIF

LastModDate = "July 15, 2015"

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


def get_assembly_id_list(iciffile):
  assembly_id_list = ['-']
  dat = {}

  mmCIF.read_mmCIF_file(iciffile,dat,focus_category=['atom_site','pdbx_struct_assembly','pdbx_struct_assembly_gen'])

  if ('atom_site' not in dat):
    return([])    
  if ('pdbx_struct_assembly' not in dat) or ('pdbx_struct_assembly_gen' not in dat):
    return(['-'])    

  S = dat['atom_site']
  asym_id_dic = {}
  for i in range(len(S['id'])):
    asym_id_dic[S['label_asym_id'][i]] = 1

  asym_id_list_asymmetric_unit = sorted(asym_id_dic.keys(),lambda x,y:cmp(x,y))

  A = dat['pdbx_struct_assembly']
  G = dat['pdbx_struct_assembly_gen']
  Nassembly = len(A['id'])
  items = ('id','details','method_details','oligomeric_details','oligomeric_count')
  for i in range(Nassembly):
    assembly_id = A['id'][i]

    for item in (items):
      if (item in A.keys()):
        sys.stdout.write(" %s"%(A[item][i]))
      else:
        sys.stdout.write("-")
    sys.stdout.write(" %s"%(G['asym_id_list'][i]))
    sys.stdout.write(" %s"%(G['oper_expression'][i]))

    oper_expression_list = mmCIF.get_oper_expression_list(G['oper_expression'][i])
    asym_id_list = sorted(G['asym_id_list'][i].split(','),lambda x,y:cmp(x,y))

    ## Check the assembly is equal to the asymmetric unit,or not
    ## (1) asym_id_list[] is identical asym_id_list_asymmetric_unit[]   
    ## (2) Length of oper_expression_list[] is 1. (Single transformation)

    if (asym_id_list==asym_id_list_asymmetric_unit) and (len(oper_expression_list)==1):
      sys.stdout.write("#assembly_id '%s' is the asymmetric_unit."%(assembly_id))
    else:
      assembly_id_list.append(assembly_id)

    sys.stdout.write("\n")

  dat = {}

  return(assembly_id_list)



def make_ogmmdir_ogmmfile_name(odir,subout,pdb_id,assembly_id):
  if (subout=='T'):
    ogmmdir = odir + '/' + pdb_id[1:3]
  else:
    ogmmdir = odir

  if (assembly_id == '') or (assembly_id == '-'):
    ogmmfile = ogmmdir + '/' + pdb_id + '.gmm' 
  else:
    ogmmfile = ogmmdir + '/' + pdb_id + '-' + assembly_id + '.gmm' 
  
  return((ogmmdir,ogmmfile))



###############
#### MAIN #####
###############

OPT = {}

OPT['A'] = 'F'
OPT['div'] = '0/1'
OPT['ng'] = '10'
OPT['odir'] = 'tmpout/'
OPT['icifdir'] = '/home/takawaba/DB/mmCIF'
OPT['ilist'] = ''
OPT['emalg'] = 'G'
OPT['delzw'] = 'T'
OPT['delid'] = 'T'
OPT['subout'] = 'T'
OPT['update'] = 'F'
OPT['odirup'] = ''
OPT['atmsel'] = 'A'
OPT['atmrw']  = 'C'
OPT['maxatm'] = '100000'
OPT['avex'] = 'F'
OPT['justcnt'] = 'F'

if (len(sys.argv)<3):
  print "mmCIF2gmm.py <options>"
  print " for making GMMs from PDB atomic models described in mmCIF file."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -icifdir : input directory for mmCIF file [%s]"%(OPT['icifdir'])
  print " -ilist   : input file for PDB_ID list [%s]"%(OPT['ilist'])  
  print " -odir    : output GMM directory [%s]"%(OPT['odir'])
  print " -subout  : output with subdirectory ('T' or 'F') [%s]"%(OPT['subout'])
  print " -A      : Action ('T' or 'F') [%s]"%(OPT['A'])
  print " -div    : Job division (bunshi)/(bunbo)[%s]"%(OPT['div'])
  print " -avex   : avoid pdb_id if its GMM files are existed (T or F)[%s]"%(OPT['avex'])
  print " -justcnt: stop after just counting number of entries for GMM conversion. (T or F)[%s]"%(OPT['justcnt']) 
  print "<options for update>"
  print " -update : make gmm only updating files (T or F) [%s]"%(OPT['update'])
  print " -odirup : output dir only for update GMMs [%s]"%(OPT['odirup'])
  print "<options for 'gmconvert'>"
  print " -ng      : Ngauss [%s]"%(OPT['ng'])
  print " -atmsel : Atom selection 'A'll atom except hydrogen, 'R'esidue-based (only ' CA ' and ' P ') [%s]"%(OPT['atmsel'])
  print " -atmrw  : Model for radius and weight. 'A':atom model, 'R':residue model, 'C':decide from content. [%s]"%(OPT['atmrw'])
  print " -maxatm : maximum allowed number of atoms for '-atmsel A'. If over '-minatmA', then change '-atmsel R'.[%s]"%(OPT['maxatm'])
  print " -emalg  : type for EM algorithm. 'P'oint_observed_EM,'G'mm_observe_EM [%s]"%(OPT['emalg'])
  print " -delzw  : Delete Zero-weight gdfs from the GMM. ('T' or 'F') [%s]"%(OPT['delzw'])
  print " -delid  : Delete identical gdfs in the GMM. ('T' or 'F') [%s]"%(OPT['delid'])
  sys.exit(1)

read_option(sys.argv,OPT)

### (1) Get pdb_id_list[] from OPT['ilist'] or the directory OPT['icifdir'] ###

pdb_id_list = []

if (OPT['ilist'] != ''):
  read_list_file(OPT['ilist'],pdb_id_list)
elif (os.path.isdir(OPT['icifdir'])):
  dirlist = os.listdir(OPT['icifdir'])
  for subdir in (dirlist):
    dir = OPT['icifdir'] + '/' + subdir
    if (os.path.isfile(dir)):
      pass
    elif (os.path.isdir(dir)):
      filelist = os.listdir(dir)
      for file in (filelist):
        filefull = OPT['icifdir'] + '/' + subdir + '/' + file 
        if (os.path.isfile(filefull)):
          if (file.endswith('.gz')) and (len(file)>=1):
# >> file example <<
#1ai0.cif.gz
#1ai1.cif.gz
#1ai2.cif.gz
            field  = file.split('.')
            pdb_id = field[0]
            pdb_id_list.append(pdb_id)

print "#len(pdb_id_list):%d"%(len(pdb_id_list))

#### (2) When OPT['avex']=='T' (avoiding preexisting GMM files) ##
if (OPT['avex']=='T'):
  pdb_id_list_orig = []
  for pdb_id in (pdb_id_list):
    pdb_id_list_orig.append(pdb_id)

  pdb_id_list = []

  for pdb_id in (pdb_id_list_orig):
    (ogmmdir,ogmmfile) = make_ogmmdir_ogmmfile_name(OPT['odir'],OPT['subout'],pdb_id,'')
    if (os.path.isfile(ogmmfile)==0):
      pdb_id_list.append(pdb_id)
      print pdb_id
  
  print "#len(pdb_id_list):%d"%(len(pdb_id_list))

### (2) for the case for update ###
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

  pdb_id_list_orig = []
  for pdb_id in (pdb_id_list):
    pdb_id_list_orig.append(pdb_id)

  pdb_id_list = []

  for pdb_id in (pdb_id_list_orig):
    (ogmmdir,ogmmfile) = make_ogmmdir_ogmmfile_name(OPT['odir'],OPT['subout'],pdb_id,'')
    ogmmfile_date = 0
    ogmmfile_dt = ''
    if (os.path.isfile(ogmmfile)):
      ogmmfile_date = os.stat(ogmmfile).st_mtime
      ogmmfile_dt = datetime.fromtimestamp(os.stat(ogmmfile).st_mtime).strftime('%Y-%m-%d')

    iciffile = OPT['icifdir'] + '/' + pdb_id[1:3] + '/' + pdb_id + '.cif.gz'
    iciffile_date = 0
    icifpfile_dt = ''
    if (os.path.isfile(iciffile)):
      iciffile_date = os.stat(iciffile).st_mtime
      iciffile_dt = datetime.fromtimestamp(os.stat(iciffile).st_mtime).strftime('%Y-%m-%d')

    diff_date = iciffile_date - ogmmfile_date
    if (diff_date>0):
      print "pdb_id %s mmcif_date %s gmm_date %s diff_date %f cif '%s' gmm '%s'"%(pdb_id,iciffile_date,ogmmfile_date,iciffile_date-ogmmfile_date, iciffile_dt, ogmmfile_dt)
      pdb_id_list.append(pdb_id)


### (3) Assign Nstart,Nend by OPT['div']  ###
(bunshi,bunbo) = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)

Npdb_id = len(pdb_id_list);
Nstart   = bunshi*int(Npdb_id/bunbo);
Nend     = (bunshi+1)*int(Npdb_id/bunbo);
if (bunshi>=(bunbo-1)):
  Nend = Npdb_id
print "#Npdb_id %d bunshi/bunbo %d/%d start %d end %d"%(Npdb_id,bunshi,bunbo,Nstart,Nend)

if (OPT['justcnt']=='T'):
  print "#STOP after counting entries to be updated."
  sys.exit(1)

### (3) do commands ###
for i in range(Nstart,Nend):
  #pdb_id = fields[0]
  pdb_id = pdb_id_list[i]
  iciffile = OPT['icifdir'] + '/' + pdb_id[1:3] + '/' + pdb_id + '.cif.gz'

  #dat = {}
  #mmCIF.read_mmCIF_file(iciffile,dat,focus_category=['pdbx_struct_assembly'])
  
  assembly_id_list =  get_assembly_id_list(iciffile)

  print "#%s %s"%(pdb_id,assembly_id_list)

  for assembly_id in (assembly_id_list):
    (ogmmdir,ogmmfile) = make_ogmmdir_ogmmfile_name(OPT['odir'],OPT['subout'],pdb_id,assembly_id)

    if (os.path.isdir(ogmmdir)==0):
      if (OPT['A']=='T'):
        os.system("mkdir %s"%(ogmmdir))
      else:
        print "#mkdir %s"%(ogmmdir)

    command = "gmconvert -icif %s"%(iciffile)
    if (assembly_id != '-'):
      command += " -assembly %s"%(assembly_id)
    command += " -model M -ng %s -atmsel %s -atmrw %s -emalg %s -delzw %s -delid %s -stdlog F"%(OPT['ng'],OPT['atmsel'],OPT['atmrw'],OPT['emalg'],OPT['delzw'],OPT['delid'])
    if (int(OPT['maxatm'])>0):
      command += " -maxatm %s"%(OPT['maxatm'])
    command += " -ogmm %s"%(ogmmfile)

    print "#%s"%(command)

    if (OPT['A']=='T'):
      os.system(command)


    if (OPT['odirup'] != ''):
      (ogmmdir_up,ogmmfile_up) = make_ogmmdir_ogmmfile_name(OPT['odirup'] + '/' + OPT['odir'],OPT['subout'],pdb_id,assembly_id)
      if (os.path.isdir(ogmmdir_up)==0):
        if (OPT['A']=='T'):
          os.system("mkdir %s"%(ogmmdir_up))
        else:
          print "#mkdir %s"%(ogmmdir_up)
      upcommand = "cp %s %s"%(ogmmfile,ogmmfile_up)
      print "#%s"%(upcommand)
      if (OPT['A']=='T'):
        os.system(upcommand)


