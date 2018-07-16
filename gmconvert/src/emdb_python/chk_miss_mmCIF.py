#!/usr/bin/env python
##
## <chk_miss_mmCIF.py>
##


import sys
import os
import random
from datetime import datetime
import math

import mmCIF

LastModDate = "Dec 19, 2014"

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




###############
#### MAIN #####
###############

OPT = {}

OPT['icifdir'] = '/DB/mmCIF'
OPT['ilist'] = ''
OPT['subout'] = 'T'
OPT['olist'] = 'miss_pdbids.list'

if (len(sys.argv)<3):
  print "chk_miss_mmCIF.py <options>"
  print " for checking mmCIF with missing GMMs."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -icifdir : input directory for mmCIF file [%s]"%(OPT['icifdir'])
  print " -ilist   : input file for PDB_ID list [%s]"%(OPT['ilist'])  
  print " -olist   : output file for list of missing PDB_IDs [%s]"%(OPT['olist']) 
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


### (3) do commands ###
MISS_PDB_ID = {}
for pdb_id in (pdb_id_list):
  iciffile = OPT['icifdir'] + '/' + pdb_id[1:3] + '/' + pdb_id + '.cif.gz'
  dat = {}
  
  mmCIF.read_mmCIF_file(iciffile,dat,focus_category=['pdbx_struct_assembly'])
  if ('pdbx_struct_assembly' in dat):
    assembly_id_list = dat['pdbx_struct_assembly']['id']
  else:
    assembly_id_list = []

  if (len(assembly_id_list)>1): 
    ogmmdir = OPT['odir'] + '/' + pdb_id[1:3]
    ogmmfile = ogmmdir + '/' + pdb_id + '-' + '1' + '.gmm' 
    if (os.path.isfile(ogmmfile)==0):
      assembly_id_list =  get_assembly_id_list(iciffile)
      print "#%s %s"%(pdb_id,assembly_id_list)

      for assembly_id in (assembly_id_list):

        ogmmdir = OPT['odir'] + '/' + pdb_id[1:3]
        if (assembly_id == '-'):
          ogmmfile = ogmmdir + '/' + pdb_id + '.gmm' 
        else:
          ogmmfile = ogmmdir + '/' + pdb_id + '-' + assembly_id + '.gmm' 

        if (os.path.isfile(ogmmfile)==0):
          print "#MISSING %s"%(ogmmfile)
          MISS_PDB_ID[pdb_id] = MISS_PDB_ID.get(pdb_id,'') + ' ' + assembly_id 


of = open(OPT['olist'],'w')
print "#write_missing_PDB_IDs() --> '%s'"%(OPT['olist'])
of.write("#COMMMAND %s\n"%(OPT['COMMAND']))
of.write("#DATE     %s\n"%(OPT['START_DATE']))
for pdb_id in (MISS_PDB_ID.keys()):
  of.write("%s %s\n"%(pdb_id,MISS_PDB_ID[pdb_id]))
of.close()

