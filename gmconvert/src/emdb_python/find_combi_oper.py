#!/usr/bin/env python
##
## <find_combi_oper.py>
##


import sys
import os
import random
from datetime import datetime
import math

import mmCIF

LastModDate = "Dec 15, 2014"

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



###############
#### MAIN #####
###############

OPT = {}

OPT['icifdir'] = '/home/takawaba/DB/mmCIF'
OPT['olist'] = 'out.list'

if (len(sys.argv)<3):
  print "find_combi_oper.py <options>"
  print " for pick up mmCIF files with combined oper_expression including ')('"
  print " such as, '(X0)(1-60)' or '(1-60)(61-88)'"
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -icifdir : input directory for mmCIF file [%s]"%(OPT['icifdir'])
  print " -olist   : output list file [%s]"%(OPT['olist'])
  sys.exit(1)

read_option(sys.argv,OPT)

read_option(sys.argv,OPT)


### (1) Get pdb_asmbl_id_list[] from OPT['ilist'] or the directory OPT['icifdir'] ###

pdbid_list = []

if (os.path.isdir(OPT['icifdir'])):
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
          dat = {}
          mmCIF.read_mmCIF_file(filefull,dat,focus_category=['pdbx_struct_assembly_gen'])
          findit = 0
          if ('pdbx_struct_assembly_gen' in dat):
            if ('oper_expression' in dat['pdbx_struct_assembly_gen']):
              for x in (dat['pdbx_struct_assembly_gen']['oper_expression']):
                if (x.find(')('))>0:
                  findit = 1
          if (findit==1):
            print "COMBI_OPER '%s'"%(filefull)
            field = file.split('.')
            pdbid_list.append(field[0])
                  
of = open(OPT['olist'],'w')
print "#output_pdbid_list() --> '%s'"%(OPT['olist'])
for pdbid in (pdbid_list):
  of.write("%s\n"%(pdbid))
of.close()



 
