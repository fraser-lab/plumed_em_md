#!/usr/bin/env python
import sys
import os
from datetime import datetime

LastModDate = "July 14, 2014"

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



def read_list_file(ifname,list,property):
# >> example <<
# EMD-5225 80s
# EMD-5326 80s
# EMD-5327 80s
# EMD-5328 80s
# EMD-5329 80s
# EMD-1042 c1
# EMD-1047 c1
# EMD-1080 c1
# EMD-1081 c1
# EMD-1095 c1
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
     if (len(fields)>=2):
       property[fields[0]] = fields[1]
  f.close()


###############
#### MAIN #####
###############

OPT = {}
OPT['A'] = 'F'
OPT['ilist'] = ''
OPT['icsv'] = ''
OPT['of'] = '-'

if (len(sys.argv)<3):
  print "suzuki2pair.py <options>"
  print " for suzuki's csv table into pair-score file."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -ilist : input list file [%s]"%(OPT['ilist'])
  print " -icsv  : input CSV file [%s]"%(OPT['icsv'])
  print " -A     : Action ('T' or 'F') [%s]"%(OPT['A'])
  print " -of    : output file [%s]"%(OPT['of'])
  sys.exit(1)

read_option(sys.argv,OPT)

### [1] Read id list ###
id_list = []
property = {} 
read_list_file(OPT['ilist'],id_list,property)
id_list_dic = {}
for id in (id_list):
  id_list_dic[id] = 1


### [2] Read icsv file ###
if (OPT['of']=='-'):
  of = sys.stdout
else:
  print "#write_output_pair_score_file() --> '%s'"%(OPT['of'])
  of = open(OPT['of'],'w')
  of.write("#COMMAND %s\n"%(OPT['COMMAND']))
  of.write("#DATE    %s\n"%(OPT['START_DATE']))
if (os.access(OPT['icsv'],os.R_OK)==0):
  print "#ERROR:can't open icsvfile '%s'"%(OPT['icsv'])
  sys.exit(1)

fi = open(OPT['icsv'])
Nline = 0
IDlist = []
for line in fi:
  line = line.rstrip('\n')
  if (Nline==0):
    list = line.split(',')
#emdb-5921,emdb-2604,pdb-3j3i,pdb-4cuv,
    IDlistCSV = []
    for x in (list):
      (head,tail) = x.split('-') 
      if (x.startswith('emdb-')):
         IDlistCSV.append('EMD-'+tail)
      if (x.startswith('pdb-')):
         IDlistCSV.append('PDB-'+tail)
  else:
    field = line.split(',')
    x = field[0]
    (head,tail) = x.split('-') 
    if (x.startswith('emdb-')):
      idA = 'EMD-'+tail
    if (x.startswith('pdb-')):
      idA = 'PDB-'+tail
    for i in range(1,len(field)):
      idB = IDlistCSV[i-1]
      if (idA in id_list_dic) and (idB in id_list_dic):
        of.write("%s %s %s\n"%(idA,idB,field[i]))
  Nline += 1
fi.close()

if (OPT['of']!='-'):
  print "#write_output_pair_score_file() --> '%s'"%(OPT['of'])
  of.close()
