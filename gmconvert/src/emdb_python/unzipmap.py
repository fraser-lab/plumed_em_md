#!/usr/bin/env python
##
## <unzipmap.py>
##


import sys
import os
import random
from datetime import datetime

LastModDate = "Oct 23, 2012"

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
OPT['idir'] = '/DB/emdb/structures'
OPT['tail'] = '.map.gz'
OPT['subdir'] = 'map'
OPT['A'] = 'F'
OPT['avex'] = 'F'


if (len(sys.argv)<3):
  print "unzipmap.py <options>"
  print " for unzipping density map *.map.gz."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -idir   : input molecular file directory [%s]"%(OPT['idir'])
  print " -subdir : subdir under 'idir' [%s]"%(OPT['subdir'])
  print " -tail   : tail of map file [%s]"%(OPT['tail'])
  print " -avex   : avoid unzip, if the unzipped file already exists. ('T' or 'F') [%s]"%(OPT['tail'])
  print " -A      : Action ('T' or 'F') [%s]"%(OPT['A'])
  sys.exit(1)

read_option(sys.argv,OPT)

dirlist = os.listdir(OPT['idir'])
for dir in (dirlist):
  dirfull = OPT['idir'] + '/' + dir + '/' + OPT['subdir']
  if (os.path.isdir(dirfull)): 
    dirlist_under = os.listdir(dirfull)
    for file in (dirlist_under):
      if (file.endswith(OPT['tail'])):
#[y]               [yfull]
#'emd_1210.map.gz' '/DB/emdb/structures/EMD-1210/map/emd_1210.map.gz'
#'emd_5002.map.gz' '/DB/emdb/structures/EMD-5002/map/emd_5002.map.gz'
        field = file.split('.')
        ungzfile = ''
        for i in range(len(field)-1):
          if (i>0):
            ungzfile += '.'
          ungzfile += field[i]
        file_full     = OPT['idir'] + '/' + dir + '/' + OPT['subdir'] + '/' + file 
        ungzfile_full = OPT['idir'] + '/' + dir + '/' + OPT['subdir'] + '/' + ungzfile 
        if (OPT['avex']!='T') or (os.path.isfile(ungzfile_full)==0): 
          comstr = "zcat %s > %s"%(file_full,ungzfile_full)
          if (OPT['A']=='F'):
            print "#%s"%(comstr)
          if (OPT['A']=='T'):
            os.system(comstr)
