#!/usr/bin/env python
##
## <emdb_vol_pca.py>
##


import sys
import os
import random
from datetime import datetime
from xml.etree.ElementTree import *
import math

LastModDate = "Aug 4, 2015"

def read_option(argv,opt_dic):
  opt_dic['start_datetime_now'] = datetime.now()
  opt_dic['START_DATE'] = opt_dic['start_datetime_now'].strftime("%Y/%m/%d %H:%M:%S")
  opt_dic['COMMAND'] = argv[0]
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



def read_XML_file(ifname,dat):
  contents = ''
  f = open(ifname)
  for line in f:
    contents += line
  f.close()
  #print "'%s'"%(contents)
  elem = fromstring(contents)
  acc = elem.get('accessCode').rstrip('\n')
  #print "#acc '%s'"%(acc)
  pdbEntryId = ''
  fittedPDBEntryId = ''
  for e in elem.getiterator():
    if (e.tag.startswith('title')):
      #print ">%s :: '%s'"%(e.tag,e.text)
      string       = e.text.replace('\n',' ')
      dat['title'] = string.replace('\t',' ')
# <numColumns>256</numColumns>
# <numRows>256</numRows>
# <numSections>256</numSections>
    if (e.tag == 'numColumns'):
      dat['numColumns'] = e.text.rstrip('\n')
    if (e.tag == 'numRows'):
      dat['numRows'] = e.text.rstrip('\n')
    if (e.tag == 'numSections'):
      dat['numSections'] = e.text.rstrip('\n')
    if (e.tag == 'resolutionByAuthor'):
      dat['resolutionByAuthor'] = e.text.rstrip('\n')
    if (e.tag == 'contourLevel'):
      dat['contourLevel'] = e.text.rstrip('\n')
    if (e.tag == 'authors'):
      dat['authors'] = e.text.rstrip('\n')
    if (e.tag == 'keywords'):
      dat['keywords'] = e.text
    if (e.tag == 'entry'):
      dat['entry'] = dat.get('entry','') + ':' + e.text.rstrip('\n')
    if (e.tag == 'pdbEntryId'):
      #print ">%s :: '%s'"%(e.tag,e.text)
      id = e.text
      id = id.rstrip('\n')
      if (pdbEntryId != ''):
        pdbEntryId += ':'
      pdbEntryId += id
    if (e.tag == 'fittedPDBEntryId'):
      #print ">%s :: '%s'"%(e.tag,e.text)
      id = e.text
      id = id.rstrip('\n')
      if (fittedPDBEntryId != ''):
        fittedPDBEntryId += ':'
      fittedPDBEntryId += id
  dat['acc'] = acc
  dat['acc'] = acc
  if (pdbEntryId==''):
    pdbEntryId = 'NO-PDB'
  if (fittedPDBEntryId==''):
    fittedPDBEntryId = 'NO-FIT'
  dat['pdbEntryId'] = pdbEntryId 
  dat['fittedPDBEntryId'] = fittedPDBEntryId 
  if (dat.has_key('resolutionByAuthor')==0):
    dat['resolutionByAuthor'] = 'NORES'

def execute_gmconvert(command):
# >> EXAMPLE of output of 'gmconvert' 
# #COMMAND gmconvert -imap /home/takawaba/DB/emdb/structures/EMD-2503/map/emd_2503.map.gz -ng 0 -ogmm temp.gmm -zth 1.7
# #MODE:V2G
# #Make_Exp_minus_X_Table(MaxX 20.000000 Width 0.010000 Ndivide 2000)
# #Read_Density_File("/home/takawaba/DB/emdb/structures/EMD-2503/map/emd_2503.map.gz" FileType 'C')
# #NC 400 Order N
# #MODE 2
# ANGLE_Alpha 90.000000
# ANGLE_Beta  90.000000
# ANGLE_Gamma 90.000000
# GRID_SIZE_X 400
# GRID_SIZE_Y 400
# GRID_SIZE_Z 400
# GRID_WIDTH 1.196000
# #Set_Small_Value_Voxel_to_Zero(threshold 1.700000 Nzero 60390921/64000000 94.36 %)
# THRESHOLD         1.700000
# NVOXEL_FOREGROUND 3609079
# NVOXEL_BACKGROUND 60390921
# VOLUME_FOREGROUND 6174331.225195
# VOLUME_BACKGROUND 103315430.127344
# #sumDensity 11394724.861821
# #0 Min -181.677473 Max 181.697001
# #1 Min -181.572385 Max 181.572318
# #2 Min -224.693745 Max 224.693814
# MOL_M        0.0097642155   -0.0000332227    0.0000346840
# MOL_CovM_X  8252.5630470519    0.0000114809    0.0021917452
# MOL_CovM_Y     0.0000114809 8242.1297166085    0.0366587065
# MOL_CovM_Z     0.0021917452    0.0366587065 12621.8236485691
# MOL_PCvar 12621.8236488770 8252.5630470509 8242.1297163017
# #ERROR(-ng 0):Cannot estimate GMM with Ngauss=0.
  dat = {}
  f = os.popen(command)
  for line in f:
    #print line
    if (len(line)>0) and (line.startswith('#')==0):
      fields = line.split()
      if (len(fields)>=2):
        key = fields[0]
        value = ''
        for i in range(1,len(fields)):
          if (value != ''):
            value += ' '
          value += fields[i]
        dat[key] = value 
      pass
  f.close()
  return(dat)
 

###############
#### MAIN #####
###############

OPT = {}
OPT['idir'] = '/home/takawaba/DB/emdb/structures'
OPT['ilist'] = ''

OPT['header_subdir'] = 'header'
OPT['map_subdir']    = 'map'
OPT['A'] = 'F'
OPT['div'] = '0/1'

OPT['O'] = 'G'


OPT['justcnt'] = 'F'
OPT['otab'] = 'out.csv'

if (len(sys.argv)<3):
  print "emdb_vol_pca.py <options>"
  print " to calculate volume and PCA for EMDB entries."
  print " using 'ContourLevel' of the emdb header XML file."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -idir   : input molecular file directory [%s]"%(OPT['idir'])
  print " -ilist  : input file for id list (optional) [%s]"%(OPT['ilist'])
  print " -header_subdir : subdir for 'header' (containing *.xml) under 'idir' [%s]"%(OPT['header_subdir'])
  print " -map_subdir : subdir for 'map'    (containing *.map) under 'idir' [%s]"%(OPT['map_subdir'])
  print " -A      : Action ('T' or 'F') [%s]"%(OPT['A'])
  print " -div    : Job division (bunshi)/(bunbo)[%s]"%(OPT['div'])
  print " -justcnt: stop after just counting number of entries for GMM conversion. (T or F)[%s]"%(OPT['justcnt'])
  print " -otab   : output tab-splited file [%s]"%(OPT['otab'])
  sys.exit(1)

read_option(sys.argv,OPT)
[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)

read_option(sys.argv,OPT)


#### (1) get acclist[] from OPT['ilist'] or OPT['idir'] ####

acclist = []
if (OPT['ilist']!=''):
  tmplist = []
  read_list_file(OPT['ilist'],tmplist)
  for id in (tmplist):
    if (id.startswith('EMD-')) or (id.startswith('emd-')):
      (head,acc) = id.split('-')
      acclist.append(acc)
    elif (id.startswith('EMD_')) or (id.startswith('emd_')):
      (head,acc) = id.split('_')
      acclist.append(acc)
else:
  dirlist = os.listdir(OPT['idir'])
  for dir in (dirlist):
    if (dir.startswith('EMD-')):
    ## dir = 'EMD-1022'
      (head,acc) = dir.split('-')
      acclist.append(acc)
  #for acc in (acclist):
  #  print acc
  #sys.exit(1)

## [2] read contourLevel from the XML file under OPT['header_subdir'] ##

print "## [1] read contourLevel from the XML file under OPT['header_subdir'] ##"
contourLevel = {}

for acc in (acclist):
  dirfull = OPT['idir'] + '/' + 'EMD-' + acc + '/' + OPT['header_subdir']
  
  if (os.path.isdir(dirfull)):
    dirlist_under = os.listdir(dirfull)
    for file in (dirlist_under):
      if (file.endswith('.xml')):
        file_full     = OPT['idir'] + '/' + 'EMD-' + acc +  '/' + OPT['header_subdir'] + '/' + file 
        #print file_full
        dat = {}
        read_XML_file(file_full,dat)
        acc = dat['acc']
        contourLevel[acc] = dat.get('contourLevel','0.0')
        ## for the case such as '2.48e+03',
        if (contourLevel[acc].find('e')!=-1):
          sys.stdout.write("%s-->"%(contourLevel[acc]))
          (head,tail) = contourLevel[acc].split('e')
          contourLevel[acc] = "%f"%(float(head) * math.pow(10.0,int(tail)))
          sys.stdout.write("%s\n"%(contourLevel[acc]))
  #        print "'%s' '%s' "%(dat['acc'],contourLevel[dat['acc']])

Nacc = len(acclist);
Nstart   = bunshi*int(Nacc/bunbo);
Nend     = (bunshi+1)*int(Nacc/bunbo);
if (bunshi>=(bunbo-1)):
  Nend = Nacc
print "#Nacclist %d bunshi/bunbo %d/%d start %d end %d"%(len(acclist),bunshi,bunbo,Nstart,Nend)


### [4] do 'gmconvert' for calculating volume and PCvar ###

if (OPT['A']=='T'):
  of = open(OPT['otab'],'w')
  print "#write_tabular_file() --> '%s'"%(OPT['otab'])
  of.write("id\tGRID_WIDTH\tGRID_SIZE_X\tGRID_SIZE_Y\tGRID_SIZE_Z\tthreshold\tvolume\tPCvar1\tPCvar2\tPCvar3\n");

for i in range(Nstart,Nend):
  acc = acclist[i]

  dirfull  = OPT['idir'] + '/' + 'EMD-' + acc + '/' + OPT['map_subdir']
  imapfile = OPT['idir'] + '/' + 'EMD-' + acc + '/' + OPT['map_subdir'] + '/' +  'emd_' + acc + '.map.gz'

  if (os.path.isfile(imapfile)):
    if (OPT['O']=='G'):
      command = "gmconvert -imap %s -ng 0 -ogmm temp.gmm -zth %s"%(imapfile,contourLevel.get(acc,'0.0'))

    print "#%s"%(command)
    if (OPT['A']=='T'):
      dat =  execute_gmconvert(command)
      #print dat
      if (dat.get('VOLUME_FOREGROUND','') != '') and (dat.get('MOL_PCvar','') != ''):
        pcvar = dat['MOL_PCvar'].split() 
        of.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(acc,dat.get('GRID_WIDHT',''),dat.get('GRID_SIZE_X',''),dat.get('GRID_SIZE_Y',''),dat.get('GRID_SIZE_Z',''),dat.get('THRESHOLD',''),dat.get('VOLUME_FOREGROUND',''),pcvar[0],pcvar[1],pcvar[2]))

if (OPT['A']=='T'):
  print "#write_tabular_file() --> '%s'"%(OPT['otab'])
  of.close()



now = datetime.now()
print "#START_DATE %s"%(OPT['START_DATE'])
print "#END_DATE   %s"%(now.strftime("%Y/%m/%d %H:%M:%S"))
print "#COMP_TIME  %s"%(now - OPT['start_datetime_now'])

