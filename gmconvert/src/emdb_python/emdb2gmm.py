#!/usr/bin/env python
##
## <emdb2gmm.py>
##


import sys
import os
import random
from datetime import datetime
from xml.etree.ElementTree import *
import math

LastModDate = "Feb 16, 2016"

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



def read_XML_file(ifname,dat):
  #print "#read_XML_file('%s')"%(ifname)
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
OPT['ng'] = '10'
OPT['maxsize'] = '64'
OPT['emalg'] = 'G'
OPT['delzw'] = 'T'
OPT['delid'] = 'T'


OPT['ogdir'] = '/DB/emdb/GMM/GMM_NG20'

OPT['update'] = 'F'
OPT['avex'] = 'F'
OPT['ogdirup'] = ''
OPT['justcnt'] = 'F'

if (len(sys.argv)<3):
  print "emdb2gmm.py <options>"
  print " to convert EMDB maps into GMM files."
  print " using 'ContourLevel' of the emdb header XML file."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -idir   : input molecular file directory [%s]"%(OPT['idir'])
  print " -ilist  : input file for id list (optional) [%s]"%(OPT['ilist'])
  print " -ogdir  : output GMM directory [%s]"%(OPT['ogdir'])
  print " -header_subdir : subdir for 'header' (containing *.xml) under 'idir' [%s]"%(OPT['header_subdir'])
  print " -map_subdir : subdir for 'map'    (containing *.map) under 'idir' [%s]"%(OPT['map_subdir'])
  print " -ng     : Ngauss [%s]"%(OPT['ng'])
  print " -maxsize: maximum voxel size of each axis(if over, reducing size).[%s]"%(OPT['maxsize'])
  print " -emalg  : type for EM algorithm. 'P'oint_observed_EM,'G'mm_observe_EM [%s]"%(OPT['emalg'])
  print " -delzw : Delete Zero-weight gdfs from the GMM. ('T' or 'F') [%s]"%(OPT['delzw'])
  print " -delid : Delete identical gdfs in the GMM. ('T' or 'F') [%s]"%(OPT['delid'])
  print " -A      : Action ('T' or 'F') [%s]"%(OPT['A'])
  print " -div    : Job division (bunshi)/(bunbo)[%s]"%(OPT['div'])
  print " -justcnt: stop after just counting number of entries for GMM conversion. (T or F)[%s]"%(OPT['justcnt'])
  print "<options for update>"
  print " -update : make gmm only updating files (T or F) [%s]"%(OPT['update'])
  print " -ogdirup: output dir only for update GMMs. Copying from the '-ogdir'. [%s]"%(OPT['ogdirup'])
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


### [2] For OPT['update']='T'  ###
print "#Nacc_all %d"%(len(acclist))

if (OPT['update']=='T'):
  print "#update '%s'"%(OPT['update'])
  acclist_orig = []
  for acc in (acclist):
    acclist_orig.append(acc)
  acclist = []

  for acc in (acclist_orig):
    imapfile = OPT['idir'] + '/' + 'EMD-' + acc + '/' + OPT['map_subdir'] + '/' +  'emd_' + acc + '.map.gz'
    ogmmfile =  "%s/emd_%s.gmm"%(OPT['ogdir'],acc)
    ogmmfile_date = 0
    ogmmfile_dt = ''
    imapfile_date = 0
    imapfile_dt = '' 
    if (os.path.isfile(imapfile)):
      imapfile_date = os.stat(imapfile).st_mtime
      imapfile_dt = datetime.fromtimestamp(os.stat(imapfile).st_mtime).strftime('%Y-%m-%d')
    if (os.path.isfile(ogmmfile)):
      ogmmfile_date = os.stat(ogmmfile).st_mtime
      ogmmfile_dt = datetime.fromtimestamp(os.stat(ogmmfile).st_mtime).strftime('%Y-%m-%d')
    diff_date = imapfile_date - ogmmfile_date
    if (diff_date>0):
      print "%s map_date %s gmm_date %s diff_date %f map %s gmm %s"%(acc,imapfile_date,ogmmfile_date,imapfile_date-ogmmfile_date, imapfile_dt, ogmmfile_dt)
      acclist.append(acc)
  pass

print "#Nacc %d"%(len(acclist))
 
if (OPT['justcnt']=='T'):
  print "#STOP after counting entries to be updated."
  sys.exit(1)

#for id in (acclist):
#  if (exists_gmm.get(id,0)==1):
#    print "%s  exists_gmm"%(id)
#  else:
#    print "%s  non_exists_gmm"%(id)
#sys.exit(1)

## [3] read contourLevel from the XML file under OPT['header_subdir'] ##

print "## [1] read contourLevel from the XML file under OPT['header_subdir'] ##"
contourLevel = {}

for acc in (acclist):
  dirfull = OPT['idir'] + '/' + 'EMD-' + acc + '/' + OPT['header_subdir']
  
  if (os.path.isdir(dirfull)):
    dirlist_under = os.listdir(dirfull)
    for file in (dirlist_under):
      #if (file.endswith('.xml')):
      if (file == 'emd-' + acc + '.xml'):
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





### [4] do 'gmconvert' ###

for i in range(Nstart,Nend):
  acc = acclist[i]

  dirfull = OPT['idir'] + '/' + 'EMD-' + acc + '/' + OPT['map_subdir']
  imapfile = OPT['idir'] + '/' + 'EMD-' + acc + '/' + OPT['map_subdir'] + '/' +  'emd_' + acc + '.map.gz'
  ogmmfile =  "%s/emd_%s.gmm"%(OPT['ogdir'],acc)

  if (os.path.isfile(imapfile)):
    command = "gmconvert -imap %s -ng %s -emalg %s -delzw %s delid %s -ogmm %s -zth %s"%(imapfile,OPT['ng'],OPT['emalg'],OPT['delzw'],OPT['delid'],ogmmfile,contourLevel.get(acc,'0.0'))
    if (int(OPT['maxsize'])>0):
      command += " -maxsize %s"%(OPT['maxsize'])

    print "#%s"%(command)
    if (OPT['A']=='T'):
      os.system(command)

    if (OPT['ogdirup'] != ''):
      ogmmfileup =  "%s/emd_%s.gmm"%(OPT['ogdirup'],acc)
      upcommand = "cp %s %s"%(ogmmfile,ogmmfileup)
      print "#%s"%(upcommand)
      if (OPT['A']=='T'):
        os.system(upcommand)
    


