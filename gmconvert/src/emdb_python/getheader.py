#!/usr/bin/env python
##
## <getheader.py>
##


import sys
import os
import random
from datetime import datetime
from xml.etree.ElementTree import *


LastModDate = "Apr 5, 2013"

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


###############
#### MAIN #####
###############

OPT = {}
OPT['idir'] = '/DB/emdb/structures'
OPT['tail'] = '.xml'
OPT['subdir'] = 'header'
OPT['A'] = 'F'
if (len(sys.argv)<3):
  print "getheader.py <options>"
  print " for making tab-splited summary of the emdb header XML file."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -idir   : input molecular file directory [%s]"%(OPT['idir'])
  print " -subdir : subdir under 'idir' [%s]"%(OPT['subdir'])
  print " -tail   : tail of map file [%s]"%(OPT['tail'])
  print " -A      : Action ('T' or 'F') [%s]"%(OPT['A'])
  sys.exit(1)

read_option(sys.argv,OPT)

dirlist = os.listdir(OPT['idir'])
print "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%('acc','resolutionByAuthor','numColumns','numRows','numSections','fittedPDBEntryId','pdbEntryId','keywords','entry','title')
for dir in (dirlist):
  dirfull = OPT['idir'] + '/' + dir + '/' + OPT['subdir']
  if (os.path.isdir(dirfull)): 
    dirlist_under = os.listdir(dirfull)
    for file in (dirlist_under):
      if (file.endswith(OPT['tail'])):
        file_full     = OPT['idir'] + '/' + dir + '/' + OPT['subdir'] + '/' + file 
        #print file_full
        dat = {}
        read_XML_file(file_full,dat)
        print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(dat['acc'],dat.get('resolutionByAuthor',''),dat.get('numColumns',''),dat.get('numRows',''),dat.get('numSections',''),dat.get('contourLevel','?'),dat.get('fittedPDBEntryId',''),dat.get('pdbEntryId',''),dat.get('keywords',''),dat.get('entry',''),dat.get('title',''))
