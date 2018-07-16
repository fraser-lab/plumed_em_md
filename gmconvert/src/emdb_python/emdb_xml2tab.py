#!/usr/bin/env python
##
## <emdb_xml2tab.py>
##


import sys
import os
import random
from datetime import datetime
from xml.etree.ElementTree import *

LastModDate = "Aug 9, 2015"

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
  root = fromstring(contents)
  acc = root.get('accessCode').rstrip('\n')
  #print "#acc '%s'"%(acc)
  pdbEntryId = ''
  fittedPDBEntryId = ''
  for c in root:
    #print "CHILD %s"%(c.tag)
    if (c.tag == 'map'):  
      for g in c:
        if (g.tag == 'contourLevel'):
          dat['contourLevel'] = g.text.rstrip('\n')

        for h in g:
          if (h.tag == 'numColumns'):
            dat['numColumns'] = h.text.rstrip('\n')
          if (h.tag == 'numRows'):
            dat['numRows'] = h.text.rstrip('\n')
          if (h.tag == 'numSections'):
            dat['numSections'] = h.text.rstrip('\n')
          if (h.tag == 'pixelX'):
            dat['pixelX'] = h.text.rstrip('\n')
          if (h.tag == 'pixelY'):
            dat['pixelY'] = h.text.rstrip('\n')
          if (h.tag == 'pixelZ'):
            dat['pixelZ'] = h.text.rstrip('\n')
          if (h.tag == 'cellAs'):
            dat['cellA'] = h.text.rstrip('\n')
          if (h.tag == 'cellB'):
            dat['cellB'] = h.text.rstrip('\n')
          if (h.tag == 'cellC'):
            dat['cellC'] = h.text.rstrip('\n')

    if (c.tag == 'sample'):  
      for g in c:
        if (g.tag == 'name'):
          dat['sample'] = g.text.rstrip('\n')
          dat['sample'] = dat['sample'].replace('\n', ' ')
          dat['sample'] = dat['sample'].replace('\r', ' ')
          dat['sample'] = dat['sample'].replace('\n', ' ')
          dat['sample'] = dat['sample'].replace('\r', ' ')
          #print " sample '%s'"%(dat['sample'])
    if (c.tag == 'processing'):  
      for g in c:
        for h in g:
          if (h.tag == 'resolutionByAuthor'):
            dat['resolutionByAuthor'] = h.text.rstrip('\n')

    if (c.tag == 'deposition'):  
      for g in c:
        if (g.tag == 'keywords'):
          dat['keywords'] = g.text

        for h in g:
          if (h.tag == 'fittedPDBEntryId'):
            id = h.text
            id = id.rstrip('\n')
            if (fittedPDBEntryId != ''):
              fittedPDBEntryId += ':'
            fittedPDBEntryId += id

    if (c.tag == 'experiment'):  
      for g in c:
        if (g.tag == 'fitting'):  
          for h in g:
            if (h.tag == 'pdbEntryIdList'):
              for i in h:
                if (g.tag == 'pdbEntryId'):
                  id = i.text
                  id = id.rstrip('\n')
                  if (pdbEntryId != ''):
                    pdbEntryId += ':'
                  pdbEntryId += id


    for g in c:
      #print ">%s :: '%s'"%(c.tag,g.tag)
      if (g.tag.startswith('title')):
        #print ">%s :: '%s'"%(e.tag,e.text)
        string       = g.text.replace('\n',' ')
        dat['title'] = string.replace('\t',' ')

      if (g.tag == 'authors'):
        dat['authors'] = g.text.rstrip('\n')
        dat['entry'] = dat.get('entry','') + ':' + g.text.rstrip('\n')


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
OPT['ixml'] = ''
OPT['idir'] = '/home/takawaba/DB/emdb/structures'
OPT['tail'] = '.xml'
OPT['subdir'] = 'header'
OPT['otab'] = 'out.tab'
if (len(sys.argv)<3):
  print "emdb_xml2tab.py <options>"
  print " for making tab-splited summary of the emdb header XML file."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -idir   : input molecular file directory [%s]"%(OPT['idir'])
  print " -subdir : subdir under 'idir' [%s]"%(OPT['subdir'])
  print " -tail   : tail of map file [%s]"%(OPT['tail'])
  print " -otab   : output tab-splitted file  [%s]"%(OPT['otab'])
  print " -ixml   : input xml file (optional. mainly for testing) [%s]"%(OPT['ixml'])
  sys.exit(1)

read_option(sys.argv,OPT)

dirlist = os.listdir(OPT['idir'])
keylist = ['acc','resolutionByAuthor','numColumns','numRows','numSections','cellA', 'cellB', 'cellC', 'pixelX','pixelY','pixelZ','fittedPDBEntryId','pdbEntryId','keywords','entry','title','sample']

Nid_out = 0
of = open(OPT['otab'],'w')
print "#write_tab_file() --> '%s'"%(OPT['otab'])
#of.write("#")
for i in range(len(keylist)):
  of.write("%s\t"%(keylist[i]))
of.write("\n")

if (OPT['ixml'] != ''):
  dat = {}
  read_XML_file(OPT['ixml'],dat)
  for i in range(len(keylist)):
    of.write("%s:%d\t"%(keylist[i],i+1))
  of.write("\n")

  dat = {}
  read_XML_file_old(OPT['ixml'],dat)
  for key in (keylist):
    of.write("%s\t"%(dat.get(key,'')))
  of.write("\n")

  of.close()
  sys.exit(1)

if (OPT['idir'] != ''):
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
          for key in (keylist):
            of.write("%s\t"%(dat.get(key,'')))
          of.write("\n")
          Nid_out += 1 

#        print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(dat['acc'],dat.get('resolutionByAuthor',''),dat.get('numColumns',''),dat.get('numRows',''),dat.get('numSections',''),dat.get('contourLevel','?'),dat.get('fittedPDBEntryId',''),dat.get('pdbEntryId',''),dat.get('keywords',''),dat.get('entry',''),dat.get('title',''))
  print "#write_tab_file(Nid_out:%d) --> '%s'"%(Nid_out,OPT['otab'])
  of.close()
