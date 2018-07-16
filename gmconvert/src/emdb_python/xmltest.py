#!/usr/bin/env python

import sys
import os
from xml.etree.ElementTree import *
from xml.etree import ElementTree 

ifname = sys.argv[1]
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
for e in elem.getiterator():
  print e.tag
  if (e.tag.startswith('title')):
    #print ">%s :: '%s'"%(e.tag,e.text)
    title = e.text.rstrip('\n')
  if (e.tag == 'pdbEntryId'):
    print ">%s :: '%s'"%(e.tag,e.text)
    id = e.text
    id = id.rstrip('\n')
    if (len(id)>=4):
      if (pdbEntryId != ''):
        pdbEntryId += ':' 
      pdbEntryId += id

print "%s\t%s\t%s"%(acc,pdbEntryId,title)
 
