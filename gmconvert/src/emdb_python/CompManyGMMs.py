#!/usr/bin/env python
##
## <CompManyGMMs.py>
##


import sys
import os
import random
from datetime import datetime
from xml.etree.ElementTree import *

LastModDate = "July 27, 2014"

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



def do_gmfit(cg,sg1,I='R',NI='10',NS='10',Action='F'):
  command = "gmfit -cg %s -sg1 %s -KV T -I %s -NI %s -NS %s"%(cg,sg1,I,NI,NS)
  print "#command:%s"%(command) 
  Etotal = '0.0'
  CorrCoeff = '0.0'
  Roverlap = '0.0'
  if (Action=='T'):
    f = os.popen(command)
    for line in f:
  # Etotal_GBEST1      -5.401453e-07
  # EcmpfitGG_GBEST1   -5.401453e-07
  # ErepulsGG_GBEST1   0.000000e+00
  # CorrCoeffGG_GBEST1 0.861865
  # RoverlapGG_GBEST1  0.891127
  # RMSDorigGG_GBEST1  13.770008

      if (line.startswith("Etotal_GBEST1")):
        line = line.rstrip('\n')
        (key,value) = line.split()             
        Etotal = value
      
      if (line.startswith("CorrCoeffGG_GBEST1")):
        line = line.rstrip('\n')
        (key,value) = line.split()             
        CorrCoeff = value

      if (line.startswith("RoverlapGG_GBEST1")):
        line = line.rstrip('\n')
        (key,value) = line.split()             
        Roverlap = value

    f.close()
  return(Etotal,CorrCoeff,Roverlap,command) 


###############
#### MAIN #####
###############

OPT = {}
OPT['M'] = 'L'
OPT['lgdir'] = 'tmpout'
OPT['tail'] = '.gmm'
OPT['A'] = 'F'
OPT['div'] = '0/1'
OPT['Q'] = ''
OPT['NI'] = '100'
OPT['NS'] = '10'
OPT['I'] = 'P'
OPT['of'] = '-'
OPT['ilist'] = ''

now = datetime.now()
OPT['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")

if (len(sys.argv)<3):
  print "CompManyGMMs.py <options>"
  print " for making GMM from the EMDB map."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -M     : Mode. 'L':one-vs-library 'A':all-vs-all"
  print " -ilist : input list file for library GMM [%s]"%(OPT['ilist'])
  print " -lgdir : input directory for library GMM [%s]"%(OPT['lgdir'])
  print " -tail : tail of map file [%s]"%(OPT['tail'])
  print " -A    : Action ('T' or 'F') [%s]"%(OPT['A'])
  print " -div  : Job division (bunshi)/(bunbo)[%s]"%(OPT['div']) 
  print " -of   : output search result file [%s]"%(OPT['of']) 
  print "<options only for '-M L'>"
  print " -Q   : input Query GMM file [%s]"%(OPT['Q'])
  print "<options for gmfit>"
  print " -NI  : number of initial configurations. [%s]"%(OPT['NI'])
  print " -I   : Initial Configuration. 'R':random, 'P':PCA, 'Y':symmetric,'F':segment-fitting"
  print "      : 'y:symmetric segment fitting, 'L':local transformation to the original config,'O':keep original config.[%s]"%(OPT['I'])
  print " -NS  : number of configurations for local search [%s]"%(OPT['NS'])
  sys.exit(1)

read_option(sys.argv,OPT)
[bunshi,bunbo] = OPT['div'].split('/')
bunshi = int(bunshi)
bunbo  = int(bunbo)

lib_gmmfile_list = []
lib_gmmfile_property = {}

### (1) read list file ###
if (OPT['ilist'] != ''):
  idlist = []
  read_list_file(OPT['ilist'],lib_gmmfile_list,lib_gmmfile_property)
else:
  ### or, scan the directory 'idir' ###
  dirlist = os.listdir(OPT['lgdir'])
  for file in (dirlist):
    if (file.endswith(OPT['tail'])):
      field = file.split('.') 
      lib_gmmfile_list.append(file)
    print lib_gmmfile_list


########################################
### MODE 'L' : one-vs-library search ###
########################################

if (OPT['M']=='L'):
  Nfile = len(lib_gmmfile_list);
  Nstart   = bunshi*int(Nfile/bunbo)
  Nend     = (bunshi+1)*int(Nfile/bunbo)
  if (bunshi>=(bunbo-1)):
    Nend = Nfile
  print "#Nfile %d bunshi/bunbo %d/%d start %d end %d"%(Nfile,bunshi,bunbo,Nstart,Nend)
  
  ### (2) do gmfit ###
  Etotal = {}
  CorrCoeff = {}
  Roverlap = {}
  command_example = ''
  for i in range(Nstart,Nend):
    libfile = lib_gmmfile_list[i]
    libfile_full     = OPT['lgdir'] + '/' + libfile
    if (libfile_full.endswith('.gmm')==0):
      libfile_full += '.gmm'

    (etotal,cc,roverlap,command_str) = do_gmfit(OPT['Q'],libfile_full,I=OPT['I'],NI=OPT['NI'],NS=OPT['NS'],Action=OPT['A'])
    Etotal[libfile]    = etotal
    CorrCoeff[libfile] = cc 
    Roverlap[libfile]  = roverlap 
    if (command_example==''):
      command_example = command_str 
 
  #### (3) Output  filelst sorted by CorrCoeff ####
  now = datetime.now()
  OPT['END_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")
  
  if (OPT['of']=='-'):
    of = sys.stdout
  else:
    of = open(OPT['of'],'w')
    print "#write_library_search_result() --> '%s'"%(OPT['of']) 
 
  of.write("#COMMAND    %s\n"%(OPT['COMMAND']))
  of.write("#START_DATE %s\n"%(OPT['START_DATE']))
  of.write("#END_DATE   %s\n"%(OPT['END_DATE']))
  of.write("#COMMAND_EXAMPLE   %s\n"%(command_example))
  of.write("#Nfile %d\n"%(Nfile))
  
  sorted_list  = sorted(lib_gmmfile_list, lambda x,y:cmp(float(CorrCoeff.get(y,'0.0')),float(CorrCoeff.get(x,'0.0'))))
  rank = 0
  
  of.write("#[rank] [gmmfile] [CorrCoeff] [Roverlap] [property]\n")
  for file in (sorted_list):
    rank += 1
    of.write("%-4d %s %s %s %s\n"%(rank,file,CorrCoeff.get(file,'0.0'),Roverlap.get(file,'0.0'),lib_gmmfile_property.get(file,'x')))
  
  if (OPT['of']!='-'):
    of.close()


########################################
### MODE 'A' : all-vs-all comparison ###
########################################
if (OPT['M']=='A'):
  Nfile = len(lib_gmmfile_list);
  Npair = Nfile * Nfile
  Nstart   = bunshi*int(Npair/bunbo)
  Nend     = (bunshi+1)*int(Npair/bunbo)
  if (bunshi==0):
    Nstart = 0
  if (bunshi==(bunbo-1)):
    Nend = Npair 

  sys.stdout.write("#div %d/%d  Nstart %d Nend %d Npair %d\n"%(bunshi,bunbo,Nstart,Nend,Npair))
  if (OPT['of']=='-'):
    of = sys.stdout
  else:
    of = open(OPT['of'],'w')
    print "#write_all_vs_all_result() --> '%s'"%(OPT['of']) 

  of.write("#COMMAND    %s\n"%(OPT['COMMAND']))
  of.write("#START_DATE %s\n"%(OPT['START_DATE']))
  of.write("#Nfile %d\n"%(Nfile))
  of.write("#div %d/%d  Nstart %d Nend %d Npair %d\n"%(bunshi,bunbo,Nstart,Nend,Npair))
#  of.write("#[gmmfileA] [gmmfileB] [propertyA] [propertyB] [CC] [Roverlap]\n")
  of.write("#[gmmfileA] [gmmfileB] [CC] [Roverlap] [propertyA] [propertyB]\n")
  npair = 0
  command_example = ''
  for a in range(Nfile):

    libfileA = lib_gmmfile_list[a]
    libfile_fullA = OPT['lgdir'] + '/' + libfileA
    if (libfile_fullA.endswith('.gmm')==0):
      libfile_fullA += '.gmm'

    for b in range(Nfile):
      libfileB = lib_gmmfile_list[b]
      libfile_fullB = OPT['lgdir'] + '/' + libfileB
      if (libfile_fullB.endswith('.gmm')==0):
        libfile_fullB += '.gmm'
      if (Nstart<=npair) and (npair <Nend):
        (etotal,cc,roverlap,command_str) = do_gmfit(libfile_fullA,libfile_fullB,I=OPT['I'],NI=OPT['NI'],NS=OPT['NS'],Action=OPT['A'])
        #of.write("%s %s %2s %2s %s %s\n"%(libfileA,libfileB,lib_gmmfile_property[libfileA],lib_gmmfile_property[libfileB],cc,roverlap))
        of.write("%s %s %s %s %s %s\n"%(libfileA,libfileB,cc,roverlap,lib_gmmfile_property[libfileA],lib_gmmfile_property[libfileB]))
        if (command_example == ''):
          command_example = command_str
      npair += 1


  of.write("#COMMAND_EXAMPLE   %s\n"%(command_example))
  now = datetime.now()
  OPT['END_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")
  of.write("#END_DATE   %s\n"%(OPT['END_DATE']))
  if (OPT['of']!='-'):
    print "#write_all_vs_all_result() --> '%s'"%(OPT['of']) 
    of.close()
