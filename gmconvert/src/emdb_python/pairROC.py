#!/usr/bin/env python

import os
import sys
from datetime import datetime


LastModDate = "Aug 1, 2014"

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




def read_pair_score_file(ifname,scdat,colA='1',colB='2',colS='5'):
  print "#read_pair_score_file('%s')"%(ifname)
  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR:Can't open listfile '%s'."%(ifname)
    sys.exit(1)

  icolA = int(colA)-1
  icolB = int(colB)-1
  icolS = int(colS)-1
  f = open(ifname)
  for line in f:
    if (line.startswith('#')==0) and (len(line)>10):
      line = line.rstrip('\n')
      field = line.split()
      Nfield = len(field)
      if (0<=icolA) and (icolA<Nfield)and (0<=icolB) and (icolB<Nfield) and (0<=icolS) and (icolS<Nfield):
        idA    = field[int(colA)-1]
        idB    = field[int(colB)-1]
        score  = float(field[int(colS)-1])
        index = idA + ':' + idB
        scdat[index] = score

###############
#### MAIN #####
###############

OPT = {}
OPT['ilist'] = ''
OPT['ipair'] = ''
OPT['cA'] = '1'
OPT['cB'] = '2'
OPT['cS'] = '3'
OPT['of'] = '-'
OPT['sctype'] = 'S'
OPT['exself'] = 'T'

if (len(sys.argv)<3):
  print "pairROC.py <options>"
  print " for making ROC curve from pairwise comparisons."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -ilist : input list file [%s]"%(OPT['ilist'])
  print " -ipair : input pair score file [%s]"%(OPT['ipair'])
  print " -cA    : column for idA   (1,2,3....) [%s]"%(OPT['cA'])
  print " -cB    : column for idB   (1,2,3....) [%s]"%(OPT['cB'])
  print " -cS    : column for Score (1,2,3....) [%s]"%(OPT['cS'])
  print " -sctype: score type. 'S'imilarity, 'D'istance [%s]"%(OPT['sctype'])
  print " -of    : output file [%s]"%(OPT['of'])
  print " -exself: excluding self matches ('T' or 'F') [%s]"%(OPT['exself'])
  sys.exit(1)

read_option(sys.argv,OPT)

### [1] read pair score file ###
ScDat = {}
read_pair_score_file(OPT['ipair'],ScDat,colA=OPT['cA'],colB=OPT['cB'],colS=OPT['cS'])
#if (OPT['sctype']=='D'):
#  pairlist = sorted(ScDat.keys(), lambda x,y:cmp(ScDat[x],ScDat[y]))
#else:
#  pairlist = sorted(ScDat.keys(), lambda x,y:cmp(ScDat[y],ScDat[x]))

### [2] read list file and make pairlise[] ###
Property = {}
id_list = []
read_list_file(OPT['ilist'],id_list,Property)

tmp_pair_list = []
for i in range(len(id_list)):
  for j in range(i,len(id_list)):
    if (i!=j) or ((i==j) and (OPT['exself'] != 'T')):
      pair_index = "%s:%s"%(id_list[i],id_list[j])
      tmp_pair_list.append(pair_index)

if (OPT['sctype']=='D'):
  pairlist = sorted(tmp_pair_list, lambda x,y:cmp(ScDat[x],ScDat[y]))
else:
  pairlist = sorted(tmp_pair_list, lambda x,y:cmp(ScDat[y],ScDat[x]))

tmp_pair_list = []



### [3] Count No[0],No[1] ###

No = [0,0]
for pair in (pairlist):
  (idA,idB) = pair.split(':')
  if (idA in Property) and (idB in Property) and (Property[idA] ==Property[idB]):
    observe = 1
  else: 
    observe = 0
  No[observe] += 1
 
### [4] Caculate recall,precision,FPR,TPR ###
if (OPT['of']=='-'):
  of = sys.stdout
else:
  of = open(OPT['of'],'w') 

of.write("#COMMAND %s\n"%(OPT['COMMAND']))
of.write("#DATE    %s\n"%(OPT['START_DATE']))
of.write("#sctype  %s\n"%(OPT['sctype']))
of.write("#exself  %s\n"%(OPT['exself']))
of.write("### Nop[O][P] : number of observation O and prediction P\n")
of.write("#recall(=sensitivity=TPR)         = Nop[1][1]/(Nop[1][0] + Nop[1][1])\n")
of.write("#precision                        = Nop[1][1]/(Nop[0][1] + Nop[1][1])\n")
of.write("#specificity(recall of the class0)= Nop[0][0]/(Nop[0][0] + Nop[0][1])\n")
of.write("#FPR (1-specificity)              = Nop[0][1]/(Nop[0][0] + Nop[0][1])\n")

Nop = [[0,0],[0,0]]

Nop[0][0] = No[0] 
Nop[1][0] = No[1] 
Nop[0][1] = 0
Nop[1][1] = 0
of.write("#[idA:1] [idB:2] [Score:3] [FPR:4] [TPR(recall):5] [precision:6] [classA:7] [classB:8]\n")

AUC_FPR_TPR = 0.0
FPR0 = -1.0
TPR0 = -1.0


for pair in (pairlist):
  (idA,idB) = pair.split(':')
  if (idA in Property) and (idB in Property) and (Property[idA] ==Property[idB]):
    observe = 1
  else: 
    observe = 0

  Nop[observe][1] += 1
  Nop[observe][0] -= 1
  recall    = float(Nop[1][1])/float(Nop[1][0] + Nop[1][1])
  precision = float(Nop[1][1])/float(Nop[0][1] + Nop[1][1])
  ## sensitivity == recall
  sensitivity = float(Nop[1][1])/float(Nop[1][0] + Nop[1][1])
  ## specificity == recall for the class 0 
  specificity = float(Nop[0][0])/float(Nop[0][0] + Nop[0][1])
  FPR = 1-specificity
  TPR = sensitivity
  of.write("%s %s %s %f %f %f %s %s\n"%(idA,idB,ScDat[pair],FPR,recall,precision,Property.get(idA,'---'),Property.get(idB,'---')))
  if (FPR0>=0.0):
    AUC_FPR_TPR += (TPR + TPR0)*(FPR-FPR0)/2.0


  FPR0 = FPR
  TPR0 = TPR

of.write("#%s AUC_FPR_TPR  %f\n"%(OPT['ipair'],AUC_FPR_TPR))
sys.stdout.write("#%s AUC_FPR_TPR  %f\n"%(OPT['ipair'],AUC_FPR_TPR))


if (OPT['of']!='-'):
  of.close()
