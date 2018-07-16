#!/usr/bin/env python
##
## <sasbdb2gmm.py>
##


import sys
import os
import random
from datetime import datetime
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



def read_header_of_cif_file(ifname):
#>>EXAMPLES OF SASBDB CIF FILES <<
# data_SASDAB5_MODEL_101
# #
# _sas_model.id                  101
# _sas_model.type_of_model       dummy
# _sas_model.software            dammif
# _sas_model.version             .
# _sas_model.symmetry            .
# _sas_model.comment             .
# _sas_model.fitting_id          93
# _sas_model.radius              1.90
# #
# loop_
# _atom_site.group_PDB
# _atom_site.id
# _atom_site.type_symbol
# _atom_site.label_atom_id
# _atom_site.label_alt_id
# _atom_site.label_comp_id
#-----------------------------------------
# data_SASDAG2_MODEL_18
# #
# _sas_model.id                  18
# _sas_model.type_of_model       dummy
# _sas_model.software            gasbor
# _sas_model.version             .
# _sas_model.symmetry            P1
# _sas_model.comment             "ran it 10 times and took the model with the lowest NSD (damsel) and visually best fit (pymol)"
# _sas_model.fitting_id          10
# _sas_model.radius              1.90
#
#-----------------------------------------
# data_SASDA96_MODEL_168
# #
# _sas_model.id                  168
# _sas_model.type_of_model       atomic
# _sas_model.software            CRYSOL
# _sas_model.version             .
# _sas_model.symmetry            .
# _sas_model.comment             .
# _sas_model.fitting_id          142
# _sas_model.radius              .
# #
#-----------------------------------------
# data_SASDA94_MODEL_50
# #
# _sas_model.id                  50
# _sas_model.type_of_model       mix
# _sas_model.software            Bunch
# _sas_model.version             .
# _sas_model.symmetry            P2
# _sas_model.comment             .
# _sas_model.fitting_id          44
# _sas_model.radius              1.90
#
  dat = {}
  if not os.access(ifname,os.R_OK):
    print "#ERROR:Can't open CIF file '%s'" %(ifname)
    return(0)
  f = open(ifname)
  for line in f:
    line = line.rstrip('\n')
    if (line.startswith('data_')):
# data_SASDAB5_MODEL_101
# data_SASDAG2_MODEL_18
# data_SASDA94_MODEL_50
      fields = line.split('_')
      dat['sasbdb_id'] = field[1] 

    if (line.startswith('_sas_model.')):
      fields = line.split()
      #print fields
      (head,key) = fields[0].split('.')
      dat[key] = fields[1]
  f.close()
  return(dat)





###############
#### MAIN #####
###############

OPT = {}
OPT['idir'] = '/home/takawaba/DB/emdb/structures'
OPT['icifdir'] = '/home/takawaba/etc/sas/splitcif'
OPT['ilist'] = ''

OPT['A'] = 'F'
OPT['div'] = '0/1'
OPT['ng'] = '20'
OPT['emalg'] = 'G'
OPT['delzw'] = 'T'
OPT['delid'] = 'T'

OPT['O'] = 'G'

OPT['ogdir']  = '/home/takawaba/etc/sas/gmm'

OPT['update'] = 'F'
OPT['avex'] = 'F'
OPT['ogdirup'] = ''
OPT['justcnt'] = 'F'

if (len(sys.argv)<3):
  print "sasbdb2gmm.py <options>"
  print " to convert EMDB maps into GMM files."
  print " using 'ContourLevel' of the emdb header XML file."
  print " code by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -icifdir  : input dir for cif file [%s]"%(OPT['icifdir'])
  print " -idir   : input molecular file directory [%s]"%(OPT['idir'])
  print " -ilist  : input file for id list (optional) [%s]"%(OPT['ilist'])
  print " -ogdir  : output GMM directory [%s]"%(OPT['ogdir'])
  print " -ng     : Ngauss [%s]"%(OPT['ng'])
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


#### (1) get model_id_list[] from OPT['ilist'] or OPT['idir'] ####

model_id_list = []
if (OPT['ilist']!=''):
  tmplist = []
  read_list_file(OPT['ilist'],tmplist)
  for id in (tmplist):
    if (id.startswith('EMD-')) or (id.startswith('emd-')):
      (head,acc) = id.split('-')
      model_id_list.append(acc)
    elif (id.startswith('EMD_')) or (id.startswith('emd_')):
      (head,acc) = id.split('_')
      model_id_list.append(acc)
else:
  fnamelist = os.listdir(OPT['icifdir'])
  for fname in (fnamelist):
    if (fname.endswith('.cif')):
      (head,tail) = fname.split('.')
      model_id_list.append(head)
for model_id in (model_id_list):
  print model_id 


### [2] For OPT['update']='T'  ###
print "#Nmodel_id_all %d"%(len(model_id_list))

if (OPT['update']=='T'):
  print "#update '%s'"%(OPT['update'])
  model_id_list_orig = []
  for acc in (model_id_list):
    model_id_list_orig.append(acc)
  model_id_list = []

  for acc in (model_id_list_orig):
    imapfile = OPT['idir'] + '/' + 'EMD-' + acc + '/' +  'emd_' + acc + '.map.gz'
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
      model_id_list.append(acc)
  pass

print "#Nacc %d"%(len(model_id_list))
 
if (OPT['justcnt']=='T'):
  print "#STOP after counting entries to be updated."
  sys.exit(1)

#for id in (model_id_list):
#  if (exists_gmm.get(id,0)==1):
#    print "%s  exists_gmm"%(id)
#  else:
#    print "%s  non_exists_gmm"%(id)
#sys.exit(1)

Nmodel_id = len(model_id_list);
Nstart   = bunshi*int(Nmodel_id/bunbo);
Nend     = (bunshi+1)*int(Nmodel_id/bunbo);
if (bunshi>=(bunbo-1)):
  Nend = Nmodel_id
print "#Nmodel_id %d bunshi/bunbo %d/%d start %d end %d"%(Nmodel_id,bunshi,bunbo,Nstart,Nend)



### [4] do 'gmconvert' ###

for i in range(Nstart,Nend):
  model_id = model_id_list[i]

  iciffile = OPT['icifdir'] + '/' + model_id + '.cif'
  ogmmfile = OPT['ogdir']   + '/' + model_id + '.gmm'



  if (os.path.isfile(iciffile)):
    CifInfo = read_header_of_cif_file(iciffile)
    command = "gmconvert -icif %s -ng %s -emalg %s -delzw %s delid %s -ogmm %s"%(iciffile,OPT['ng'],OPT['emalg'],OPT['delzw'],OPT['delid'],ogmmfile)
    if (CifInfo.get('type_of_model','.') == 'atomic') or (CifInfo.get('type_of_model','.') == 'mix'):
      command += " -atmrw A"
    elif (CifInfo.get('radius','.') != '.') or (CifInfo.get('type_of_model','.') == 'mix'):
      command += " -atmrw U -raduni %s"%(CifInfo['radius'])
    else:
      command += " -atmrw U -raduni 1.9"
    print "#%s"%(command)
    if (OPT['A']=='T'):
      os.system(command)

    if (OPT['ogdirup'] != ''):
      ogmmfileup =  "%s/%s.gmm"%(OPT['ogdirup'],model_id)
      upcommand = "cp %s %s"%(ogmmfile,ogmmfileup)
      print "#%s"%(upcommand)
      if (OPT['A']=='T'):
        os.system(upcommand)
    


