
## <mmCIF.py>                       ##
## functions for reading mmCIF file ##

import sys
import os
import math
import re
from datetime import datetime

LastModDate = "Dec 17, 2014"

def read_option(argv,opt_dic):
  now = datetime.now()
  opt_dic['START_DATE'] = now.strftime("%Y/%m/%d %H:%M:%S")
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]


def read_mmCIF_file(ifname,dat,focus_category=[]):
#>> STRUCTURE OF VARIABLE dat{}
#   --> dat is "list of double dictionary".
#
# dat[TABLE_NAME][COLUMN_NAME] = []   ## list of "double dictionary" 
# dat[TABLE_NAME][COLUMN_NAME][0], dat[TABLE_NAME][COLUMN_NAME][1], dat[TABLE_NAME][COLUMN_NAME][2],...
#
# * All the data has "string" type.
#
# >examples <
#    _entry.id              -> dat['entry']['id']              = [] 
#    _database_PDB_rev.num  -> dat['database_PDB_rev']['num']  = []
#    _database_PDB_rev.date -> dat['database_PDB_rev']['date'] = []

  contents = read_file_contents_as_one_string(ifname,focus_category=focus_category)
  #print "#len_contents %d"%(len(contents))
  #print "CONTENTS '%s'"%(contents)
  fields = split_line_into_fields(contents)
  #for x in fields:
  #  print "[%s]"%(x) 

  translate_mmCIF_fields_into_dic(dat,fields)
  #for x in (dat.keys()):
  #  print "key:'%s' value:'%s'"%(x,dat[x].keys())
  #sys.exit(1)
  contents = ''
  fields = ''




def split_line_into_fields(line):
  #print "#split_line_into_fields('%s'):"%(line)
#
# * Quatated-words              'abc def ghi' 
# * Double-quatated-words       "abc def ghi"
# * line-head-semicolonr-words \fabc def ghi\f are treated as one single field.
#
# >> example <<
# "_database_2.database_id     PDB"  --> ["_database_2.database_id","PDB"]
# "_database_PDB_rev_record.details   ?" --> ["_database_PDB_rev_record.details","?"]
# "_exptl.method            'X-RAY DIFFRACTION'" -- >["_exptl.method",  "X-RAY DIFFRACTION"]
# "_struct_keywords.pdbx_keywords   'COMPLEX (TRANSFERASE/CYCLIN)'" --> ["_struct_keywords.pdbx_keywords",   "'COMPLEX (TRANSFERASE/CYCLIN)'"]
  state = 'A'   # 'A'lphnum, ' ':Spacer, ''' quatation_in, '"':double-quatation_in, ';':semicolon-in 
  Nsta = 0 
  fields = []
  xnext = ''
  for i in range(len(line)):
    x = line[i]

    xnext = ''
    if (i<(len(line)-1)):
      xnext = line[i+1]

    #print "%d %s %s"%(i,x,state) 
    if (state=='A'):
      if (x==' '):
        state = ' '
        fields.append(line[Nsta:i]) 
    elif (state == ' '):
      if (x != ' ') and (x != '\t') and (x != '"'):
        state = 'A'
        Nsta = i 
      if (x == "'"):
        state = "'" 
        Nsta = i+1 
      if (x == '"'):
        state =  '"' 
        Nsta = i+1 
      if (x == '\f'):
        state =  ';' 
        Nsta = i+1 

    elif (state == "'"):
      if (x == "'") and ((xnext == ' ') or (xnext == '')):
        fields.append(line[Nsta:i]) 
        state = ' '

    elif (state == '"'):
      if (x == '"') and ((xnext == ' ') or (xnext == '')):
        fields.append(line[Nsta:i]) 
        state = ' '

    elif (state == ';'):
      if (x == '\f'):
        fields.append(line[Nsta:i]) 
        state = ' '

  if (state=='A'):
    fields.append(line[Nsta:]) 

  return(fields)


def read_file_contents_as_one_string_original(ifname):

  if not os.access(ifname,os.R_OK):
    if  os.access(ifname+'.gz',os.R_OK):
      f = os.popen('gunzip -c ' + ifname + '.gz')
    else:
      print "#WARNING:Can't open filename '%s' and '%s'" %(ifname,ifname+'.gz')
      print "#ERROR(read_file_contents_as_one_string):Can't open file '%s'"%(ifname)
      return("")
  else:
    if (ifname.endswith('.gz')):
      f = os.popen('gunzip -c ' + ifname)
      #print "#'gunzip -c %s'"%(ifname)
    else:
      f = open(ifname)

  contents = ''
  for line in f:
    if (line.startswith('#')==0):
      line = line.replace('\n',' ')
      line = line.replace('\t',' ')
      ## linehead-semicolon ';' is replace to '\f'
      if (line.startswith(';')):
        line = '\f' + line[1:]
        pass

      contents += line
  f.close()
  #print "#len_contents %d"%(len(contents))
  return(contents)



def read_file_contents_as_one_string(ifname,focus_category=[]):
  #print "#read_file_contents_as_one_string('%s',focus_category=(%s) len %d)"%(ifname,focus_category,len(focus_category))
  if not os.access(ifname,os.R_OK):
    if  os.access(ifname+'.gz',os.R_OK):
      f = os.popen('gunzip -c ' + ifname + '.gz')
    else:
      print "#WARNING:Can't open filename '%s' and '%s'" %(ifname,ifname+'.gz')
      print "#ERROR(read_file_contents_as_one_string):Can't open file '%s'"%(ifname)
      return("")
  else:
    if (ifname.endswith('.gz')):
      f = os.popen('gunzip -c ' + ifname)
      #print "#'gunzip -c %s'"%(ifname)
    else:
      f = open(ifname)

  contents = ''
  table_str = ''
  category     = 'initial'
  category_pre = ''
  for line in f:
    if (line.startswith('#')==0):
      line = line.replace('\n',' ')
      line = line.replace('\t',' ')
      ## linehead-semicolon ';' is replace to '\f'
      if (line.startswith(';')):
        line = '\f' + line[1:]
        pass

      if (line.startswith('loop_')):
        if (table_str != '') and (category != ''):
          if (len(focus_category)==0) or (category in focus_category):
            #print "CATEGORY1 '%s' TABLE_STR '%s'"%(category,table_str)
            contents += table_str
          table_str = '' 
          category = ''
      if (line.startswith('_')):
        fields = line.split('.')
        category = fields[0][1:]
        if (category != category_pre) and (category != '') and (category_pre != ''):
          if (len(focus_category)==0) or (category_pre in focus_category):
            #print "CATEGORY2 '%s' TABLE_STR '%s'"%(category_pre,table_str)
            contents += table_str
          table_str = '' 

      #print "'%s' '%s' '%s'"%(category,category_pre,line)
      category_pre = category  
      table_str += line  

  if (len(focus_category)==0) or (category in focus_category):
    contents += table_str
  f.close()
  #print "#len_contents %d"%(len(contents))
  return(contents)




def translate_mmCIF_fields_into_dic(dat,fields):
#>> STRUCTURE OF VARIABLE dat{}
# dat[TABLE_NAME][COLUMN_NAME] = []
# dat[TABLE_NAME][COLUMN_NAME][0], dat[TABLE_NAME][COLUMN_NAME][1], dat[TABLE_NAME][COLUMN_NAME][2],...
# >example <
#    _entry.id              -> dat['entry']['id'] 
#    _database_PDB_rev.num  -> dat['database_PDB_rev']['num']
#    _database_PDB_rev.date -> dat['database_PDB_rev']['date']
#
# >>> example of input fields[] <<<                      [status]
# data_4HHB                                              'SCALAR'
# _entry.id   4HHB                                       'SCALAR'
# _audit_conform.dict_name       mmcif_pdbx.dic          'SCALAR'
# _audit_conform.dict_version    4.007                   'SCALAR'
# _audit_conform.dict_location   http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx.dic
# _database_2.database_id     PDB                        'SCALAR'
# _database_2.database_code   4HHB                       'SCALAR'
# loop_                                                  'LOOP_HEAD'
# _database_PDB_rev.num                                  'LOOP_KEY'
# _database_PDB_rev.date                                 'LOOP_KEY'
# _database_PDB_rev.date_original                        'LOOP_KEY'
# _database_PDB_rev.status                               'LOOP_KEY'
# _database_PDB_rev.replaces                             'LOOP_KEY'
# _database_PDB_rev.mod_type                             'LOOP_KEY'
# 1 1984-07-17 1984-03-07 ? 4HHB 0                       'LOOP_VALUE'
# 2 1989-10-15 ?          ? 4HHB 3                       'LOOP_VALUE'
# 3 2003-04-01 ?          ? 4HHB 1                       'LOOP_VALUE'
# 4 2009-02-24 ?          ? 4HHB 1                       'LOOP_VALUE'
# 5 2011-07-13 ?          ? 4HHB 1                       'LOOP_VALUE'
# loop_                                                  'LOOP_HEAD'
# _database_PDB_rev_record.rev_num                       'LOOP_KEY'
# _database_PDB_rev_record.type                          'LOOP_KEY'
# _database_PDB_rev_record.details                       'LOOP_KEY'
# 2 MTRIX ?                                              'LOOP_VALUE'
# 3 JRNL  ?                                              'LOOP_VALUE'
# 4 VERSN ?                                              'LOOP_VALUE'
# 5 VERSN ?                                              'LOOP_VALUE'
# _pdbx_database_PDB_obs_spr.id               SPRSDE     'SCALAR'
# _pdbx_database_PDB_obs_spr.date             1984-07-17 'SCALAR'
# _pdbx_database_PDB_obs_spr.pdb_id           4HHB       'SCALAR'
# _pdbx_database_PDB_obs_spr.replace_pdb_id   1HHB       'SCALAR'

  #print "#translate_mmCIF_fields_into_dic(dat,fields):"
  state = 'SCALAR'
  loop_keycols = []
  Nfields = len(fields)
  i = 0
  while (i<Nfields):
    #print "[%d] '%s' status '%s'"%(i,fields[i],state)
    if fields[i].startswith('loop_'):
      state = 'LOOP_HEAD'
      loop_keycols = []
      i = i + 1
    elif (state == 'SCALAR') and fields[i].startswith('_'):
      key   = fields[i][1:]
      value = fields[i+1]
      #print "key '%s' value '%s'"%(key,value)
      (keytab,keycol) = key.split('.') 
      if (dat.has_key(keytab)==0):
        dat[keytab] = {} 
      if (dat[keytab].has_key(keycol)==0):
        dat[keytab][keycol] = []
      dat[keytab][keycol].append(value)
      i = i + 2
    elif ((state == 'LOOP_HEAD') or (state == 'LOOP_KEY')) and fields[i].startswith('_'):
      state = 'LOOP_KEY'
      key   = fields[i][1:]
      (keytab,keycol) = key.split('.') 
      loop_keycols.append(keycol)
      if (dat.has_key(keytab)==0):
        dat[keytab] = {} 
      if (dat[keytab].has_key(keycol)==0):
        dat[keytab][keycol] = []
      i = i + 1 
    elif (state == 'LOOP_KEY') and (fields[i].startswith('_')==0):
      state = 'LOOP_VALUE'
      Nloopkey = len(loop_keycols)
      # repeat reading [Nloopkey] value-fields #
      j = i
      proper_Nloopkey_fields = 1
      while (j<Nfields) and (proper_Nloopkey_fields==1):
        proper_Nloopkey_fields = 1

        if ((j+Nloopkey)>Nfields):
          proper_Nloopkey_fields = 0
        elif (fields[j].startswith('_')):
          proper_Nloopkey_fields = 0
        elif (fields[j]=='loop_'):
          proper_Nloopkey_fields = 0

        if (proper_Nloopkey_fields==1):
          for k in range(Nloopkey):
            if ((j+k)<Nfields):
              dat[keytab][loop_keycols[k]].append(fields[j+k])
              #print "keytab '%s' loop_keycol '%s' field '%s' "%(keytab,loop_keycols[k],fields[j+k])
          j = j + Nloopkey
      i = j
      state = 'SCALAR' 
    else:
      i = i + 1 
 


 
def make_entity_from_mmCIFdic(M,ENTITY,ExHOH_type):
  if ('entity' in M) and ('id' in M['entity']):
    for i in range(len(M['entity']['id'])):
      #print i
      type    = M['entity']['type'][i]
      if (ExHOH_type !='T') or (type !='water'):
        entity = {}
        entity['type'] = type
        entity['type'] = entity['type'][0:16]

        print "'%s' '%s' '%s' '%s'"%(M['entity']['id'][i],M['entity']['type'][i],M['entity']['src_method'][i],M['entity']['pdbx_description'][i])
        entity['entity_id'] = int(M['entity']['id'][i])

        pdbx_description = M['entity']['pdbx_description'][i]
        entity['pdbx_description'] = pdbx_description.replace("'","''")
        entity['pdbx_description'] = entity['pdbx_description'][0:512] 

        entity['pdb_id'] = M['entry']['id'][0].lower()

        entity['poly_type'] = ''
        entity['uniprot_id'] = ''
        entity['uniprot_acc'] = ''
        ENTITY.append(entity)

  if ('entity_poly' in M):
    for i in range(len(M['entity_poly']['type'])):
      entity_id = M['entity_poly']['entity_id'][i]
     # print "entity_id '%s' poly_type '%s'"%(M['entity_poly']['entity_id'][i],M['entity_poly']['type'])
      for entity in (ENTITY):
        if (entity['entity_id']==int(entity_id)):
          entity['poly_type'] = M['entity_poly']['type'][i]
          entity['poly_type'] = entity['poly_type'][0:64]


  if ('struct_ref' in M):
      for i in range(len(M['struct_ref']['db_name'])):
        entity_id = M['struct_ref']['entity_id'][i]
        db_name   = M['struct_ref']['db_name'][i]
        if (db_name == 'UNP'):
          for entity in (ENTITY):
            if (entity['entity_id']==int(entity_id)):
              if ('db_code' in M['struct_ref']):
                db_code           = M['struct_ref']['db_code'][i]
              else:
                db_code = ''
              if ('pdbx_db_accession' in M['struct_ref']):
                pdbx_db_accession = M['struct_ref']['pdbx_db_accession'][i]
              else:
                pdbx_db_accession = ''
              if (entity['uniprot_id'] != ''):
                entity['uniprot_id'] += ':'
              entity['uniprot_id'] += db_code 
              if (entity['uniprot_acc'] != ''):
                entity['uniprot_acc'] += ':'
              entity['uniprot_acc'] +=  pdbx_db_accession 
      
  pass 




def make_unitmol_from_mmCIFdic(M,UNITMOL,ENTITY):
  AfrmTRI = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','MSE':'m','DA':'a','DT':'t','DG':'g','DC':'c','A':'a','T':'t','G':'g','C':'c','U':'u'} 

  Natom      = {}
  Nheavyatom = {}
  entity_id_from_asym_id = {}
  auth_asym_id_from_asym_id = {}
  comp_id = {}

  if ('atom_site' not in M):
    return(0)
  ## (1) checking "atom_site" table ##
  Matm = M['atom_site']
  for i in range(len(Matm['id'])):
    one_model_out = 0
    if ('pdbx_PDB_model_num' in Matm):
      if (Matm['pdbx_PDB_model_num'][i]=='1'):
        one_model_out = 1
    else:
      one_model_out = 1

    if (one_model_out==1):
      asym_id   = Matm['label_asym_id'][i]
      if (asym_id not in Natom):
        Natom[asym_id]   = 0
        Nheavyatom[asym_id]   = 0
        comp_id[asym_id]  = {} 
      Natom[asym_id] += 1
      if (Matm['type_symbol'][i] != 'H'):
        Nheavyatom[asym_id] += 1

      entity_id_from_asym_id[asym_id] =  Matm['label_entity_id'][i]
      auth_asym_id_from_asym_id[asym_id] = Matm['auth_asym_id'][i]
      seq_id = Matm['label_seq_id'][i]
      if (seq_id == '.'):
        seq_id = Matm['auth_seq_id'][i]
      comp_id[asym_id][seq_id] = Matm['label_comp_id'][i]

  ### (2) checking "entity_poly" table ###
  seq_entity = {}
  if ('entity_poly' in M):
    Menpoly = M['entity_poly']
    for i in range(len(Menpoly['entity_id'])):
      seq_entity[Menpoly['entity_id'][i]] = Menpoly['pdbx_seq_one_letter_code'][i].replace(' ','')

    #for entity_id in (seq_entity.keys()):
    #  print "entity '%s' seq '%s'"%(entity_id,seq_entity[entity_id])

  ## (3) Finaly making entity and append into ENTITY[] ##
  for asym_id in Natom.keys():
    entity_id  = entity_id_from_asym_id[asym_id]
    entity_exist = 0
    entity_cor = {} ## entity corresponding the unitmol
    for e in (ENTITY):
      if (int(e['entity_id'])==int(entity_id)):
        entity_exist = 1
        entity_cor = e

    if (entity_exist==1):
      unitmol = {}
      unitmol['asym_id'] = asym_id
      unitmol['natom']   = Natom[asym_id]
      unitmol['nheavyatom']   = Nheavyatom[asym_id]
      unitmol['entity_id']  = entity_id_from_asym_id[asym_id]
      unitmol['auth_asym_id']  = auth_asym_id_from_asym_id[asym_id]
      seq_id_list = sorted(comp_id[asym_id].keys(),lambda x,y:cmp(int(x),int(y))) 
      unitmol['nresidue']   = len(seq_id_list)
      unitmol['sequence']   = ''

      if (unitmol['entity_id']  in seq_entity):
        for i in range(len(seq_entity[unitmol['entity_id']])):
          if str(i+1) in comp_id[asym_id]: 
            unitmol['sequence'] += AfrmTRI.get(comp_id[asym_id][str(i+1)],'X') 
          else:
            #unitmol['sequence'] += seq_entity[unitmol['entity_id']][i].lower()
            unitmol['sequence'] += 'x'

      unitmol['pdb_id'] = M['entry']['id'][0].lower()
      unitmol['type'] = ''
      unitmol['poly_type'] = ''
      unitmol['pdbx_description'] = ''

      unitmol['type']             = entity_cor['type']
      unitmol['poly_type']        = entity_cor['poly_type']
      unitmol['pdbx_description'] = entity_cor['pdbx_description']
      unitmol['uniprot_id']       = entity_cor['uniprot_id']
      unitmol['uniprot_acc']      = entity_cor['uniprot_acc']

      unitmol['comp_id'] = '' 
      if (unitmol['type']!='polymer') and (unitmol['nresidue']==1):        
        unitmol['comp_id'] = comp_id[asym_id][seq_id_list[0]] 

      #print unitmol
      UNITMOL.append(unitmol) 
   
  return(1)



def make_asmblmol_from_mmCIFdic(M,ASMBLMOL,UNITMOL):
  pdb_id = M['entry']['id'][0].lower()
  print "#make_asmblmol_from_mmCIFdic(pdb_id '%s')"%(pdb_id)

  if ('pdbx_struct_assembly_gen' in M) and ('assembly_id' in M['pdbx_struct_assembly_gen']):
    for i in range(len(M['pdbx_struct_assembly_gen']['assembly_id'])): 
      assembly_id = M['pdbx_struct_assembly_gen']['assembly_id'][i]
      oper_expression_list = get_oper_expression_list(M['pdbx_struct_assembly_gen']['oper_expression'][i])
      asym_id_list = []
      if (M['pdbx_struct_assembly_gen']['asym_id_list'][i].find(',')):
        asym_id_list = M['pdbx_struct_assembly_gen']['asym_id_list'][i].split(',')
      else:
        asym_id_list.append(M['pdbx_struct_assembly_gen']['asym_id_list'][i])
      for oper in (oper_expression_list):
        for asym_id in (asym_id_list):
          asmblmol_exist = 0
          for a in (ASMBLMOL):
            if (a['asym_id']==asym_id) and (a['oper_expression']==oper):
              asmblmol_exist = 1
              a['assembly_ids'].append(assembly_id) 
          if (asmblmol_exist==0):
            ## checking asym_id in UNITMOL[] to delete 'HOH' waters ##
            for u in (UNITMOL):
              if (u['asym_id']==asym_id):
                a = {}
                a['pdb_id']  = pdb_id
                a['asym_id'] = asym_id
                a['oper_expression'] = oper 
                a['assembly_ids'] = []
                a['assembly_ids'].append(assembly_id)
                a['entity_id'] = u['entity_id'] 
                ASMBLMOL.append(a)








def make_assembly_from_mmCIFdic(M,ASSEMBLY,UNITMOL):
  pdb_id = M['entry']['id'][0].lower()

  if ('pdbx_struct_assembly' in M):
    for i in range(len(M['pdbx_struct_assembly']['id'])): 
      assembly_id = M['pdbx_struct_assembly']['id'][i]
      a = {}
      a['details'] = 'NULL'
      a['method_details'] = 'NULL'
      a['oligomeric_count'] = 'NULL'
      a['absa'] = 'NULL'
      a['ssa']  = 'NULL'
      a['more'] = 'NULL'

      a['pdb_id'] = pdb_id
      a['assembly_id'] = assembly_id

      if ('details' in M['pdbx_struct_assembly']):
        a['details']        = M['pdbx_struct_assembly']['details'][i]
      if ('method_details' in M['pdbx_struct_assembly']):
        a['method_details'] = M['pdbx_struct_assembly']['method_details'][i]
      if ('oligomeric_count' in M['pdbx_struct_assembly']):
        a['oligomeric_count'] = M['pdbx_struct_assembly']['oligomeric_count'][i]
        if (a['oligomeric_count'].isdigit() !=True):
          a['oligomeric_count'] = 'NULL'

      ASSEMBLY.append(a)

  if ('pdbx_struct_assembly_gen' in M) and ('assembly_id' in M['pdbx_struct_assembly_gen']):
    for i in range(len(M['pdbx_struct_assembly_gen']['assembly_id'])): 
      assembly_id = M['pdbx_struct_assembly_gen']['assembly_id'][i]
      for a in (ASSEMBLY):
        if (a['assembly_id']==assembly_id):
          oper_expression_list = get_oper_expression_list(M['pdbx_struct_assembly_gen']['oper_expression'][i])
          asym_id_list = []
          if (M['pdbx_struct_assembly_gen']['asym_id_list'][i].find(',')):
            asym_id_list = M['pdbx_struct_assembly_gen']['asym_id_list'][i].split(',')
          if ('asym_ids' not in a):
            a['asym_ids'] = []
          if ('oper_expressions' not in a):
            a['oper_expressions'] = []

          for asym_id in (asym_id_list):
            ## checking asym_id in UNITMOL[] to delete 'HOH' waters ##
            unitmol_exist = 0
            for u in (UNITMOL):
              if (u['asym_id']==asym_id):
                unitmol_exist = 1
            if (unitmol_exist == 1):
              for oper in (oper_expression_list):
                a['asym_ids'].append(asym_id)
                a['oper_expressions'].append(oper)

  if ('pdbx_struct_assembly_prop' in M):
    for i in range(len(M['pdbx_struct_assembly_prop']['biol_id'])):
      assembly_id = M['pdbx_struct_assembly_prop']['biol_id'][i]
      for a in (ASSEMBLY):
        if (a['assembly_id']==assembly_id):
          if (M['pdbx_struct_assembly_prop']['type'][i]=='ABSA (A^2)'):
            a['absa']     = M['pdbx_struct_assembly_prop']['value'][i]
            a['absa'] = check_float_str(a['absa'])
          if (M['pdbx_struct_assembly_prop']['type'][i]=='SSA (A^2)'):
            a['ssa']     = M['pdbx_struct_assembly_prop']['value'][i]
            a['ssa'] = check_float_str(a['ssa'])
          if (M['pdbx_struct_assembly_prop']['type'][i]=='MORE'):
            a['more']     = M['pdbx_struct_assembly_prop']['value'][i]
            a['more'] = check_float_str(a['more'])

  ### assign 'nasmblmol' ###
  for a in (ASSEMBLY):
    a['nasmblmol'] = len(a['asym_ids'])
    #print "#assembly_id %s nasmblmol %d"%(a['assembly_id'],len(a['asym_ids']))


def make_pdb_dic_from_mmCIFdic(M,PDB):
  PDB['pdb_id']            = ''
  PDB['title']             = ''
  PDB['exptl_method']      = ''
  PDB['resolution_high']   = 'NULL'
  PDB['rev_date_original'] = ''
  PDB['rev_date_latest']   = ''
  PDB['date_insert']       = ''

  if ('entry' in M) and ('id' in M['entry']):
    PDB['pdb_id']  = M['entry']['id'][0].lower() 

  if ('struct' in M) and ('title' in M['struct']):
    PDB['title']  = M['struct']['title'][0][0:255].replace("'","''") 

  if ('exptl' in M) and ('method' in M['exptl']):
    PDB['exptl_method']    = M['exptl']['method'][0].replace("'","''")

  if ('refine' in M) and ('ls_d_res_high' in M['refine']):
    if (re.match(r'^[0-9]+[\.]*[0-9]*$',M['refine']['ls_d_res_high'][0])):
      PDB['resolution_high'] = "'" + M['refine']['ls_d_res_high'][0] + "'"

  if ('database_PDB_rev' in M):
    PDB['rev_date_original'] =  M['database_PDB_rev']['date_original'][0] 
    PDB['rev_date_latest']   =  M['database_PDB_rev']['date'][len(M['database_PDB_rev']['date'])-1]
  
  now = datetime.now()
  PDB['date_insert'] = now.strftime("%Y-%m-%d") 

  pass



###########################################
##### FUNCTIONS TO WRITE IN PDB FORMAT ####
###########################################



def write_asym_id_molecule_in_PDB(ofname,opt_dic,M,asym_id,one_alt_id='T',residue_select='F'):
#def write_asym_id_molecule_in_PDB(ofname,opt_dic,M,asym_id,entity_id,auth_asym_id,one_alt_id='T'):
  print "#write_mmCIF_asym_id_molecule_in_PDB()-->'%s'"%(ofname)
  ### (1) Setup table between auth_seq_id and label_seq_id ###
  auth_seq_id_from_label_seq_id = {}
  label_seq_id_list = []
  auth_asym_id = '' 
  Natom = len(M['atom_site']['group_PDB'])
  A   = M['atom_site']
  for i in range(Natom):
    if (A['label_asym_id'][i] == asym_id):
      auth_asym_id = A['auth_asym_id'][i]
      if (A['label_seq_id'][i] not in auth_seq_id_from_label_seq_id):
        label_seq_id_list.append(A['label_seq_id'][i])
      auth_seq_id_from_label_seq_id[A['label_seq_id'][i]] = A['auth_seq_id'][i]
  ### (2) write asym_id_molecule in PDB ###
  if (ofname=='-'):
    of = sys.stdout
  else:
    of = open(ofname,"w")
  mol_id = "%s_%s"%(M['entry']['id'][0].lower(),asym_id)
  of.write("HEADER    %-11s                             %-9s   %s\n"%(mol_id,PDB_HEADER_date(M['database_PDB_rev']['date_original'][0]),M['entry']['id'][0].upper()))
#           HEADER    ELECTRON TRANSPORT                      25-FEB-94   1CCH
#           HEADER    COMPLEX (TRANSFERASE/CYCLIN)            14-JUL-96   1FIN
  of.write("REMARK  [pdb_id] %s\n"%(M['entry']['id'][0].lower()))
  of.write("REMARK  [asym_id]      %s\n"%(asym_id))
  of.write("REMARK  [auth_asym_id] %s\n"%(auth_asym_id))
  #of.write("REMARK  [entity_id]    %s\n"%(entity_id))
  of.write("REMARK  [COMMAND] %s\n"%(opt_dic.get('COMMAND','')))
  of.write("REMARK  [DATE]    %s\n"%(opt_dic.get('START_DATE','')))
  of.write("REMARK  [FILENAME] %s\n"%(ofname))
  of.write("REMARK SEQ_ID[label  auth][label auth][label auth][label auth][label auth]")
           #REMARK SEQ_ID    22   705    23   706    24   707    25   708    26   709
  n = 0
  for label_seq_id in (label_seq_id_list):
    if ((n%5)==0):
      of.write("\nREMARK SEQ_ID")
    of.write(" %5s %5s"%(label_seq_id,auth_seq_id_from_label_seq_id[label_seq_id]))
    n += 1
  of.write("\n")
  of.write("REMARK  [1]:auth_asym_id [2]:label_seq_id [3] type_symbol\n")
  of.write("REMARK              [1] [2]                                                [3]\n")
  #         HETATM 1237 FE   HEM A 155      15.271  27.962   0.622  1.00  7.86          FE

  A = M['atom_site']
  for i in range(len(A['label_asym_id'])):
    if (A['label_asym_id'][i] == asym_id):
      if (decide_output_ATOM_HETATM_line(A,i,one_model_out='T',one_alt_id='T',ExHOH='T',residue_select=residue_select)==1):
        write_ATOM_HETATM_line(of,A,i,A['Cartn_x'][i],A['Cartn_y'][i],A['Cartn_z'][i],one_alt_id='T')
  if (ofname!='-'):
    of.close()







def write_asym_id_oper_expresssion_molecule_in_PDB(ofname,opt_dic,M,asym_id,oper_expression,assembly_ids=[],one_alt_id='T',residue_select='F'):
  print "#write_mmCIF_asym_id_oper_expression_molecule_in_PDB(asym_id:%s oper_expression:%s)-->'%s'"%(asym_id,oper_expression,ofname)
  Natom = len(M['atom_site']['group_PDB'])
  print "#Natom %d"%(Natom) 
  Matm   = M['atom_site']
  Masbl  = M['pdbx_struct_assembly']
  Masblg = M['pdbx_struct_assembly_gen']
  Moper  = M['pdbx_struct_oper_list']
  R = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
  T = [0.0,0.0,0.0]

  ## (1) Find number of assembly, assemly_gen, oper_list, and setup R matrix and T vector ###

  setup_Rmatrix_and_Tvector_from_oper_expression(R,T, M, oper_expression)

  ### (2) Setup table between auth_seq_id and label_seq_id ###
  auth_seq_id_from_label_seq_id = {}
  label_seq_id_list = []
  auth_asym_id = '' 
  for i in range(Natom):
    if (Matm['label_asym_id'][i] == asym_id):
      auth_asym_id = Matm['auth_asym_id'][i]
      if (Matm['label_seq_id'][i] not in auth_seq_id_from_label_seq_id):
        label_seq_id_list.append(Matm['label_seq_id'][i])
      auth_seq_id_from_label_seq_id[Matm['label_seq_id'][i]] = Matm['auth_seq_id'][i]

  ### (4) write the transformed XYZ in PDB format ###
  of = open(ofname,"w")
  asmblmol_id = "%s_%s_%s"%(M['entry']['id'][0].lower(),asym_id,oper_expression)
  of.write("HEADER    %-11s                             %-9s   %s\n"%(asmblmol_id,PDB_HEADER_date(M['database_PDB_rev']['date_original'][0]),M['entry']['id'][0].upper()))
  #         HEADER    KINASE                                  10-JUN-03   1OI9
  #         HEADER    COMPLEX (TRANSFERASE/CYCLIN)            14-JUL-96   1FIN
  #         HEADER    2cch_1_G_1                              06-JAN-16   2CCH
  #         REMARK 300 THIS ENTRY. THE REMARK MAY ALSO PROVIDE INFORMATION ON
  #         REMARK 300 BURIED SURFACE AREA.
  of.write("REMARK     [pdb_id] %s\n"%(M['entry']['id'][0].lower()))
  of.write("REMARK     [asym_id]         %s\n"%(asym_id))
  of.write("REMARK     [auth_asym_id]    %s\n"%(auth_asym_id)) 
  of.write("REMARK     [oper_expression] %s\n"%(oper_expression))
  if (len(assembly_ids)>0):
    of.write("REMARK     [assembly_id]")
    for assembly_id in (assembly_ids):
      of.write(" %s"%(assembly_id))
    of.write("\n")

  of.write("REMARK     [COMMAND] %s\n"%(opt_dic['COMMAND']))
  of.write("REMARK     [DATE]    %s\n"%(opt_dic['START_DATE']))
  of.write("REMARK     [FILENAME] %s\n"%(ofname))
  of.write("REMARK     [R0] %+6.3f %+6.3f %+6.3f [T0] %.3f\n"%(R[0][0],R[0][1],R[0][2],T[0]))
  of.write("REMARK     [R1] %+6.3f %+6.3f %+6.3f [T1] %.3f\n"%(R[1][0],R[1][1],R[1][2],T[1]))
  of.write("REMARK     [R2] %+6.3f %+6.3f %+6.3f [T2] %.3f\n"%(R[2][0],R[2][1],R[2][2],T[2]))

  of.write("REMARK SEQ_ID[label  auth][label auth][label auth][label auth][label auth]")
           #REMARK SEQ_ID    22   705    23   706    24   707    25   708    26   709
  n = 0
  for label_seq_id in (label_seq_id_list):
    if ((n%5)==0):
      of.write("\nREMARK SEQ_ID")
    of.write(" %5s %5s"%(label_seq_id,auth_seq_id_from_label_seq_id[label_seq_id]))
    n += 1
  of.write("\n")
  of.write("REMARK     [1]:auth_asym_id [2]:label_seq_id [3] type_symbol\n")
  of.write("REMARK              [1] [2]                                                [3]\n")
#           HETATM 1237 FE   HEM A 155      15.271  27.962   0.622  1.00  7.86          FE

  for i in range(Natom):
    if (Matm['label_asym_id'][i] == asym_id):
      if (decide_output_ATOM_HETATM_line(Matm,i,one_model_out='T',one_alt_id='T',ExHOH='T',residue_select=residue_select)==1):
        x = float(Matm['Cartn_x'][i])
        y = float(Matm['Cartn_y'][i])
        z = float(Matm['Cartn_z'][i])
        nx = R[0][0]*x + R[0][1] * y + R[0][2] * z + T[0]
        ny = R[1][0]*x + R[1][1] * y + R[1][2] * z + T[1]
        nz = R[2][0]*x + R[2][1] * y + R[2][2] * z + T[2]
        write_ATOM_HETATM_line(of,Matm,i,"%8.3f"%(nx), "%8.3f"%(ny), "%8.3f"%(nz),tFactor='',one_alt_id='T')

  of.write("TER\n")
  of.close()


def write_asymmetric_unit_in_PDB(ofname,M,one_alt_id='T',Rfit=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],Tfit=[0.0,0.0,0.0],tFactor='',ExHOH='T',residue_select='F'):
  print "#write_asymmetric_unit_in_PDB()-->'%s'"%(ofname)
  if (ofname=='-'):
    of = sys.stdout
  else:
    print "#write_asymmetric_unit_in_PDB()-->'%s'"%(ofname)
    of = open(ofname,"w")

  Natom = 0;
  if ('atom_site' in M):
    Natom = len(M['atom_site']['group_PDB'])
  else:
    print "#Natom %d"%(Natom)
    return(0)

  of.write("HEADER    %-11s                             %-9s   %s\n"%(M['entry']['id'][0].lower(),PDB_HEADER_date(M['database_PDB_rev']['date_original'][0]),M['entry']['id'][0].upper()))
  of.write("REMARK  [1]:auth_asym_id [2]:label_seq_id [3] type_symbol\n")
  of.write("REMARK              [1] [2]                                                [3]\n")

  A = M['atom_site']
  for i in range(Natom):
    if (decide_output_ATOM_HETATM_line(A,i,one_model_out='T',one_alt_id='T',ExHOH='T',residue_select=residue_select)==1):
      x = float(A['Cartn_x'][i])
      y = float(A['Cartn_y'][i])
      z = float(A['Cartn_z'][i])
      nx = Rfit[0][0]*x + Rfit[0][1] * y + Rfit[0][2] * z + Tfit[0]
      ny = Rfit[1][0]*x + Rfit[1][1] * y + Rfit[1][2] * z + Tfit[1]
      nz = Rfit[2][0]*x + Rfit[2][1] * y + Rfit[2][2] * z + Tfit[2]
      write_ATOM_HETATM_line(of,A,i,"%8.3f"%(nx), "%8.3f"%(ny), "%8.3f"%(nz),tFactor=tFactor,one_alt_id='T')

  if (ofname!='-'):
    of.close()
  return(1)






def write_assembly_in_PDB(ofname,M,assembly_id,Rfit=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],Tfit=[0.0,0.0,0.0],tFactor='',ExHOH='T',one_alt_id='T',residue_select='F'):
  if (ofname=='-'):
    of = sys.stdout
  else:
    print "#write_assembly_in_PDB(assembly_id:%s residue_select:'%s')-->'%s'"%(assembly_id,residue_select,ofname)
    of = open(ofname,"w")

  ## (1) make asym_id_list[] and oper_expression_list[]  
  asym_id_list = []
  oper_expression_list = []

  if ('pdbx_struct_assembly_gen' in M) and ('assembly_id' in M['pdbx_struct_assembly_gen']):
    for i in range(len(M['pdbx_struct_assembly_gen']['assembly_id'])):
      if (M['pdbx_struct_assembly_gen']['assembly_id'][i]==assembly_id):
        oper_expressions = []
        asym_ids         = []
        oper_expressions = get_oper_expression_list(M['pdbx_struct_assembly_gen']['oper_expression'][i])
        if (M['pdbx_struct_assembly_gen']['asym_id_list'][i].find(',')):
          asym_ids = M['pdbx_struct_assembly_gen']['asym_id_list'][i].split(',')

        for asym_id in (asym_ids):
          for oper in (oper_expressions):
            asym_id_list.append(asym_id)
            oper_expression_list.append(oper)
           
  ### (2) write REMARK in assembly information ###
  Masbl  = M['pdbx_struct_assembly']

  nasbl  = -1
  for i in range(len(Masbl['id'])):
    if (Masbl['id'][i] == assembly_id):
      nasbl = i

  of.write("REMARK   pdbx_struct_assembly.id                  %s\n"%(Masbl['id'][nasbl]))
  if ('details' in Masbl):
    of.write("REMARK   pdbx_struct_assembly.details             %s\n"%(Masbl['details'][nasbl]))
  if ('oligomeric_details' in Masbl):
    of.write("REMARK   pdbx_struct_assembly.oligomeric_details  %s\n"%(Masbl['oligomeric_details'][nasbl]))
  if ('oligomeric_count' in Masbl):
    of.write("REMARK   pdbx_struct_assembly.oligomeric_count    %s\n"%(Masbl['oligomeric_count'][nasbl]))
   
  model_num = 0
  for i in range(len(asym_id_list)):
    model_num += 1
    of.write("REMARK   MODEL %4d asym_id %2s oper_expression %2s\n"%(model_num,asym_id_list[i],oper_expression_list[i]))


  ### write transformed coordinates in PDB ###
  Roper = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
  Toper = [0.0,0.0,0.0]
  A = M['atom_site']
  model_num = 0
  
  for i in range(len(asym_id_list)):
    asym_id = asym_id_list[i] 
    oper    = oper_expression_list[i] 
    model_num += 1
    of.write("MODEL %8d\n"%(model_num)) 
    
    setup_Rmatrix_and_Tvector_from_oper_expression(Roper,Toper, M, oper)
    (R,T) = combine_two_Rmat_and_Tvec(Rfit,Tfit,Roper,Toper)

    for j in range(len(A['id'])):
      if (A['label_asym_id'][j] == asym_id):
        if (decide_output_ATOM_HETATM_line(A,j,one_model_out='T',one_alt_id='T',ExHOH='T',residue_select=residue_select)==1):
          x = float(A['Cartn_x'][j])
          y = float(A['Cartn_y'][j])
          z = float(A['Cartn_z'][j])
          nx = R[0][0]*x + R[0][1] * y + R[0][2] * z + T[0]
          ny = R[1][0]*x + R[1][1] * y + R[1][2] * z + T[1]
          nz = R[2][0]*x + R[2][1] * y + R[2][2] * z + T[2]
          write_ATOM_HETATM_line(of,A,j,"%8.3f"%(nx), "%8.3f"%(ny), "%8.3f"%(nz),tFactor=tFactor,one_alt_id='T')
    of.write("TER\n")
    of.write("ENDMDL\n")

  if (ofname != '-'):
    of.close()


def PDB_HEADER_date(datestr):
#'2006-01-16' --> '16-JAN-06'
#'1996-07-14' --> '14-JUL-96'   

  month = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
  field = datestr.split('-')
  year  = field[0][2:]
  m = int(field[1])-1
  if (m>=0) and (m<12):
    month = month[m]
  else:
    month = 'XXX'
  day   =  field[2]
  return(day + '-' + month + '-' + year)



def decide_output_ATOM_HETATM_line(A,i,one_model_out='T',one_alt_id='T',ExHOH='T',residue_select='F'):
  ## A : M['atom_site']

  out= 1

  if (one_model_out == 'T'):
    if ('pdbx_PDB_model_num' in A):
      if (A['pdbx_PDB_model_num'][i]!='1'):
        out = 0
  
  my_label_alt_id = A['label_alt_id'][i]
  if (one_alt_id == 'T') and ((A['label_alt_id'][i] == '.') or (A['label_alt_id'][i] == 'A') or (A['label_alt_id'][i] == '1')):
    my_label_alt_id = ' '

  if (one_alt_id == 'T') and (my_label_alt_id != ' '):
    out = 0

  if (ExHOH=='T') and ((A['auth_comp_id'][i]=='HOH') or (A['auth_comp_id'][i]=='DOD')):
    out = 0
 
  if (residue_select=='T'):
#  >> residue select accept condition <<
#    (atom_id='CA' and type='C') or (atom_id='P' and type='P')
#  >> residue select reject condition <<
#    not [(atom_id='CA' and type='C') or (atom_id='P' and type='P')]
# =  not [(atom_id='CA' and type='C')] and not[(atom_id='P' and type='P')]
# =  (atom_id!='CA' or type!='C')] and (atom_id!='P' or type!='P')
    if ((A['label_atom_id'][i]!='CA') or (A['type_symbol'][i]!='C')) and ((A['label_atom_id'][i]!='P') or (A['type_symbol'][i]!='P')):
      out = 0

  #print "#out %d i %d '%s' '%s' '%s'"%(out,i,A['label_atom_id'][i],A['type_symbol'][i],A['auth_comp_id'][i]) 

  return(out)



def write_ATOM_HETATM_line(of,A,i,Cartn_x,Cartn_y,Cartn_z,tFactor='',one_alt_id='F'):
  ## A : M['atom_site']

  my_label_seq_id = A['label_seq_id'][i]
  if (my_label_seq_id=='.'):
    my_label_seq_id = '1'

  my_label_alt_id = A['label_alt_id'][i]
  if (one_alt_id == 'T') and ((A['label_alt_id'][i] == '.') or (A['label_alt_id'][i] == 'A') or (A['label_alt_id'][i] == '1')):
    my_label_alt_id = ' '

  my_auth_atom_id = mk_pdb_atom_id_from_mmCIF_auth_atom_id(A['auth_atom_id'][i],A['type_symbol'][i])
 
  of.write("%-6s%5d"%(A['group_PDB'][i],int(A['id'][i])%100000))
  of.write(" %-4s%1s%3s %s%4s "%(my_auth_atom_id,my_label_alt_id,A['auth_comp_id'][i],A['auth_asym_id'][i][0],my_label_seq_id))
  of.write("   %8.3f%8.3f%8.3f%6.2f"%(float(Cartn_x),float(Cartn_y),float(Cartn_z),float(A['occupancy'][i])))
  if (tFactor==''):
    of.write("%6.2f"%(float(A['B_iso_or_equiv'][i])))
  else:
    of.write("%6.2f"%(float(tFactor)))
  of.write("          %2s"%(A['type_symbol'][i]))
  of.write("\n")



def mk_pdb_atom_id_from_mmCIF_auth_atom_id(auth_atom_id,type_symbol):
#  [auth_atom_id]
#
#  >> RULE <<
#   if (length(auth_atom_id)<4){
#     if (length(type_symbol)==2) auth_atom_id starts at 1st position, otherwise starts at 2nd position.
#   }
#
#>> 3icb.pdb <<
#          1         2         3         5         5         6         7
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#ATOM    586  CA  SER A  74       4.597   6.256   1.954  1.00 37.65           C
#ATOM    592  CA  GLN A  75       2.756   4.149   4.526  1.00 45.45           C
#HETATM  602 CA    CA A  76      22.527  19.883   5.062  1.00  3.59          CA
#HETATM  603 CA    CA A  77      25.690   9.260   6.305  1.00  4.99          CA
  L            = len(auth_atom_id)
  Ltype_symbol = len(type_symbol)

  if (L==4):
    return(auth_atom_id)
  elif (L==3):
   if (Ltype_symbol==2):
     return("%s "%(auth_atom_id))
   else:
     return(" %s"%(auth_atom_id))
  elif (L==2):
   if (Ltype_symbol==2):
     return("%s  "%(auth_atom_id))
   else:
     return(" %s "%(auth_atom_id))
  elif (L==1):
   if (Ltype_symbol==2):
     return("%s   "%(auth_atom_id))
   else:
     return(" %s  "%(auth_atom_id))
    
  return('')


def check_float_str(x):
  try:
    retstr = '%f'%(float(x))
  except ValueError:
    retstr = 'NULL'
  return(retstr)



def get_oper_expression_list_orig(str0):
# >> example <<
# '1'                -> ['1']
# 'P'                -> ['P']
#  '1,2,3'           -> ['1','2','3']
# '(1-5)'            -> ['1','2','3','4','5']
# '(1-60)'           -> ['1','2',...,'59','60']
# '(1,2,6,10,23,24)' -> ['1','2','6','10','23','24']
# '(X0)(1-60)'       -> ['X0-1','X0-2',...,'X0-59','X0-60']
#  taken from PDBcode:'1m4x' 
# '(1-60)(61-88)'       -> ['1-61','1-62',...,'60-87','60-88']
  str = str0

  str = str.replace(")(",",")
  str = str.replace("(",",")
  str = str.replace(")",",")
  list0 = str.split(",") 
  list = []
  p = re.compile('^[0-9]+\-[0-9]+$') 
  for x in (list0):
    if (p.match(x)!='') and (p.match(x)):
      #print "match '%s' '%s'"%(x,p.match(x))
      (sta,end) = x.split('-')
      for i in range(int(sta),int(end)+1):
        s = "%d"%(i) 
        list.append(s)
    elif (x != ''):
      list.append(x)
  return(list)



def get_oper_expression_list(string):
# >> example <<
# '1'                -> ['1']
# 'P'                -> ['P']
#  '1,2,3'           -> ['1','2','3']
# '(1-5)'            -> ['1','2','3','4','5']
# '(1-60)'           -> ['1','2',...,'59','60']
#  taken from PDBcode:'2buk'
# '(1,2,6,10,23,24)' -> ['1','2','6','10','23','24']
# '(X0)(1-60)'       -> ['X0-1','X0-2',...,'X0-59','X0-60']
#  taken from PDBcode:'1m4x'
# '(1-60)(61-88)'       -> ['1-61','1-62',...,'60-87','60-88']
## >> caution <<
## '[t2]-[t1]' means that firstly, the transformation 't1' is done, 
##  next, the transformation 't2' is done.
 
  field = string.split(')(')
  oper_list = [[] for i in range(len(field))]

  for i in range(len(field)):
    #print "[%d] '%s'"%(i,field[i])
    field[i] = field[i].lstrip("(")
    field[i] = field[i].rstrip(")")
    item_list  = field[i].split(",")
    oper_list[i] = []
    p = re.compile('^[0-9]+\-[0-9]+$')
    for x in (item_list):
      if (p.match(x)!='') and (p.match(x)):
        (sta,end) = x.split('-')
        for k in range(int(sta),int(end)+1):
          s = "%d"%(k)
          oper_list[i].append(s)
      elif (x != ''):
        oper_list[i].append(x)

  if (len(field)==1):
    return(oper_list[0])
  elif (len(field)==2):
    all_oper_list = []
    for x in (oper_list[0]):
      for y in (oper_list[1]):
        all_oper_list.append("%s-%s"%(x,y))
    return(all_oper_list)
  elif (len(field)==3):
    all_oper_list = []
    for x in (oper_list[0]):
      for y in (oper_list[1]):
        for z in (oper_list[2]):
          all_oper_list.append("%s-%s-%s"%(x,y,z))
    return(all_oper_list)
  else:
    printf("#ERROR:too many ')(' in oper_expression_string '%s'."%(string))




def setup_Rmatrix_and_Tvector_from_oper_expression(R,T, M, oper_expression):
  field = oper_expression.split('-')

### (1) Single operation ###
  if (len(field)==1):
    setup_Rmatrix_and_Tvector_from_single_oper_expression(R,T, M, oper_expression)

### (2) Double operation ###
  elif (len(field)==2):
# '[t2]-[t1]' = (R2,T2)-(R2,T1)
#
# x1  = R1*x0 + T1
# x2  = R2*x1 + T2
#     = R2*(R1*x0 + T1) + T2
#     = R2*R1*x0 + R2*T1 + T2
#     = R*    x0 + T
# Threfore,
#    R = R2*R1,   T = R2*T1 + T2

    R1 = [[0.0,0.0,0.0] ,[0.0,0.0,0.0] ,[0.0,0.0,0.0] ]
    T1 = [0.0,0.0,0.0]
    R2 = [[0.0,0.0,0.0] ,[0.0,0.0,0.0] ,[0.0,0.0,0.0] ]
    T2 = [0.0,0.0,0.0]
    t2_string = field[0]
    t1_string = field[1]
    setup_Rmatrix_and_Tvector_from_single_oper_expression(R1,T1, M,t1_string)
    setup_Rmatrix_and_Tvector_from_single_oper_expression(R2,T2, M,t2_string)
    (R,T) = combine_two_Rmat_and_Tvec(R2,T2,R1,T1)


 

def setup_Rmatrix_and_Tvector_from_single_oper_expression(R,T, M, oper_expression):
  Moper  = M['pdbx_struct_oper_list']
  nop = -1
  for i in (range(len(Moper['id']))):
    if (Moper['id'][i]==oper_expression):
      nop = i
  if (nop>=0):
    for i in range(3):
      T[i] = float(Moper["vector[%d]"%(i+1)][nop])
      for j in range(3):
        R[i][j] = float(Moper["matrix[%d][%d]"%(i+1,j+1)][nop])


def combine_two_Rmat_and_Tvec(R2,T2,R1,T1):
# '[t2]-[t1]' = (R2,T2)-(R2,T1)
#
# x1  = R1*x0 + T1
# x2  = R2*x1 + T2
#     = R2*(R1*x0 + T1) + T2
#     = R2*R1*x0 + R2*T1 + T2
#     = R*    x0 + T
# Threfore,
#    R = R2*R1,   T = R2*T1 + T2
  R = [[1.0,0.0,0.0] ,[0.0,1.0,0.0] ,[0.0,0.0,1.0] ]
  T = [0.0,0.0,0.0]

#   R = R2*R1
  for i in range(3):
    for j in range(3):
      R[i][j] = 0.0
      for k in range(3):
        R[i][j] += R2[i][k]*R1[k][j]
#   T = R2*T1 + T2
  for i in range(3):
    T[i] = T2[i] 
    for k in range(3):
      T[i] += R2[i][k]*T1[k]
  return(R,T)




#############
### MAIN ####
#############


def _main():

  OPT = {}
  OPT['O'] = 'assembly'
  OPT['res'] = 'F'
  if (len(sys.argv)<2):
    print "mmCIF.py [input_mmCIF_file]"
    print "  coded by T.Kawabata. LastModDate:%s"%(LastModDate)
    print "<option>"
    print "-O   : 'assembly','unitmol','asmblmol'  [%s]"%(OPT['O'])
    print "-res : residue_select (T or F)  [%s]"%(OPT['res'])
    sys.exit(1)

  ifname = sys.argv[1]
  read_option(sys.argv,OPT)
  M = {}
 
  #read_mmCIF_file(ifname,dat,focus_category=['pdbx_struct_assembly','pdbx_struct_assembly_gen'])
  read_mmCIF_file(ifname,M)
  ENTITY = []
  ASSEMBLY = []
  UNITMOL = []
 
  make_entity_from_mmCIFdic(M,ENTITY,'T')
  make_unitmol_from_mmCIFdic(M,UNITMOL,ENTITY)
  make_assembly_from_mmCIFdic(M,ASSEMBLY,UNITMOL)
  for a in (ASSEMBLY):
    print "%s %s %s"%(a['nasmblmol'],a['oper_expressions'],a['asym_ids'])


  Nasym_id = {}
  comp_id_of_asym_id  = {}
  asym_id_list = []
  a = M['atom_site']
  for i in range(len(M['atom_site']['id'])):
    Nasym_id[a['label_asym_id'][i]] = Nasym_id.get(a['label_asym_id'][i],0) + 1
    if (a['label_asym_id'][i] not in comp_id_of_asym_id.keys()):
      comp_id_of_asym_id[a['label_asym_id'][i]] = a['auth_comp_id'][i]
      asym_id_list.append(a['label_asym_id'][i])
    else:
      comp_id_of_asym_id[a['label_asym_id'][i]] = ''

  ### ['assembly']: output assembly ###

  if (OPT['O'] == 'assembly'):
    if ('pdbx_struct_assembly' in M):
      A = M['pdbx_struct_assembly']
      G = M['pdbx_struct_assembly_gen']
      Nassembly = len(A['id'])
      items = ('id','details','method_details','oligomeric_details','oligomeric_count')
      for i in range(Nassembly):
        for item in (items):
          if (item in A.keys()):
            sys.stdout.write(" %s"%(A[item][i]))
          else: 
            sys.stdout.write("-")
        for j in range(len(G['assembly_id'])):
          if (A['id'][i]==G['assembly_id'][j]):
            sys.stdout.write(" (%s)[%s]"%(G['oper_expression'][j],G['asym_id_list'][j]))
        sys.stdout.write("\n")
   
      tar_assembly_id = raw_input("Input assembly  id >:") 
      opdbfile        = raw_input("Input outputPDBfile>:") 
      if (tar_assembly_id=='asym') or (tar_assembly_id=='asymmetric') or (tar_assembly_id=='asymmetric unit'):
        #write_asymmetric_unit_in_PDB(opdbfile,M)
        write_asymmetric_unit_in_PDB(opdbfile,M,residue_select=OPT['res'])
      else:
        #write_assembly_in_PDB(opdbfile,M,tar_assembly_id)
        write_assembly_in_PDB(opdbfile,M,tar_assembly_id,residue_select=OPT['res'])
    sys.exit(1) 

  ### ['unitmol']: output unitmol ###
  if (OPT['O'] == 'unitmol'):
    for asym_id in (asym_id_list):
      print "ASYM_ID '%-2s' comp_id '%3s' Natom %5d"%(asym_id,comp_id_of_asym_id[asym_id],Nasym_id[asym_id])
   
    tar_asym_id = raw_input("Input asym_id>:") 
    opdbfile    = raw_input("Input outputPDBfile>:") 
    write_asym_id_molecule_in_PDB(opdbfile,OPT,M,tar_asym_id,residue_select=OPT['res'])

  ### ['asmblmol']: output asmblmol ###
  if (OPT['O'] == 'asmblmol'):
    a = M['atom_site']
    oper_expression_list = []
    oper_dic = {}

    for i in range(len(M['pdbx_struct_assembly_gen']['assembly_id'])):
      oper_str   = M['pdbx_struct_assembly_gen']['oper_expression'][i]
      oper_list = get_oper_expression_list(oper_str)
      for oper in (oper_list):
        oper_dic[oper] = 1
    oper_expression_list = sorted(oper_dic.keys(), lambda x,y:cmp(x,y))

    for asym_id in (asym_id_list):
      print "ASYM_ID '%-2s' comp_id '%3s' Natom %5d"%(asym_id,comp_id_of_asym_id[asym_id],Nasym_id[asym_id])
    for id in (oper_expression_list):
      print "OPER    '%-2s'"%(id)


    tar_asym_id         = raw_input("Input asym_id>:") 
    tar_oper_expression = raw_input("Input oper_expression>:") 
    opdbfile    = raw_input("Input outputPDBfile>:") 
    write_asym_id_oper_expresssion_molecule_in_PDB(opdbfile,OPT,M,tar_asym_id,tar_oper_expression,residue_select=OPT['res'])
    pass

if __name__ == '__main__':_main()
                                      


