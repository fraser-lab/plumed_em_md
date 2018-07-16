## 
## <jvxl.py>
##
##  functions for dealing with Jmol JVXL file
##

import os
import sys
import math
from datetime import datetime


LastModDate = "May 20, 2014"

# >> How to interpret jvxlVertexData <<
# Example:
# >> simple.obj<<
# v 0 0 0
# v 1 0 0
# v 0 1 0
# v 0 0 1
# f 1 3 2
# f 1 2 4
# f 4 3 1
# f 4 2 3
#
# ** "*.obj" file consists of two sections : "vertex" section ('v') and "face" section ('f').
#
#
# >> simple.jvxl <<
#<jvxlTriangleData count="4" encoding="jvxltdiff" data="!]]Z^]!Z[_[["></jvxlTriangleData>
#<jvxlVertexData count="4" min="{0.0 0.0 0.0}" max="{1.0 1.0 1.0}" encoding="base90xyz2" data="####|#|####|####|#|####|"></jvxlVertexData>
#
# data of triangle ="!]]Z^]!Z[_[[" should be interpreted as [[1 3 2], [1 2 4], [4 3 1], [4 2 3]]
# data of vertexes ="####|#|####|####|#|####|" should be interpreted as [[0 0 0], [1 0 0], [0 1 0], [0 0 1]].
#
# The first important point is that the order of the vertexes are changed by the "f"ace information.
#
#  --------------------------------------------------------------------------------------------------
# | The order of the vertexes is decided by the appearing order of the vertexes in the face section. |
#  --------------------------------------------------------------------------------------------------
#
# In this case, the original face section is described as:
# f 1 3 2
# f 1 2 4
# f 4 3 1
# f 4 2 3
# By re-ordering the vertex number, this is interpreted as follows:
# f 1 2 3
# f 1 3 4
# f 4 2 1
# f 4 3 2
# It means the numbers the vertexes are replaced as follows:
#   [orig]   [new]
#  v  1   ->  1
#  v  2   ->  3
#  v  3   ->  2
#  v  4   ->  4
#
# Therefore, the obj file shown previously is essentialy interpreted as follows:
# >> simple_iden.obj (essentialy identical to original one)<<
# v 0 0 0
# v 0 1 0
# v 1 0 0
# v 0 0 1
# f 1 2 3
# f 1 3 4
# f 4 2 1
# f 4 3 2
# >>  How to interpret  Vertex Data <<
# data of vertexes ="####|#|####|####|#|####|" should be interpreted as [[0 0 0], [0 1 0], [1 0 0], [0 0 1]].
#
# [V0] "count" is the number of vertexes. In this case, [count]=4. 
#
# [V1] Basically, the string is two-figure('futa-keta') base90 numbers.
#   <base 90>
#   '#':35 => 0
#   '$':36 => 1
#   '%':37 => 2 
#    :
#   '\':92 -> '!' => 57
#   :
#   'z':123 => 87
#   '{':123 => 88
#   '|':124 => 89
#
#   Note that '\' (92=>57) is replaced with '!'.
#   Following this rule, the string "####|#|####|####|#|####|" can be conveted as :
#                   [0,0,0,0,89,0,89,0,0,0,0,89,0,0,0,0,89,0,89,0,0,0,0,89].
#
#   The first [count]*3(=12) number are 'upper' figure, the last [count]*3(=12) numbers are 'lower' figures.
#      The uppers = [0,0,0,0,89,0,89,0,0,0,0,89]
#      The lowers = [0,0,0,0,89,0,89,0,0,0,0,89].
#   The [upper] and [lower] values can be conveted to the float value as follows:
#
#       [value] = ([MAX]-[MIN])(90*[upper]+[lower])/(90*90-1) + [MIN] 
#
#     (0   0,   0  0,  0  0) => (1,0,0)
#     (0   0,  89 89,  0  0) => (0,1,0)
#     (89 89,   0  0,  0  0) => (1,0,0)
#     (0   0,   0  0, 89 89) => (0,0,1)
#
#  [V2] 4 or more continuous repeated symbols can be described in a shorter way usint '~'.
#
#     A x [n] ==> A~[n-1][space]
#
#    ex) 
#    'A'      => 'A~0 '  
#    'AA'     => 'A~1 '
#    'AAA'    => 'A~2 '
#    'AAAA'   => 'A~3 ' same length
#    'AAAAA'  => 'A~4 ' shorter !
#    'AAAAAA' => 'A~5 ' shorter !
#  
#  [V3] Special cares for '<' and '&'
#
#     '<<<<' => '~;3 '
#     '<<<', => '~;2 '
#     '<<',  => '~;1 '
#     '<',   => '~;0 '
#     '&&&&' => '~%3 '
#     '&&&', => '~%2 '
#     '&&',  => '~%1 '
#     '&',   => '~%0 '
#
# >> How to interpret jvxlTriangleData <<
#
#<jvxlTriangleData count="4" encoding="jvxltdiff" data="!]]Z^]!Z[_[["></jvxlTriangleData>
# data of triangle ="!]]Z^]!Z[_[[" should be interpreted as [[1 2 3], [1 3 4], [4 2 1], [4 3 2]]
# 
# [T0] "count" is the number of faces (triangles). In this case, [count]=4. 
#
# [T1] Vertex numbers are described as the "difference" from the previous vertex number.
# 
# For example,
# data of triangle ="!]]Z^]!Z[_[[" is converted as [57, 58, 58, 55, 59, 58, 57, 55, 56, 60, 56, 56]. 
# using the same rules described in [V1],[V2],V[3].
#
# These numbers represents the "difference" from the previous number.
# [57, 58, 58, 55, 59, 58, 57, 55, 56, 60, 56, 56]. 
#   0   1  1   -2   2   1   0  -2  -1   3  -1  -1
# Therefore, it can be interpreted as follows:
# (0 ) 1 (+1) 2 (+1) 3
# (-2) 1 (+2) 3 (+1) 4
# ( 0) 4 (-2) 2 (-1) 1
# (+3) 4 (-1) 3 (-1) 2
#
# However, this strategy only can described -32:(25)-- 0(57)-- +32 (89) 
# When the diffrence value is out of the range [-32:32], 
# it suddenly uses normal "decimal" number, such as '+123' or '-49' following
# the ASCII table.
#[ascii] [80] 
#  35 '#'  0 : __not_use__  
#  36 '$'  1 : __not_use__ 
#  37 '%'  2 : __not_use__
#  38 '&'  3 : __not_use__
#  39 '''  4 : __not_use__
#  40 '('  5 : __not_use__ 
#  41 ')'  6 : __not_use__
#  42 '*'  7 : __not_use__
#  43 '+'  8 : PLUS indicator 
#  44 '.'  9 : __not_use__ 
#  45 '-' 10 : MINUS indicator 
#  46 '.' 11 : __not_use__ 
#  47 '/' 12 : __not_use__ 
#  48 '0' 13 :0 (decimal) 
#  49 '1' 14 :1 (decimal)
#  50 '2' 15 :2 (decimal)
#  51 '3' 16 :3 (decimal)
#  52 '4' 17 :4 (decimal) 
#  53 '5' 18 :5 (decimal) 
#  54 '6' 19 :6 (decimal)
#  55 '7' 20 :7 (decimal)
#  56 '8' 21 :8 (decimal)
#  57 '9' 22 :9 (decimal)
#  58 ':' 25 : __not_use__ 
#  59 ';' 25 : __not_use__
#  60 '<' 25 :-32 
#  61 '=' 25 :-31 
# :
#  90 'Z' 55 :-2 
#  91 '[' 56 :-1 
#  92 '!' 57 : 0 
#  93 ']' 58 :+1 
#  94 '^' 59 :+2 
#  :
# 124 '|' 89 ':+32

#This information is described in a folowing WEB sites:
#http://jmol.sourcearchive.com/documentation/12.0.40-1ubuntu1/classorg_1_1jmol_1_1jvxl_1_1data_1_1JvxlCoder_afc85bacfa66cd8b8d93301cf2c1dc856.html
#
#-------------------------------------------------------------------------------------------------------------------
#static boolean org::jmol::jvxl::data::JvxlCoder::appendXmlTriangleData 	( 	StringBuffer  	sb,
#		int  	triangles[][],
#		int  	nData,
#		int[]  	vertexIdNew,
#		boolean  	escapeXml 
#	) 		[inline, static, private]
#
#encode triangle data -- [ia ib ic] [ia ib ic] [ia ib ic] ... algorithm written by Bob Hanson, 11/2008. The principle is that not all vertices may be represented -- we only need the used vertices here. Capitalizing on the fact that triangle sets tend to have common edges and similar numbers for sequential triangles.
#
#a) Renumbering vertices as they appear in the triangle set
#
#[2456 2457 2458] [2456 2459 2458]
#
#becomes
#
#[ 1 2 3] [ 1 4 3]
#
#b) This allows efficient encoding of differences, not absolute numbers.
#
#0 1 2 -2 3 -1
#
#c) Which can then be represented often using a single ASCII character. I chose \ to be 0, and replace that with !.
#
#ASCII: -30 -20 -10 0 +10 +20 +30 <=>?[\]^_`abcdefghijklmnopqrstuvwxyz{|
#
#So the above sequence would simply be:
#
#!]^Z_[
#
#When the range falls outside of +/-32, we simply use a number. When a positive number follows another number, we add a "+" to it.
#
#!]^Z_[-33+250]230-210]]
#
#Preliminary trials indicated that on average a triangle can be encoded in about 7 bytes, or roughly half the 12 bytes necessary for standard binary encoding of integers. The advantage here is that we have an ASCII-readable file and no little-/big-endian issue.
#
#-------------------------------------------------------------------------------------------------------------------

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



def rot_angle_axis_to_matrix(angle,wx0,wy0,wz0):
  lenW = math.sqrt(wx0*wx0+wy0*wy0+wz0*wz0)
  wx = wx0/lenW
  wy = wy0/lenW
  wz = wz0/lenW
  s_th = math.sin(angle/180.0*math.pi)
  c_th = math.cos(angle/180.0*math.pi)
  R  = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
  S  = [[0.0,-wz, wy], [ wz,0.0,-wx], [-wy, wx,0.0]]
  SS = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
  I  = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]
  multiply_matrix3(SS,S,S)
  for i in range(3):
    for j in range(3):
      R[i][j] = I[i][j] + s_th * S[i][j] + (1.0-c_th)*SS[i][j]
  return(R)



def multiply_matrix3(C,A,B):
  for i in range(3):
    for j in range(3):
      C[i][j] = 0.0
      for k in range(3):
        C[i][j] += A[i][k] * B[k][j]


def translate_into_key_value(line,dat):
# >> Example <<
# line = '<jvxlVertexData count="4" min="{0.0 0.0 0.0}" max="{1.0 2.0 3.0}" encoding="base90xyz2" data="#~4 ||###|#~5 ||###|#"></jvxlVertexData>'
# ==> dat['count'] = '4' 
#     dat['min']='{0.0 0.0 0.0}' 
#     dat['max']='{1.0 2.0 3.0}'
#     dat['encoding']='base90xyz2'
#     dat['data']='#~4 ||###|#~5 ||###|#'
  status = 'i'
  key = ''
  value = ''
  for i in range(len(line)):
    s = line[i]
    #print "symbol '%s' status '%s' key '%s' value '%s'"%(s,status,key,value)
    if (status=='i'):
      if (s==' '):
        status = ' '
    elif (status==' '):
      if (s!=' '):
        status = 'K'
        key = s
    elif (status=='='):
      if (s=='"'):
        status = 'V'
        value = ''

    elif (status=='K'):
      if (s=='='):
        status = '='
      elif (s!=' '):
        key += s

    elif (status=='V'):
      if (s=='"'):
        if (key != ''):
          dat[key] = value
          value = ''
          key   = ''
          status = ' '
      else:
        value += s
 

def translate_base90_string_into_pos3Dlist(base90str,VertexCount,min,max):
#data="#~4 ||###|#~5 ||###|#"

  num90list =  translate_base90_string_into_num90list(base90str)
  ##  num90list[] --> pos3Dlist[][] ## 
  upp = [0,0,0]
  low = [0,0,0]
  pos3Dlist    = [[0.0,0.0,0.0] for v in range(VertexCount)]
  for v in range (VertexCount):
    for k in range(3):
      upp[k] = num90list[3*v+k]
      low[k] = num90list[3*VertexCount+3*v+k]
      num8100 = 90*upp[k]+low[k]
      pos3Dlist[v][k] = float(num8100)*(max[k]-min[k])/8099.0 + min[k]  
  return(pos3Dlist)


def translate_base90_string_into_triangle(base90str):
#data="#~4 ||###|#~5 ||###|#"
#data="!]^Z_[-33+250]230-210]]"
#[ascii] [80] 
#  43 '+'  8 : PLUS indicator 
#  45 '-' 10 : MINUS indicator 
#  48 '0' 13 :0 (decimal) 
#  49 '1' 14 :1 (decimal)
#  50 '2' 15 :2 (decimal)
# :
#  56 '8' 21 :8 (decimal)
#  57 '9' 22 :9 (decimal)
#  60 '<' 25 :-32 
#  61 '=' 26 :-31 
# :
#  90 'Z' 55 :-2 
#  91 '[' 56 :-1 
#  92 '!' 57 : 0 
#  93 ']' 58 :+1 
#  94 '^' 59 :+2 
#  :
# 124 '|' 89 ':+32

#This information is described in a folowing WEB sites:
#http://jmol.sourcearchive.com/documentation/12.0.40-1ubuntu1/classorg_1_1jmol_1_1jvxl_1_1data_1_1JvxlCoder_afc85bacfa66cd8b8d93301cf2c1dc856.html
#

  num90list =  translate_base90_string_into_num90list(base90str)

  facelist = []

  v0 = 1
  i = 0
  Nver = 0  

  while (i<len(num90list)):
    sym = base90str[i]
    nn = num90list[i]
    diffstr = base90str[i]
    if (25<=nn) and (nn<=89):  
      diff = nn - 57
    elif (sym=='+') or (sym=='-') or ((13<=nn) and (nn<=22)):
      if (sym=='+'): 
        diffstr = '+'
      elif (sym=='-'): 
        diffstr = '-'
      else:
        diffstr = '+%s'%(base90str[i])

      end = 0
      i += 1
      while (end==0):
        nn = num90list[i]
        diffstr += base90str[i] 
        if (i>=len(num90list)):
          end = 1
        elif ((base90str[i+1]=='+') or (base90str[i+1]=='-')):
          end = 1
        elif ((25<=num90list[i+1]) and (num90list[i+1]<=89)):
          end = 1 

        if (end==0):
          i += 1
      #print "diffstr '%s'"%(diffstr)
      diff = int(diffstr)

    v = v0 + diff
    #print "diffstr '%s' diff %d --> v %d"%(diffstr,diff,v)
    v0 = v

    if ((Nver%3)==0):    
      tri = [0,0,0]
      t = 0
    else:
      t += 1
    tri[t] = v 

    if ((Nver%3)==2):    
      facelist.append(tri)
 
    Nver += 1   
    i += 1
   
  return(facelist)




def translate_base90_string_into_num90list(base90str):
#data="#~4 ||###|#~5 ||###|#"
# <base 90>
# '#':35 => 0
# '$':36 => 1
# '%':37 => 2 
# :
# '\':92 -> '!' => 57
# :
# '|':124 => 89

  num90list = []

  ## (1) making character->number dic 'A2N' ##
  A2N = {}
  for i in range(35,126):
    A2N[chr(i)] = i - 35
  A2N['!'] = 57
  #print A2N

  ## (2) base90str --> num90list[] ##
  i = 0
  while (i<(len(base90str))):
    s = base90str[i]
    #print "%d '%s':%d"%(i,s,A2N.get(s,-1))
    if (s in A2N):
      num90list.append(A2N[s])
    i += 1 

  return(num90list)








def translate_pos3Dlist_into_base90_string(pos3Dlist,min,max):
  #print "#translate_pos3Dlist_into_base90_string(pos3Dlist,min %f %f %f max %f %f %f)"%(min[0],min[1],min[2],max[0],max[1],max[2])
  num8100list = []
  ### (1) float --> int with range 90*90 = 8100 ##
  for pos in (pos3Dlist):
    num8100pos = [0,0,0]
    for k in range(3):
      num8100pos[k] = int(8099*(pos[k]-min[k])/(max[k]-min[k]))
      if (num8100pos[k]>=8100):
        num8100pos[k] = 8099

    num8100list.append(num8100pos)

  ### (2) numlist --> upperlist and lowerlist
  upper90list = []
  lower90list = []

  for num8100pos in (num8100list):
    upper90 = [0,0,0]
    lower90 = [0,0,0]
    for k in range(3):
      lower90[k] = num8100pos[k]%90
      upper90[k] = (num8100pos[k]-lower90[k])/90
      #if (upper90[k]==90):
      #  print "%d num8100 %d (up %d low %d)"%(k,num8100pos[k],upper90[k],lower90[k])

    lower90list.append(lower90) 
    upper90list.append(upper90) 

  ### (3) upper90list and lower90list ==> str_base90  
# <base 90>
# 0  => '#':35
# 1  => '$':36
# 2  => '%':37 
# :
# 57 =>'\':92 --> '!':92
# :
# 89 => '|':124
  N2A = {}
  for i in range(90):
    N2A[i] = chr(i+35)
  N2A[57] = '!'

  str_base90 = ''
  for upper90 in (upper90list):
    for k in range(3):
      str_base90 += N2A[upper90[k]]
      #str_base90 += N2A.get(upper90[k],'')
  for lower90 in (lower90list):
    for k in range(3):
      str_base90 += N2A[lower90[k]]
      #str_base90 += N2A.get(lower90[k],'')
  return(str_base90)
  pass


def recover_repeat_base90_string(str):
# A~[n][space] ==> A x [n]

  newstr = '' 
  L = len(str)
  i = 0
  while (i<L):
    a = str[i]
    b = '' 
    c = '' 
    d = '' 
    if (i<(L-1)):
      b = str[i+1]
    if (i<(L-2)):
      c = str[i+2]
    if (i<(L-3)):
      d = str[i+3]
    #if (b=='~'):
    #  print "'%s%s%s%s'"%(a,b,c,d)
    if (b=='~') and (c.isdigit()):
      repstr = a 
      for k in range(int(c)):
        repstr += a
      newstr += repstr
      i += 4
    else: 
      newstr += a 
      i += 1
  return(newstr)
 
def decode_arrows_ands(str):
  str = str.replace('~;5 ','<<<<<<')
  str = str.replace('~;4 ','<<<<<')
  str = str.replace('~;3 ','<<<<')
  str = str.replace('~;2 ','<<<')
  str = str.replace('~;1 ','<<')
  str = str.replace('~;0 ','<')
  str = str.replace('~%5 ','&&&&&&')
  str = str.replace('~%4 ','&&&&&')
  str = str.replace('~%3 ','&&&&')
  str = str.replace('~%2 ','&&&')
  str = str.replace('~%1 ','&&')
  str = str.replace('~%0 ','&')
  return(str)

def encode_arrows_ands(str):
  str = str.replace('<<<<<<','~;5 ')
  str = str.replace('<<<<<', '~;4 ')
  str = str.replace('<<<<',  '~;3 ')
  str = str.replace('<<<',   '~;2 ')
  str = str.replace('<<',    '~;1 ')
  str = str.replace('<',     '~;0 ')
  str = str.replace('&&&&&&','~%5 ')
  str = str.replace('&&&&&', '~%4 ')
  str = str.replace('&&&&',  '~%3 ')
  str = str.replace('&&&',   '~%2 ')
  str = str.replace('&&',    '~%1 ')
  str = str.replace('&',     '~%0 ')
  return(str)

 
def transform_repeat_base90_string(str):
# A x [n] ==> A~[n][space]
# ex) 
# 'A'      => 'A~0 '  
# 'AA'     => 'A~1 '
# 'AAA'    => 'A~2 '
# 'AAAA'   => 'A~3 ' same length
# 'AAAAA'  => 'A~4 ' shorter !
# 'AAAAAA' => 'A~5 ' shorter !

  newstr = '' 
  L = len(str)
  i = 0
  while (i<L):
    a = str[i]
    j = 0
    while ((i+j)<L) and (str[i+j]==a):
      j += 1
    if (j>=5):
      newstr += "%s~%d "%(a,j-1)
      i += (j-1)
    else:
      newstr += a 
    i += 1
  return(newstr)


def string_to_xyz(str):
  #str = "{0.0 0.0 0.0}" =>  return((0.0,0.0,0.0))
  #str = "{1.0 2.0 3.0}" =>  return((1.0,2.0,3.0))
  str = str.lstrip('{')
  str = str.rstrip('}')
  field = str.split(' ')
  x = float(field[0])
  y = float(field[1])
  z = float(field[2])
  return((x,y,z))






def read_jvxl_file(ifname):
  #print  "#read_jvxl_file('%s')"%(ifname)
  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR:Can't open '%s'."%(ifname)
    sys.exit(1)

  f  = open(ifname)
  str_trans_jvxl = ''
  POS = []
  for line in f:
    if (line.startswith('<jvxlVertexData')==0):
      str_trans_jvxl += line 
      #of.write("%s"%(line)) 

#<jvxlVertexData count="4" min="{0.0 0.0 0.0}" max="{1.0 1.0 1.0}" encoding="base90xyz2" data="####|#|####|####|#|####|"></jvxlVertexData>
    if (line.startswith('<jvxlVertexData')):
      line = line.rstrip('\n')
      dat = {}
      translate_into_key_value(line,dat)
      #print dat
     
      ### (1) Read VertexCount,MIN[] and MAX[] ## 
      VertexCount = 0
      MIN = [0.0, 0.0, 0.0] 
      MAX = [0.0, 0.0, 0.0] 
      if (dat.get('count','')!=''):
        VertexCount = int(dat['count'])
      #print "#VertexCount %d"%(VertexCount)
  
      (MIN[0],MIN[1],MIN[2]) = string_to_xyz(dat['min'])
      (MAX[0],MAX[1],MAX[2]) = string_to_xyz(dat['max'])
      #print "#MIN %f %f %f MAX %f %f %f"%(MIN[0],MIN[1],MIN[2],MAX[0],MAX[1],MAX[2])    
      
      ### (2) Translate base90_string into POS[][] ## 
      origstr_base90 = dat.get('data','')
      str_base90 = origstr_base90
      str_base90 = recover_repeat_base90_string(str_base90)
      str_base90 = decode_arrows_ands(str_base90)
      #print "origstr_base90 '%s'"%(origstr_base90)
      #print "str_base90     '%s'"%(str_base90)
      
      poslist = translate_base90_string_into_pos3Dlist(str_base90,VertexCount,MIN,MAX)

#<jvxlTriangleData count="4" encoding="jvxltdiff" data="!]]Z^]!Z[_[["></jvxlTriangleData>
    if (line.startswith('<jvxlTriangleData')):
      line = line.rstrip('\n')
      dat = {}
      translate_into_key_value(line,dat)
      #print dat
      origstr_base90 = dat.get('data','')
      str_base90 = origstr_base90
      str_base90 = recover_repeat_base90_string(str_base90)
      str_base90 = decode_arrows_ands(str_base90)
      #print "origstr_base90 '%s'"%(origstr_base90)
      #print "str_base90     '%s'"%(str_base90)
      num90list =  translate_base90_string_into_num90list(str_base90)
      facelist = translate_base90_string_into_triangle(str_base90)
      print "#Triangle count %s"%(dat.get("count",0))
      print "#len(str_base90) %d"%(len(str_base90))
      print "#len(num90list) %d"%(len(num90list))
      print "#len(facelist) %d"%(len(facelist))
      #print num90list 
      #print facelist 
  f.close()

  return(poslist,facelist)


def write_obj_file(ofname,poslist,facelist):
  print "#write_obj_file() --> '%s'"%(ofname)
  of = open(ofname,"w")
  for pos in (poslist):
    of.write("v %f %f %f\n"%(pos[0],pos[1],pos[2]))
  for face in (facelist):
    of.write("f %d %d %d\n"%(face[0],face[1],face[2]))
  of.close()




###############
#### MAIN #####
###############
def _main():

  OPT = {}
  OPT['ijvxl'] = ''
  OPT['ojvxl'] = 'out.jvxl'
  OPT['R'] = '1:0:0:0:1:0:0:0:1'
  OPT['raxis'] = '0:1:0:0'
  OPT['T'] = '0:0:0'
  
  #str = sys.argv[1]
  #newstr = recover_repeat_base90_string(str)
  #print "'%s'-->'%s'"%(str,newstr)
  #newstr = transform_repeat_base90_string(str)
  #print "'%s'-->'%s'"%(str,newstr)
  #sys.exit(1)
  
  if (len(sys.argv)<2):
    print "jvxl.py [jvxl_file]"
    print " coded by T.Kawabata. LastModDate:%s"%(LastModDate)
    print " <options>"
    print "   -ijvxl : input jvxl file    [%s]"%(OPT['ijvxl'])
    print "   -ojvxl : output jvxl file   [%s]"%(OPT['ojvxl'])
    print "   -R     : rotation matrix    [%s]"%(OPT['R'])  
    print "   -raxis : rotation angle(degree) and axis [%s]"%(OPT['R'])  
    print "   -T     : translation vector [%s]"%(OPT['T'])  
    sys.exit(1)
  
  
  read_option(sys.argv,OPT)
  (poslist,facelist) =  read_jvxl_file(OPT['ijvxl'])

  write_obj_file("out.obj",poslist,facelist)

  sys.exit(1) 
  Rmat = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
  Tvec = [0.0,0.0,0.0]
  Tstr = OPT['T'].split(':')
  for i in range(3):
    Tvec[i] = float(Tstr[i])
 
  if (OPT['R'] != '1:0:0:0:1:0:0:0:1'):
    Rstr = OPT['R'].split(':')
    for i in range(3):
      for j in range(3):
        Rmat[i][j] = float(Rstr[3*i+j])

  if (OPT['raxis'] != '1:0:0:0:1:0:0:0:1'):
    W = OPT['raxis'].split(':')
    Rmat = rot_angle_axis_to_matrix(float(W[0]),float(W[1]),float(W[2]),float(W[3]))

 
  if (OPT['ijvxl'] != ''):
    print "#read '%s'"%(OPT['ijvxl'])
    tra_jvxl_str = read_jvxl_and_return_xyz_tranformed_jvxl(OPT['ijvxl'],Rmat,Tvec)
    if (OPT['ojvxl']!=''):
      print "#write to (T:%s R:%s) --> '%s'"%(Tvec,Rmat,OPT['ojvxl'])
      of = open(OPT['ojvxl'],'w')
      of.write("%s"%(tra_jvxl_str))
      of.close()
  #  read_jvxl_write_tranformed_xyz(OPT['ijvxl'],Rmat,Tvec,OPT['ojvxl'])


if __name__ == '__main__': _main()

 
