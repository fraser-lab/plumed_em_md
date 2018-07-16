#!/usr/bin/env python

import sys
import os
from datetime import datetime
import math

LastModDate = "Oct 22, 2013"

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



def rot_angle_axis_to_matrix(angle,wx,wy,wz):
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


def rot_matrix_to_angle_axis(R):
  #print "#rot_matrix_to_angle_axis()"
  traceR = R[0][0] + R[1][1] + R[2][2]
  theta = math.acos((traceR-1.0)/2.0)
  sin_th = math.sin(theta)
  S = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
  for i in range(3):
    for j in range(3):
      S[i][j] = (R[i][j] - R[j][i])/2.0/sin_th
  #print "S %f %f %f"%(S[0][0],S[0][1],S[0][2])
  #print "S %f %f %f"%(S[1][0],S[1][1],S[1][2])
  #print "S %f %f %f"%(S[2][0],S[2][1],S[2][2])
  wx = S[2][1]
  wy = S[0][2]
  wz = S[1][0]
  #print "wx %f wy %f wz %f"%(wx,wy,wz)
  len_w = math.sqrt(wx*wx + wy*wy + wz*wz)
  wx /= len_w
  wy /= len_w
  wz /= len_w
  theta = theta/math.pi * 180.0
  return(theta,wx,wy,wz)

def multiply_matrix3(C,A,B):
  for i in range(3):
    for j in range(3):
      C[i][j] = 0.0
      for k in range(3):
        C[i][j] += A[i][k] * B[k][j]


def multiply_matrix3_AxTrB(C,A,B):
  for i in range(3):
    for j in range(3):
      C[i][j] = 0.0
      for k in range(3):
        C[i][j] += A[i][k] * B[j][k]

def cross_vec3D(c,a,b): ## c = a x b ## 
  c[0] = a[1]*b[2] - b[1]*a[2];
  c[1] = a[2]*b[0] - b[2]*a[0];
  c[2] = a[0]*b[1] - b[0]*a[1];

def dprod_vec3D(a,b):
  return(a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

def normalize_vec3D(a):
  len = math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
  a[0] /= len
  a[1] /= len
  a[2] /= len


def read_GMM_file(ifname,gdflist):
  print "#read_GMM_file('%s')"%(ifname)
  
#REMARK NGAUSS 2
#REMARK GAUSS   1 W 0.5169505265
#REMARK GAUSS   1 M 40.249309 -77.808283 21.467092
#REMARK GAUSS   1 CovM  xx  119.0299750247 xy   30.5634422043 xz   22.8729372627
#REMARK GAUSS   1 CovM  yy  172.7122406834 yz  -31.1567962479 zz   80.8062596679
#REMARK GAUSS   1 PCvar   190.283289 126.947102 55.318085
#REMARK GAUSS   1 PCaxis1 0.333874 0.922682 -0.192836
#REMARK GAUSS   1 PCaxis2 0.824277 -0.186542 0.534574
#REMARK GAUSS   1 PCaxis3 0.457270 -0.337431 -0.822827
#REMARK GAUSS   2 W 0.4830494735
#REMARK GAUSS   2 M 40.504943 -76.504505 56.432250
#REMARK GAUSS   2 CovM  xx  131.5632484805 xy   10.6845806101 xz   -5.2276828546
#REMARK GAUSS   2 CovM  yy  147.5195516343 yz   55.8430262458 zz   84.1230814480
#REMARK GAUSS   2 PCvar   180.930583 131.877349 50.397950
#REMARK GAUSS   2 PCaxis1 0.134594 0.861461 0.489662
#REMARK GAUSS   2 PCaxis2 0.983424 -0.055540 -0.172604
#REMARK GAUSS   2 PCaxis3 -0.121495 0.504777 -0.854657
  if (os.access(ifname,os.R_OK)==0):
    print "#ERROR:Can't open gmmfile '%s'."%(ifname)
    sys.exit(1)
  f = open(ifname,'r')  
  gnum0 = 0
  for line in f:
    if (line.startswith('REMARK GAUSS')):
      line = line.rstrip('\n')
      field =  line.split()
      gnum = int(field[2])
      type = field[3]
      if (gnum != gnum0):
        if (gnum0 != 0):
          gdflist.append(gdf) 
        gdf = {}

      if (type=='W'):
        gdf['W'] = float(field[4]) 
      elif (type=='M'):
        gdf['M'] = [0.0,0.0,0.0]
        gdf['M'][0] = float(field[4])
        gdf['M'][1] = float(field[5])
        gdf['M'][2] = float(field[6])
      elif (type=='CovM') and (field[4]=='xx'):
        gdf['CovM'] = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
        gdf['CovM'][0][0] = float(field[5])
        gdf['CovM'][0][1] = gdf['CovM'][1][0] = float(field[7])
        gdf['CovM'][0][2] = gdf['CovM'][2][0] = float(field[9])
      elif (type=='CovM') and (field[4]=='yy'):
        gdf['CovM'][1][1] = float(field[5])
        gdf['CovM'][1][2] = gdf['CovM'][2][1] = float(field[7])
        gdf['CovM'][2][2] = float(field[9])
      elif (type=='PCvar'):
        gdf['PCvar'] = [0.0,0.0,0.0]
        gdf['PCvar'][0] = float(field[4])
        gdf['PCvar'][1] = float(field[5])
        gdf['PCvar'][2] = float(field[6])
      elif (type=='PCaxis1'):
        gdf['PCaxis1'] = [0.0,0.0,0.0]
        gdf['PCaxis1'][0] = float(field[4])
        gdf['PCaxis1'][1] = float(field[5])
        gdf['PCaxis1'][2] = float(field[6])
      elif (type=='PCaxis2'):
        gdf['PCaxis2'] = [0.0,0.0,0.0]
        gdf['PCaxis2'][0] = float(field[4])
        gdf['PCaxis2'][1] = float(field[5])
        gdf['PCaxis2'][2] = float(field[6])
      elif (type=='PCaxis3'):
        gdf['PCaxis3'] = [0.0,0.0,0.0]
        gdf['PCaxis3'][0] = float(field[4])
        gdf['PCaxis3'][1] = float(field[5])
        gdf['PCaxis3'][2] = float(field[6])

      gnum0 = gnum

  if (gnum0 != 0):
    gdflist.append(gdf) 
  f.close()

def set_Vmat_from_PCaxis123(V,gdf,rottype):
  if (rottype==0):
    sign = [1,1,1]
  if (rottype==1):
    sign = [1,-1,-1]
  if (rottype==2):
    sign = [-1,1,-1]
  if (rottype==3):
    sign = [-1,-1,1]

  V[0][0] = sign[0]*gdf['PCaxis1'][0]
  V[1][0] = sign[0]*gdf['PCaxis1'][1]
  V[2][0] = sign[0]*gdf['PCaxis1'][2]

  V[0][1] = sign[1]*gdf['PCaxis2'][0]
  V[1][1] = sign[1]*gdf['PCaxis2'][1]
  V[2][1] = sign[1]*gdf['PCaxis2'][2]

  V[0][2] = sign[2]*gdf['PCaxis3'][0]
  V[1][2] = sign[2]*gdf['PCaxis3'][1]
  V[2][2] = sign[2]*gdf['PCaxis3'][2]


def show_matrix3(M,title):
  print "#%s %f %f %f"%(title,M[0][0],M[0][1],M[0][2]) 
  print "#%s %f %f %f"%(title,M[1][0],M[1][1],M[1][2]) 
  print "#%s %f %f %f"%(title,M[2][0],M[2][1],M[2][2]) 


def output_displace_vector_in_pdb(ofname,displace_list,anum_offset):
  print "#output_displace_vector_in_pdb() --> '%s'"%(ofname)
  of = open(ofname,'w')
  of.write("REMARK  PDB FILE FOR DISPLACEMENT POINTS\n")
  of.write("REMARK  COMMAND %s\n"%(OPT['COMMAND']))
  of.write("REMARK  DATE    %s\n"%(OPT['START_DATE']))
  a = anum_offset
  g = 0
  Ngauss = len(displace_list)
  for disp in (displace_list):
    g += 1
    a += 1
    of.write("HETATM%5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(a,"STA","DIS","X",g, disp['sta'][0],disp['sta'][1], disp['sta'][2],disp['llen'],disp['llen']))
    a += 1
    of.write("HETATM%5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(a,"END","DIS","X",g, disp['end'][0],disp['end'][1], disp['end'][2],disp['tlen'],disp['tlen']))
    a += 1
    of.write("HETATM%5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(a,"RAX","DIS","X",g, disp['rax'][0],disp['rax'][1], disp['rax'][2],disp['rangle'],disp['rangle']))
    a += 1
    of.write("HETATM%5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(a,"RST","DIS","X",g, disp['rst'][0],disp['rst'][1], disp['rst'][2],disp['rangle'],disp['rangle']))
    a += 1
    of.write("HETATM%5d %4s %3s %1s %3d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(a,"REN","DIS","X",g, disp['ren'][0],disp['ren'][1], disp['ren'][2],disp['rangle'],disp['rangle']))
 
#      an->Anum,an->Atom,an->Resi,chainID,an->Rnum,
#      an->Pos[0],an->Pos[1],an->Pos[2],an->R,an->tFactor);
#      if (an->AHtype=='A')  fprintf(fp,"ATOM  ");
#  else if (an->AHtype=='H')  fprintf(fp,"HETATM");
#  else fprintf(fp,"ATOM? ");
# fprintf(fp,"%5s %4s %3s %c%5s   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
#      an->Anum,an->Atom,an->Resi,chainID,an->Rnum,
#      an->Pos[0],an->Pos[1],an->Pos[2],an->R,an->tFactor);
  for g in range(Ngauss):
    of.write("CONECT %4d %4d\n"%(anum_offset+5*g+1,anum_offset+5*g+2))
    of.write("CONECT %4d %4d\n"%(anum_offset+5*g+2,anum_offset+5*g+3))
    of.write("CONECT %4d %4d %4d\n"%(anum_offset+5*g+3,anum_offset+5*g+4,anum_offset+5*g+5))


 
  of.close() 

##############
#### MAIN ####
##############

OPT = {}

OPT['A'] = ''
OPT['B'] = ''
OPT['oa2b'] = ''

if (len(sys.argv)<2):
  print "cmp_pairgmm.py <options>"
  print " for comparison of a pair of GMM with same number of gdf."
  print " coded by T.Kawabata. LastModified:%s"%(LastModDate)
  print "<options>"
  print " -A: input file for GMM for molecule A [%s]"%(OPT['A'])
  print " -B: input file for GMM for molecule B [%s]"%(OPT['B'])
  print " -oa2b:output PDB file for a2b rotation[%s]"%(OPT['oa2b']) 
  sys.exit(1)

read_option(sys.argv,OPT)

gmmA = []
read_GMM_file(OPT['A'],gmmA)
print "#NgdfA %d"%(len(gmmA))

gmmB = []
read_GMM_file(OPT['B'],gmmB)
print "#NgdfB %d"%(len(gmmB))

if (len(gmmA)!=len(gmmB)):
  print "#ERROR:NgdfA %d is not equal to NgdfB %d."%(len(gmmA), len(gmmB))  
  sys.exit(1)

VmatA = [ [0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
VmatB = [ [0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
Rmat  = [ [0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
tvec  = [0.0,0.0,0.0]
lvec  = [0.0,0.0,0.0]
displace_list = []

for g in range(len(gmmA)):
  gA = gmmA[g]
  gB = gmmB[g]

  lvec[0] = math.sqrt(gB['PCvar'][0]) - math.sqrt(gA['PCvar'][0])
  lvec[1] = math.sqrt(gB['PCvar'][1]) - math.sqrt(gA['PCvar'][1])
  lvec[2] = math.sqrt(gB['PCvar'][2]) - math.sqrt(gA['PCvar'][2])
  llen = math.sqrt(lvec[0]*lvec[0] + lvec[1]*lvec[1] + lvec[2]*lvec[2]) 

  tvec[0] = gB['M'][0] - gA['M'][0]
  tvec[1] = gB['M'][1] - gA['M'][1]
  tvec[2] = gB['M'][2] - gA['M'][2]
  tlen = math.sqrt(tvec[0]*tvec[0] + tvec[1]*tvec[1] + tvec[2]*tvec[2]) 
  set_Vmat_from_PCaxis123(VmatA,gA,0)
  rangle = 1000
  raxis  = [0.0,0.0,0.0] 
  for i in range(4):
    set_Vmat_from_PCaxis123(VmatB,gB,i)
    multiply_matrix3_AxTrB(Rmat,VmatB,VmatA)
    (theta,rx,ry,rz) = rot_matrix_to_angle_axis(Rmat)
    #print "#permu_%d theta %f rot_axis (%f %f %f)"%(i,theta,rx,ry,rz)
    if (math.fabs(theta) < math.fabs(rangle)):
      rangle = theta
      raxis[0] = rx 
      raxis[1] = ry 
      raxis[2] = rz
  print "#[%d] lvec %+6.3f (%+6.3f %+6.3f %+6.3f) tvec %+6.3f (%+6.3f %+6.3f %+6.3f) rot angle %+6.3f axis (%+6.3f %+6.3f %+6.3f)"%(g,llen,lvec[0],lvec[1],lvec[2],tlen,tvec[0],tvec[1],tvec[2],rangle,raxis[0],raxis[1],raxis[2])
  dAB = [0.0, 0.0, 0.0]
  nX  = [0.0, 0.0, 0.0]
  nY  = [0.0, 0.0, 0.0]
  for i in range(3):
    dAB[i] = gB['M'][i] - gA['M'][i] 
  for i in range(3):
    nX[i] = dAB[i] - dprod_vec3D(dAB,raxis)*raxis[i] 
  normalize_vec3D(nX)
  cross_vec3D(nY,raxis,nX)
  normalize_vec3D(nY)
  #print "#DPROD raxis nX %f"%(dprod_vec3D(raxis,nX))
  #print "#DPROD raxis nY %f"%(dprod_vec3D(raxis,nY))
  #print "#DPROD nX    nY %f"%(dprod_vec3D(nX,nY))
  disp = {}
  disp['sta'] = [0.0, 0.0, 0.0]
  disp['end'] = [0.0, 0.0, 0.0]
  disp['rax'] = [0.0, 0.0, 0.0]
  disp['rst'] = [0.0, 0.0, 0.0]
  disp['ren'] = [0.0, 0.0, 0.0]
  disp['llen'] = llen 
  disp['tlen'] = tlen 
  disp['rangle'] =  rangle 
  raxis_len_disp = 5.0 
  raxis_len_rot  = math.sqrt(gA['PCvar'][0])

  for i in range(3):
    disp['sta'][i] = gA['M'][i] 
    disp['end'][i] = gB['M'][i] 
    disp['rax'][i] = gB['M'][i] + raxis_len_disp * raxis[i]
    disp['rst'][i] = gB['M'][i] + raxis_len_disp * raxis[i] + raxis_len_rot * (nX[i])
    rang_rad =  rangle/180.0*math.pi
    #rang_rad =  90.0/180.0*math.pi
    disp['ren'][i] = gB['M'][i] + raxis_len_disp * raxis[i] + raxis_len_rot * (nX[i]*math.cos(rang_rad) + nY[i]*math.sin(rang_rad))
    #disp['ren'][i] = gB['M'][i] + raxis_len_disp * raxis[i] + raxis_len_rot * (nY[i])
 
  displace_list.append(disp)


if (OPT['oa2b']!=''):
  output_displace_vector_in_pdb(OPT['oa2b'],displace_list,len(gmmA)*7)
