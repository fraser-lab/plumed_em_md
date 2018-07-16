#!/usr/bin/env python

import sys
import os
import math
import random

LastModDate = "Oct 17, 2014"


if (len(sys.argv)<3):
  print "testRg.py [Npnt] ['S'phere|'B'ox'|'G'aussian] (output pdbfile)"
  print " LastModDate:%s"%(LastModDate)
  sys.exit(1)

Npnt = int(sys.argv[1])
type = sys.argv[2][0]
opdbfile = ''
if (len(sys.argv)>3):
  opdbfile = sys.argv[3]


N = 0
sumRR = 0.0
RRg = 0.0
X = []
Y = []
Z = []

for n in range(Npnt):
  if (type=='B') or (type=='S'):
    x = 2.0*random.random() - 1.0
    y = 2.0*random.random() - 1.0
    z = 2.0*random.random() - 1.0
    RR = x*x + y*y + z*z
    if (type=='B') or ((type=='S') and (RR<=1.0)):
      sumRR += RR
      N += 1
      if (opdbfile != ''):
        X.append(x)
        Y.append(y)
        Z.append(z)
  if (type=='G'):
    U1 = random.random()
    U2 = random.random()
    x = math.sqrt(-2*math.log(U1))*math.cos(2*math.pi*U2)
    U1 = random.random()
    U2 = random.random()
    y = math.sqrt(-2*math.log(U1))*math.cos(2*math.pi*U2)
    U1 = random.random()
    U2 = random.random()
    z = math.sqrt(-2*math.log(U1))*math.cos(2*math.pi*U2)
    RR = x*x + y*y + z*z
    sumRR += RR
    N += 1
    if (opdbfile != ''):
      X.append(x)
      Y.append(y)
      Z.append(z)
    

if (N>0):
  RRg = sumRR/N

print "type %s Npnt %d N %d RRg %lf"%(type,Npnt,N,RRg)

if (opdbfile != ''):
  of = open(opdbfile,'w')
  print "output_pdb_file() --> '%s'\n"%(opdbfile)
  for i in range(N):
    of.write("ATOM  %5d %4s %3s %c%5s   %8.3f%8.3f%8.3f\n"%(i+1 ,"CA","PNT",' ',"1",100*X[i],100*Y[i],100*Z[i]))
  of.close()
