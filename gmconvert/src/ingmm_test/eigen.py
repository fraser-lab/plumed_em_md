##
## <eigen.py>
##  functions for calculating eigen values and eigen vectors 
##

import sys
import os
import math
import random

LastModDate = "Aug 3, 2014"


def cal_eigen_vector_by_Jacobi_Wilkinson(N,A,U):
  ## A and U are N x N matrices
  ##
  ##      [before]           [after] 
  ##  A : input matrix  --> diagonal matrix for eigen value
  ##  U : ------------  --> eigen vector matrix  
  ##    
  ##  Each column vector of U corresponds to the eigen vector for A.
  ## evec0 = |Umat[0][0]|  eval0 = A[0][0]
  ##         |Umat[1][0]|
  ##         |Umat[2][0]|
  ## evec1 = |Umat[0][1]|  eval1 = A[1][1]
  ##         |Umat[1][1]|
  ##         |Umat[2][1]|
  ## evec2 = |Umat[0][2]|  eval2 = A[2][2]
  ##         |Umat[1][2]|
  ##         |Umat[2][2]|


  R     = [[0.0 for i in range(N)] for j in range(N)] 
  TR    = [[0.0 for i in range(N)] for j in range(N)] 
  BUFF  = [[0.0 for i in range(N)] for j in range(N)] 
  sqrt2 = math.sqrt(2.0);

  set_identiry_matrix(N,U)
  max_mini = 0.00000001
  c_max  = 100  
  [mi,mj,max] =  find_max_abs_nondiagonal(N,A)
  c = 0;
  while ((max>max_mini) and (c<c_max)):
    #print "mi %d mj %d max %.10f"%(mi,mj,max)
    #print_matrix(3,A,"A")
    #print_matrix(3,U,"U")
    ## --------- SET of cos sin ----------- ##
    wa = (A[mi][mi] + A[mj][mj])/2.0
    sa = (A[mi][mi] - A[mj][mj])/2.0
    r  = math.sqrt(sa*sa + A[mi][mj]*A[mi][mj])
    if (sa>0.0):
      co =  math.sqrt(1.0+sa/r)/sqrt2
      si =  A[mi][mj]/2.0/r/co;
    else:
      co =  math.sqrt(1.0-sa/r)/sqrt2
      si = - A[mi][mj]/2.0/r/co

    ## -------- SET of Rot Matrix R and TR-----## 
    set_identiry_matrix(N,R)
    R[mi][mi] =  co
    R[mi][mj] = -si
    R[mj][mi] =  si
    R[mj][mj] =  co

    set_transpose_matrix(N,TR,R)

    prod_matrix(N,BUFF,A,R)
    prod_matrix(N,A,TR,BUFF)
    prod_matrix(N,BUFF,U,R)
    copy_matrix(N,U,BUFF)

    [mi,mj,max] = find_max_abs_nondiagonal(N,A)
    c+= 1




def find_max_abs_nondiagonal(N,A):   ## find maximum |A[mi][mj]| = max, and return [mi,mj,max] ## 
  max = math.fabs(A[0][1])
  mi = 0
  mj = 1
  for i in range(N):
    for j in range(N):
      if (i!=j):
        absA = math.fabs(A[i][j]) 
        if (absA>max):
          mi = i
          mj = j
          max = absA
  return([mi,mj,max])



def diagonize_CovM(CovM,eigval,eigvec):
# CovM[3][3]      : Covariance Matrix (input)
# eigval[3]     : eigen values ( to be calculated)
# eigvec[3][3] : eigen values ( to be calculated)
#
# eigval[0] : largest, g.eigval[1]:middle, g.eigval[2]:smallest
# eigvec[0:largest ][0:x, 1:y, 2:z]
# eigvec[1:middle  ][0:x, 1:y, 2:z]
# eigvec[2:smallest][0:x, 1:y, 2:z]

  A = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]
  U = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]
  y = [0.0, 0.0, 0.0]
  copy_matrix(3,A,CovM)
  cal_eigen_vector_by_Jacobi_Wilkinson(3,A,U)
  index = [0,1,2]
  sindex = sorted(index,lambda x,y:cmp(A[y][y],A[x][x]))
  for i in range(3):
    eigval[i] = A[sindex[i]][sindex[i]]
    for j in range(3):
      eigvec[i][j] = U[j][sindex[i]]


def set_MatStd_for_random_point_generation(G):
  ##  U : eigvec matrix 
  ##  L : eigval diagonal matrix 
  ##  CovM U^T = U^T L
  ##  CovM  = U^T L U
  ##  Therefore, if 
  ##   ----------------------------   
  ##     A = MatStd = U^T sqrt(L)
  ##   ----------------------------   
  ##  then 
  ##    A A^T = U^T sqrt(L)sqrt(L) U = U^T L U = CovM
  ##  
  ##  If y follws the standard Gauss distribution N(y|0,1), x = Ay + m
  ##   <(x-m)(x-m)^T> =  <(Ay)(Ay)^T> = <Ayy^TA^T> = A<yy^T>A^T  = AA^T = CovM
  ##  Therefore, x=Ay+m follows N(x|m,CovM).

  print "#set_transform_matrix_for_standard_normal_dist()"
  for g in (G.gauss):
    g.MatStd = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
    ##  g.MatStd = U^T sqrt(L) 
    for i in range(3):
      for j in range(3):
        g.MatStd[i][j] = 0.0
        for k in range(3):
          if (k==j):
            g.MatStd[i][j] += g.eigvec[k][i] * math.sqrt(g.eigval[k])
    print g.MatStd

    ## Check the calculation ##
    #CovMrecal = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
    #for i in range(3):
    #  for j in range(3):
    #    CovMrecal[i][j] = 0.0
    #    for k in range(3):
    #      CovMrecal[i][j] += g.MatStd[i][k] * g.MatStd[j][k] 
    #print "eigval %f %f %f"%(g.eigval[0],g.eigval[1],g.eigval[2])
    #eigen.print_matrix(3,g.CovM,"CovM")
    #eigen.print_matrix(3,CovMrecal,"CovMrecal")



def prod_matrix(N,C,A,B): ## C = A * B ## 
  for i in range(N):
    for j in range(N):
      C[i][j] = 0.0
      for k in range(N):
        C[i][j] += A[i][k] * B[k][j]


def prod_matrix_vec(N,y,A,x): ## y = A * x ## 
  for i in range(N):
    y[i] = 0.0
    for k in range(N):
      y[i] += A[i][k] * x[k]


def rotate_vec3_around_center(x,R,c): ## x' = R(x-c)+c ##
  xc = [0.0, 0.0, 0.0]
  y  = [0.0, 0.0, 0.0]
  for i in range(3):
    xc[i] = x[i] - c[i] 
  prod_matrix_vec(3,y,R,xc)
  for i in range(3):
    x[i] = y[i] + c[i] 




def equal_vec(N,y,x): ## y := x ## 
  for i in range(N):
    y[i] = x[i]


def equal_matrix(N,A,B): ## A := B ## 
  for i in range(N):
    for j in range(N):
      A[i][j] = B[i][j]


def prod_transposed_matrix(N,C,A,B): ## C = tA * B ## 
  for i in range(N):
    for j in range(N):
      C[i][j] = 0.0
      for k in range(N):
        C[i][j] += A[k][i] * B[k][j]

def prod_tR_matrix_R(N,C,R,A): ## C = tR A R ## 
  AR = [[0.0 for i in range(N)] for i in range(N)]

  for i in range(N):
    for j in range(N):
      AR[i][j] = 0.0
      for k in range(N):
        AR[i][j] += A[i][k] * R[k][j]

  for i in range(N):
    for j in range(N):
      C[i][j] = 0.0
      for k in range(N):
        C[i][j] += R[k][i] * AR[k][j]


def prod_R_matrix_tR(N,C,R,A): ## C = R A tR ## 
  RA = [[0.0 for i in range(N)] for i in range(N)]

  for i in range(N):
    for j in range(N):
      RA[i][j] = 0.0
      for k in range(N):
        RA[i][j] += R[i][k] * A[k][j]

  for i in range(N):
    for j in range(N):
      C[i][j] = 0.0
      for k in range(N):
        C[i][j] += RA[i][k] * R[j][k]



def copy_matrix(N,A,B): ## A = B ## 
  for i in range(N):
    for j in range(N):
      A[i][j] = B[i][j]

def set_identiry_matrix(N,A):
  for i in range(N):
    for j in range(N):
      if (i==j):
        A[i][j] = 1.0
      else:
        A[i][j] = 0.0

def set_transpose_matrix(N,TA,A):
  for i in range(N):
    for j in range(N):
      TA[j][i] = A[i][j];

def print_matrix(N,A,comment):
  print ">%s"%(comment)
  for i in range(N):
    for j in range(N):
      sys.stdout.write(" %f"%(A[i][j]))
    sys.stdout.write("\n")
  sys.stdout.write("\n")

def length_vec(N,x):
  norm = 0.0
  for i in range(N):
    norm += x[i]*x[i]
  return(math.sqrt(norm))

def Normalize_vec(N,x):
  len = length_vec(N,x)
  nx = [0.0 for i in range(N)]
  for i in range(N):
    nx[i]  = x[i]/len
  return(nx)



def change_rotation_matrix_to_axis_angle(R):
# R = I + sin_t S + (1-cos_t) S^2
# w = (a,b,c), a*a + b*b + c*c = 1
#    | 0 -c  b|       | 0  c -b|        |-b*b-c*c ab        ac     |
# S =| c  0 -a|, S^T =|-c  0  a|, S^2 = | ab      -a*a-c*c  bc     |
#    |-b  a  0|       | b -a  0|        | ac      bc       -a*a-b*b|
#
# traceR = 1 + 2*cos_t
# cos_t = (traceR-1)/2
# R - R^T = (2sin_t)S

  traceR = R[0][0] + R[1][1] + R[2][2]
  cos_t = (traceR-1.0)/2.0
  theta = math.acos(cos_t)
  sin_t = math.sqrt(1.0 - cos_t*cos_t)
  w = [1.0, 0.0, 0.0]
  S = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0] ]

  if (math.fabs(sin_t)<0.000001):
    theta = 0.0
    return(w,theta)

  for i in range(3):
    for j in range(3):
      S[i][j] = (R[i][j] - R[j][i])/(2.0*sin_t)
#  eigen.print_matrix(3,S,"S")
  w[0] = S[2][1]
  w[1] = S[0][2]
  w[2] = S[1][0]
  lenw = w[0]*w[0] + w[1]*w[1] + w[2]*w[2]
  #print "theta %f w %f %f %f lenw %f"%(theta*180.0/math.pi,w[0],w[1],w[2],lenw)
  return(w,theta)


def change_axis_angle_to_rotation_matrix(w,theta):
# R = I + sin_t S + (1-cos_t) S^2
  R  = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0] ]
  I  = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0] ]
  S  = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0] ]
  SS = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0] ]
  cos_t = math.cos(theta)
  sin_t = math.sin(theta)
  S[0][0] =  0.0
  S[0][1] = -w[2] 
  S[0][2] =  w[1] 
  S[1][0] =  w[2] 
  S[1][1] =  0.0 
  S[1][2] = -w[0]
  S[2][0] = -w[1] 
  S[2][1] =  w[0] 
  S[2][2] =  0.0
  prod_matrix(3,SS,S,S)
  for i in range(3): 
    for j in range(3): 
      R[i][j] = I[i][j] + sin_t * S[i][j] + (1.0 - cos_t) * SS[i][j]
  return(R) 



def mk_rotmatrix_from_quarternion(q):
  R = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  R[0][0] = 2*q[0]*q[0]+2*q[1]*q[1] - 1
  R[0][1] = 2*q[1]*q[2]-2*q[0]*q[3]
  R[0][2] = 2*q[1]*q[3]+2*q[0]*q[2]
  R[1][0] = 2*q[1]*q[2]+2*q[0]*q[3]
  R[1][1] = 2*q[0]*q[0]+2*q[2]*q[2] - 1
  R[1][2] = 2*q[2]*q[3]-2*q[0]*q[1]
  R[2][0] = 2*q[1]*q[3]-2*q[0]*q[2]
  R[2][1] = 2*q[2]*q[3]+2*q[0]*q[1]
  R[2][2] = 2*q[0]*q[0]+2*q[3]*q[3] - 1
  #print_matrix(3,R,"R")
  return(R)


def _main():
  if (len(sys.argv)<2):
    print "python eigen.py "
    print "  coded by T.Kawabata. LastModified:%s"%(LastModDate)
    sys.exit(1)


  q = [0.0, 0.0, 0.0, 0.0]
  qlen = 0.0
  for i in range(4):
    q[i] = random.random()
    qlen += q[i]*q[i]
  qlen = math.sqrt(qlen)
  for i in range(4):
    q[i] /= qlen

  R = mk_rotmatrix_from_quarternion(q)
  print_matrix(3,R,"R")
  (w,theta) = change_rotation_matrix_to_axis_angle(R)
  print "w",w,"theta",theta
  Ragain = change_axis_angle_to_rotation_matrix(w,theta)
  print_matrix(3,Ragain,"Ragain")


if __name__ == '__main__':_main()



