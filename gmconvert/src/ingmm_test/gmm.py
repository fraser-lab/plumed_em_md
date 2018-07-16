import sys
import os
import math
import eigen


LastModDate = 'July 28, 2014'

class GAUSS3D:

  def __init__(self):
    self.num = 0
    self.M = [0.0, 0.0, 0.0]  ## Mean Position ##
    self.CovM  = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]] ## Covariance Matrix (Sigma) ##
    self.iCovM = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]] ## Inverce of Sigma Matrix  ##
    self.det  = 1.0  ## Determinant of CovM ##
    self.Cons = 1.0  ## 1/{2pi)**3/2 * sqrt(det)} ##
    self.Weight = 1.0 ## Weight for Gaussian Mixture ##


  def pdf(self,X): 

    D = [0.0, 0.0, 0.0]
    D[0] = X[0] - self.M[0];
    D[1] = X[1] - self.M[1];
    D[2] = X[2] - self.M[2];
    xSx = 0.0;
    xSx += D[0]*(self.iCovM[0][0]*D[0] + self.iCovM[0][1]*D[1] + self.iCovM[0][2]*D[2]);
    xSx += D[1]*(self.iCovM[1][0]*D[0] + self.iCovM[1][1]*D[1] + self.iCovM[1][2]*D[2]);
    xSx += D[2]*(self.iCovM[2][0]*D[0] + self.iCovM[2][1]*D[1] + self.iCovM[2][2]*D[2]);
    val = self.Cons * math.exp(-0.5*xSx);
    return(val);
  
  def copy(self,gsource):
    self.num    = gsource.num
    self.det    = gsource.det
    self.Cons   = gsource.Cons
    self.Weight = gsource.Weight

    for i in range(3):
      self.M[i] = gsource.M[i] 
    for i in range(3):
      for j in range(3):
        self.CovM[i][j]  = gsource.CovM[i][j] 
        self.iCovM[i][j] = gsource.iCovM[i][j] 


  def __str__(self):
      s = ''
      s = "GAUSS3D(num %d M %f %f %f)"%(self.num,self.M[0], self.M[1], self.M[2])
      return s





class GAUSS_MOLECULE:
  def __init__(self):
    self.ifname = ''    ## File name ## 
    self.Ngauss = 0     ## Number of gaussians ## 
    self.gauss  = []    ## list of class GAUSS3D  (malloc later) ## 
    self.Weight = 1.0   ## Weight for the Gauss molecule ##
    self.M = [0.0, 0.0, 0.0]  ## Mean Position og GAUSS_MOLECULE ##
    self.CovM  = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]] ## Covariance Matrix (Sigma) of GAUSS_MOLECULE ##

  def read(self,ifname):
    if (os.access(ifname,os.R_OK)==0):
      print "#ERROR:Can't open '%s'."%(ifname)
      sys.exit(1)
    #print "#GAUSS_MOLECULE.read('%s')"%(ifname)
    f = open(ifname)

#HEADER 3D Gaussian Mixture Model
#REMARK COMMAND gmconvert -imap emd_2155.map -ng 10 -ogmm emd_2155_ng10.gmm
#REMARK START_DATE May 16,2014 9:5:26
#REMARK END_DATE   May 16,2014 9:5:48
#REMARK COMP_TIME_SEC  21.734008 2.173401e+01
#REMARK FILENAME emd_2155_ng10.gmm
#REMARK NGAUSS 10
#REMARK COMMENT Corr.Coeff. 0.910096
#HETATM    1  GAU GAU I   1     234.627 251.790 271.785 0.108 0.108
#REMARK GAUSS   1 W 0.1081987630
#REMARK GAUSS   1 M 234.626781 251.789643 271.785209
#REMARK GAUSS   1 CovM  xx  105.7822093213 xy  -34.5835273670 xz    3.0064228373
#REMARK GAUSS   1 CovM  yy  191.9283542354 yz   84.4964299230 zz  185.3724771412
#HETATM    2  GAU GAU I   2     264.409 316.186 275.657 0.102 0.102
#REMARK GAUSS   2 W 0.1017944528
#REMARK GAUSS   2 M 264.408933 316.186174 275.656658
#REMARK GAUSS   2 CovM  xx  261.2889015181 xy  -41.5952633176 xz   77.8715172938
#REMARK GAUSS   2 CovM  yy   78.2254798714 yz   -0.0826440874 zz  186.4984286734

    for line in f:
      if (len(line)>1) and (line.startswith('#')==0):
        line = line.rstrip('\n')
        if (line.startswith('REMARK NGAUSS')):
          field = line.split()
          self.Ngauss = int(field[2])
          for g in range(self.Ngauss):
            self.gauss.append(GAUSS3D())
        if (line.startswith('REMARK GAUSS ')):
          #print line
          field = line.split()
          n = int(field[2])
          #print "n %d"%(n)
          if (1<=n) and (n<=self.Ngauss):
            self.gauss[n-1].num = n
            if (field[3]=='W'):    
              self.gauss[n-1].Weight = float(field[4])
            if (field[3]=='M'):    
              self.gauss[n-1].M[0] = float(field[4])
              self.gauss[n-1].M[1] = float(field[5])
              self.gauss[n-1].M[2] = float(field[6])
            if (field[3]=='CovM'):    
              if (field[4]=='xx'):    
                self.gauss[n-1].CovM[0][0] = float(field[5])
                self.gauss[n-1].CovM[0][1] = float(field[7])
                self.gauss[n-1].CovM[0][2] = float(field[9])
                self.gauss[n-1].CovM[1][0] = self.gauss[n-1].CovM[0][1]
                self.gauss[n-1].CovM[2][0] = self.gauss[n-1].CovM[0][2]
              if (field[4]=='yy'):    
                self.gauss[n-1].CovM[1][1] = float(field[5])
                self.gauss[n-1].CovM[1][2] = float(field[7])
                self.gauss[n-1].CovM[2][2] = float(field[9])
                self.gauss[n-1].CovM[2][1] = self.gauss[n-1].CovM[1][2]
 
    f.close()


    for g in (self.gauss):
      g.det =  cal_inverse_matrix3D_by_Cramer_rule(g.iCovM,g.CovM)
      g.Cons = math.sqrt(g.det)
      g.Cons = 1.0/math.pow(2.0*math.pi,3.0/2.0) * math.sqrt(g.det)

    self.set_M_CovM()


  def write(self,ofname):
    print "#GAUSS_MOLECULE.write()-->'%s'"%(ofname)
    of = open(ofname,"w")
    of.write("REMARK NGAUSS %d\n"%(self.Ngauss))
    for g in (self.gauss):
      of.write("HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n"%(g.num,"GAU","GAU",'A', g.num, g.M[0], g.M[1],g.M[2], g.Weight,g.Weight))
      of.write("REMARK GAUSS%4d W %.10lf\n"%(g.num,g.Weight))
      of.write("REMARK GAUSS%4d M %lf %lf %lf\n"%(g.num,g.M[0], g.M[1], g.M[2]))
      of.write("REMARK GAUSS%4d CovM  xx %15.10lf xy %15.10lf xz %15.10lf\n"%(g.num,g.CovM[0][0],g.CovM[0][1],g.CovM[0][2]))
      of.write("REMARK GAUSS%4d CovM  yy %15.10lf yz %15.10lf zz %15.10lf\n"%(g.num,g.CovM[1][1],g.CovM[1][2],g.CovM[2][2]))
    of.write("TER\n")
    of.close()

  def write_as_one_gdf(self,ofname):
    print "#GAUSS_MOLECULE.write_as_one_gdf()-->'%s'"%(ofname)
    of = open(ofname,"w")
    of.write("REMARK NGAUSS 1\n")
    of.write("HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n"%(1,"GAU","GAU",'A', 1, self.M[0], self.M[1],self.M[2],1.0,1.0))
    of.write("REMARK GAUSS%4d W %.10lf\n"%(1,1.0))
    of.write("REMARK GAUSS%4d M %lf %lf %lf\n"%(1,self.M[0], self.M[1], self.M[2]))
    of.write("REMARK GAUSS%4d CovM  xx %15.10lf xy %15.10lf xz %15.10lf\n"%(1,self.CovM[0][0],self.CovM[0][1],self.CovM[0][2]))
    of.write("REMARK GAUSS%4d CovM  yy %15.10lf yz %15.10lf zz %15.10lf\n"%(1,self.CovM[1][1],self.CovM[1][2],self.CovM[2][2]))
    of.write("TER\n")
    of.close()


  def copy(self,Gsource):
    self.gauss = []
    for gs in (Gsource.gauss):
      newg = GAUSS3D()
      newg.copy(gs)
      self.gauss.append(newg)
    self.Ngauss = Gsource.Ngauss
    self.Weight = Gsource.Weight


  def set_M_CovM(self):
    ## Initialize
    for i in range(3):
      self.M[i] = 0.0
      for j in range(3):
        self.CovM[i][j] = 0.0


    ## M =  \sum_(g) { g.W * g.M  }
    for g in (self.gauss):
      for i in range(3):
        self.M[i] += g.Weight * g.M[i]
  
    ## Sigma =  \sum_(g) { g.W * (g.CovM + g.M*g.M^T) } - M*M^{T}
    for g in (self.gauss):
      for i in range(3):
        for j in range(3):
          self.CovM[i][j] += g.Weight * (g.CovM[i][j] + g.M[i]*g.M[j]) 

    for i in range(3):
      for j in range(3):
        self.CovM[i][j] = self.CovM[i][j] -  self.M[i] * self.M[j] 

    #print "#M %f %f %f"%(self.M[0], self.M[1], self.M[2])
    #print "#CovM0 %f %f %f"%(self.CovM[0][0], self.CovM[0][1], self.CovM[0][2])
    #print "#CovM1 %f %f %f"%(self.CovM[1][0], self.CovM[1][1], self.CovM[1][2])
    #print "#CovM2 %f %f %f"%(self.CovM[2][0], self.CovM[2][1], self.CovM[2][2])



  def __str__(self):
      s = ''
      s = "GAUSS_MOLECULE(Ngauss %d)"%(self.Ngauss)
      return s





def cal_inverse_matrix3D_by_Cramer_rule(invA,A):
# <Cramer's Rule>
#
# inv[A] = 1/|A| transpose[Delta_ij]
#
# Delta_ij = (-1)**(i+j) * |Aij|     --> Yoinshi
#
# Aij = matrix removing i-th row and j-th column.
#
# --> return(Determinant of A)

#  det  = 0.0
#  det +=  (  A[0][0]*A[1][1]*A[2][2] + A[0][2]*A[1][0]*A[2][1])
#  det +=  (  A[0][1]*A[1][2]*A[2][0] - A[0][2]*A[1][1]*A[2][0])
#  det +=  (- A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2])
  det =  determinant_matrix3D(A)

  invA[0][0] =  (A[1][1]*A[2][2] - A[1][2]*A[2][1])
  invA[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1])
  invA[0][2] =  (A[0][1]*A[1][2] - A[0][2]*A[1][1])

  invA[1][0] = -(A[1][0]*A[2][2] - A[1][2]*A[2][0])
  invA[1][1] =  (A[0][0]*A[2][2] - A[0][2]*A[2][0])
  invA[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0])

  invA[2][0] =  (A[1][0]*A[2][1] - A[1][1]*A[2][0])
  invA[2][1] = -(A[0][0]*A[2][1] - A[0][1]*A[2][0])
  invA[2][2] =  (A[0][0]*A[1][1] - A[0][1]*A[1][0])

  for i in range(3):
    for j in range(3):
      invA[i][j] /= det
 
  return(det)


def determinant_matrix3D(A):
  det  = 0.0
  det +=  (  A[0][0]*A[1][1]*A[2][2] + A[0][2]*A[1][0]*A[2][1])
  det +=  (  A[0][1]*A[1][2]*A[2][0] - A[0][2]*A[1][1]*A[2][0])
  det +=  (- A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2])
  return(det)
 
def add_matrix3D(C,A,B):
  for i in range(3):
    for j in range(3):
      C[i][j] = A[i][j] + B[i][j]
  pass


def prod_matrix3D(C,A,B):  ## C = A*B
  for i in range(3):
    for j in range(3):
      C[i][j] = 0.0
      for k in range(3):
        C[i][j] += A[i][k] * B[k][j]


def print_matrix3D(A,title):
  print "#>%s"%(title)
  print "#%f %f %f"%(A[0][0],A[0][1],A[0][2])
  print "#%f %f %f"%(A[1][0],A[1][1],A[1][2])
  print "#%f %f %f"%(A[2][0],A[2][1],A[2][2])


def overlap_GAUSS3D(gA,gB):
  CovMAB  = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
  iCovMAB = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
  bmat    = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
  MAB = [0.0,0.0,0.0]
  add_matrix3D(CovMAB,gA.CovM, gB.CovM)
  detAB = cal_inverse_matrix3D_by_Cramer_rule(iCovMAB,CovMAB)
  for i in range(3):
    MAB[i] = gA.M[i] - gB.M[i]
  xSx = 0.0
  xSx += MAB[0]*(iCovMAB[0][0]*MAB[0] + iCovMAB[0][1]*MAB[1] + iCovMAB[0][2]*MAB[2])
  xSx += MAB[1]*(iCovMAB[1][0]*MAB[0] + iCovMAB[1][1]*MAB[1] + iCovMAB[1][2]*MAB[2])
  xSx += MAB[2]*(iCovMAB[2][0]*MAB[0] + iCovMAB[2][1]*MAB[1] + iCovMAB[2][2]*MAB[2])
  ov = 1.0/math.pow(2*math.pi,3.0/2.0)/math.sqrt(detAB) * math.exp(-0.5*xSx)
  # print "#overlap gA %f %f %f gB %f %f %f ov %e"%(gA.M[0],gA.M[1],gA.M[2],gB.M[0],gB.M[1],gB.M[2],ov) 
  return(ov)


def corr_coeff_bwn_two_GAUSS_MOLECULE(GA,GB):
  VAA = 0.0
  VBB = 0.0
  VAB = 0.0

  for gi in (GA.gauss):
    for gj in (GA.gauss):
      VAA += gi.Weight * gj.Weight * overlap_GAUSS3D(gi,gj)

  for gi in (GB.gauss):
    for gj in (GB.gauss):
      VBB += gi.Weight * gj.Weight * overlap_GAUSS3D(gi,gj)

  for gi in (GA.gauss):
    for gj in (GB.gauss):
      VAB += gi.Weight * gj.Weight * overlap_GAUSS3D(gi,gj)

  CC = VAB/math.sqrt(VAA)/math.sqrt(VBB)
  #print "#VAB %e VAA %e VBB %e CC %f"%(VAB,VAA,VBB,CC)
  return(CC)



def diagonize_CovM_of_GAUSS_MOLECULE(G):
  G.eigval = [0.0,0.0,0.0]
  G.eigvec = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
  eigen.diagonize_CovM(G.CovM,G.eigval,G.eigvec)

  for g in (G.gauss):
    g.eigval = [0.0,0.0,0.0]
    g.eigvec = [[0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
    eigen.diagonize_CovM(g.CovM,g.eigval,g.eigvec)
    #print g.eigval,g.eigvec


def read_options(argv,opt_dic):
  opt_dic['COMMAND'] = argv[0]
  for i in range(1,len(argv)):
    opt_dic['COMMAND'] += ' '
    opt_dic['COMMAND'] += argv[i]
    if (argv[i][0]=='-'):
      if ((i+1)<len(argv)):
        if (argv[i+1][0]!='-'):
          opt_dic[argv[i][1:]] = argv[i+1]




def _main():
  OPT = {}
  OPT['oeigpdb'] = ''

  if (len(sys.argv)<2):
    print "python gmm.py [input_gmm_file] <options>"
    print "  coded by T.Kawabata. LastModified:%s"%(LastModDate)
    print "<options>"
    print " -oeigpdb : output pdb file for eigen vectors [%s]"%(OPT['oeigpdb'])
    sys.exit(1)

  read_options(sys.argv,OPT)

  G = GAUSS_MOLECULE()
  G.read(sys.argv[1])
  print G
  for g in (G.gauss):
    print g
  G.write("out.gmm")
  G.write_as_one_gdf("one.gmm")
  diagonize_CovM_of_GAUSS_MOLECULE(G)
  print "#G.eigval %f %f %f"%(G.eigval[0], G.eigval[1], G.eigval[2])
  if (OPT['oeigpdb'] != ''):
    print "#write_eigen_vectors_in_pdb --> '%s'"%(OPT['oeigpdb'])
    of = open(OPT['oeigpdb'],'w')
    of.write("HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n"%(1001,"C  ","AXS",'O', 1,G.M[0], G.M[1],G.M[2], 0.0,0.0))
  
    ax = [ [0.0,0.0,0.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
 
    for i in range(3):
      for j in range(3):
        ax[i][j] = G.M[j] + math.sqrt(G.eigval[i])*G.eigvec[i][j]

    of.write("HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n"%(1002,"N  ","AXS",'O', 2, ax[0][0],ax[0][1],ax[0][2],0.0,0.0))
    of.write("HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n"%(1003,"O  ","AXS",'O', 3, ax[1][0],ax[1][1],ax[1][2],0.0,0.0))
    of.write("HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f%6.3f%6.3f\n"%(1004,"S  ","AXS",'O', 4, ax[2][0],ax[2][1],ax[2][2],0.0,0.0))
    of.write("CONECT 1001 1002 1003 1004\n")
    pass


if __name__ == '__main__':_main()
