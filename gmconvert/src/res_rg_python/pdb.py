##
## <pdb.py>
##
## 

import sys
import os 
import math 

LastModDate = 'Aug 26, 2014'
   
class Molecule:
  
  def __init__(self):
      self.filename = ''   
      self.Natom = 0
      self.Nres  = 0
      self.anum  = []  # [0..Natom-1]
      self.atom  = []  # [0..Natom-1]
      self.resi  = []  # [0..Natom-1]
      self.rnum  = []  # [0..Natom-1]
      self.Pos   = []  # [0..Natom-1][3]
      #self.posX  = []  # [0..Natom-1]
      #self.posY  = []  # [0..Natom-1]
      #self.posZ  = []  # [0..Natom-1]
      self.chain = []  # [0..Natom-1]
      self.AHtype = [] # [0..Natom-1]
      self.line   = [] # [0..Natom-1] 'entire line for pdb'
      self.res_num = [] # [0..Natom-1] residue number 0... Nres
 
  def read(self,filename,AHtype="A",Chain="-"):
      print "#pdb.read(%s AHtype %s)"%(filename,AHtype)
  
      if not os.access(filename,os.R_OK):
        print "#ERROR:Can't open '%s'" % filename
        sys.exit()  
     
      f = open(filename)
      self.filename = filename
  
      resi0  = ''
      rnum0  = ''
      chain0 = ''
       
      self.Nres = 0 
      for line in f: 
        line = line.rstrip('\n')
    
    #          1         2         3         4         5         6         7
    #01234567890123456789012345678901234567890123456789012345678901234567890123456789
    #ATOM    676  CB  GLU A  85      10.440  29.552  12.788  6.00 16.96           C
    #ATOM    680  OE2 GLU A  85      10.230  30.451  16.374  8.00 41.03           O
    #ATOM    682  CA  LEU A  86       7.618  29.487   9.238  6.00 12.23           C
    #HETATM 1236  O4  SO4   154      33.810  28.815  -4.624  8.00 14.90           O
    #HETATM 1237 FE   HEM   155      15.271  27.962   0.622 24.00  7.86          FE
     
        if line:
          read_it = 0
          if (line[0:6]=='ATOM  ') and (AHtype=='A'):
            read_it = 1
          if (line[0:6]=='HETATM') and (AHtype=='H'):
            read_it = 1
          #if (((line[0:6]=='ATOM  ') or (line[0:6]=='HETATM')) and (AHtype=='B')):
          #  read_it = 1
          if (line[0:6]=='ATOM  ') and (AHtype=='B'):
            read_it = 1
          if (line[0:6]=='HETATM') and (AHtype=='B'):
            read_it = 1
          if (read_it==1):
            chain = line[21:22] 
            if (Chain=="-") or (chain == Chain):
              anum  = line[6:11] 
              atom  = line[12:16] 
              resi  = line[17:20] 
              chain = line[21:22] 
              rnum  = line[22:27] 
              x     = line[31:38] 
              y     = line[38:46] 
              z     = line[46:54] 
              if (line[0:6]=='ATOM  '):
                self.AHtype.append('A')
              if (line[0:6]=='HETATM'):
                self.AHtype.append('H')
            #print "'%s' '%s' '%s' '%s' '%s' %s %s %s\n" %(anum,atom,resi,chain,rnum,x,y,z)
              if (resi != resi0) or (rnum != rnum0) or (chain != chain0):
                self.Nres += 1 
              self.line.append(line)
              self.anum.append(anum)
              self.atom.append(atom)
              self.resi.append(resi)
              self.rnum.append(rnum)
              self.chain.append(chain)
              pos = [0.0, 0.0,0.0]
              pos[0] = float(x) 
              pos[1] = float(y) 
              pos[2] = float(z) 
              self.Pos.append(pos)
              #self.posX.append(float(x))
              #self.posY.append(float(y))
              #self.posZ.append(float(z))
              self.res_num.append(self.Nres)
              #sys.stdout.write("%s %d\n"%(line,self.Nres))
              resi0 = resi
              rnum0 = rnum 
              chain0 = chain 
      f.close
      self.Natom = len(self.anum)
      print "#Natom %d"%(self.Natom)

  def radius_of_gylation_DD(self):
    sumDD = 0.0
    for i in range(self.Natom):
      for j in range(self.Natom):
        dx = self.Pos[i][0] - self.Pos[j][0]
        dy = self.Pos[i][1] - self.Pos[j][1]
        dz = self.Pos[i][2] - self.Pos[j][2]
        sumDD += dx*dx + dy*dy + dz*dz 
    return(math.sqrt(sumDD/2.0/self.Natom/self.Natom))

  def radius_of_gylation_gcen(self):
    gx = 0.0 
    gy = 0.0 
    gz = 0.0 
    for i in range(self.Natom):
      gx += self.Pos[i][0]
      gy += self.Pos[i][1]
      gz += self.Pos[i][2]
    gx /= self.Natom
    gy /= self.Natom
    gz /= self.Natom

    sumDD = 0.0
    for i in range(self.Natom):
      dx = self.Pos[i][0] - gx
      dy = self.Pos[i][1] - gy
      dz = self.Pos[i][2] - gz
      sumDD += dx*dx + dy*dy + dz*dz 
    return(math.sqrt(sumDD/self.Natom))


 
  def __str__(self):
      s = '' 
      s = "Natom %d"%(self.Natom)
      return s  
  
 
def _main():
    print "woops!!"
    if (len(sys.argv)<2):
      print "#python pdb.py [pdbfile]"
      sys.exit()
    p = Molecule()
    print p
    p.read(sys.argv[1])
    print p
#    print "Rg %f"%(p.radius_of_gylation_DD())   
    print "Rg %f"%(p.radius_of_gylation_gcen())   
if __name__ == '__main__':_main()
