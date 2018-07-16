/*

 <gmconvert.c>

 for getting Gaussian Mixture Model from atoms /3D density map

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdbool.h>
#include "TabExpMinX.h"
#include "globalvar.h"
#include "pdbstruct.h"
#include "PdbIO.h"
#include "Radius.h"
#include "gauss.h"
#include "GaussIO.h"
#include "Matrix3D.h"
#include "Voxel.h"
#include "MCubeFunc.h"
#include "MCubeIO.h"
#include "MCubeVolSur.h"
#include "mc_verface.h"
#include "PointEM.h"
#include "GmmAtomEM.h"
#include "AtomKmean.h"
#include "MapCCP4.h"
#include "3DmapEM.h"
#include "Gmm3DmapEM.h"
#include "3DmapKmean.h"
#include "Atom2Vox.h"
#include "GridMap.h"
#include "BasicIO.h"
#include "qRMS.h"
#include "PDB_from_mmCIF.h"
#include "ATOMs_gmconvert_from_PDBmmCIF.h"
#include "Ellipsoid.h"

struct PARAMETERS    PAR;
struct TAB_EXP_MIN_X TABEXP;

/*** Functions (LOCAL) ***/

int main(argc,argv)
 int argc;
 char **argv;
{
 int   i,k,g,n,L,emok,x,y,z; 
 char   MODE[8];
 char ipdbfile[MAX_FILENAME],imapfile[MAX_FILENAME],igmmfile[MAX_FILENAME];
 char iciffile[MAX_FILENAME],assembly_id[MAX_FILENAME];
 char corename[MAX_FILENAME];
 char  ogmmfile[MAX_FILENAME],opcgmmfile[MAX_FILENAME],ovoxfile[MAX_FILENAME],omapfile[MAX_FILENAME],oimapfile[MAX_FILENAME],oobjfile[MAX_FILENAME];
 char  ogmmfile_init[MAX_FILENAME]; 
 char  iobjfileA[MAX_FILENAME],iobjfileB[MAX_FILENAME];
 char  owrlfile[MAX_FILENAME],opdbfile[MAX_FILENAME],oradwgtfile[MAX_FILENAME],owrlfile_ellipsoid[MAX_FILENAME];
 char  igmmfile2[MAX_FILENAME],itarpdbfile[MAX_FILENAME];
 char  ompdbfile[MAX_FILENAME], impdbfile[MAX_FILENAME];
 char  GMMlistfile[MAX_FILENAME], PDBlistfile[MAX_FILENAME];
 char  ERRlistfile[MAX_FILENAME], omaperrfile[MAX_FILENAME];
 char  ovfile[MAX_FILENAME], errfile[MAX_FILENAME];
 FILE  *ptr_file, *fp_ovfile, *fp_errfile;
 char  GMMbuf[1000], *token;
 double PDBw;
 double min[3], max[3], Min[3], Max[3];
 int   Ngauss,Ngauss_malloc,Ngauss2, Ngauss_malloc2;
 struct ATOM     HeadAtom;
 struct RESIDUE  HeadResidue;
 struct RADIUS   HeadRadius;
 struct ATOM HeadAtomTar; 
 char   buff[MAX_FILENAME],output_comment[MAX_FILENAME],ChainID,chain;
 char   chainA,chainB,residueA[5],residueB[5],InitType,CalCC_AtmGmm,AtomRadWgt,AtomSelect,RadiusType;
 int    Natom,Natomref,Nrepeat,Srand,Nrepeat_kmeans;
 struct GAUSS3D *GaussArray,*GaussArray2;
 double Mean_Mol[3],CovM_Mol[3][3],iCovM_Mol[3][3],Det_Mol,Min_Mol[3],Max_Mol[3],PCvar_Mol[3],PCaxis_Mol[3][3];
 double logLike_fin,CC,CCingmm,ConvThreEM,VoxelVol,sumDensity;
 double isumDensity, osumDensity;
 double scale_n, scale_d;
 struct VOXEL ivox,ovox,ovox_GMM; 
 char   Atom2VoxType, JustCount;
 float  SigmaAtom2Vox,RadiusUniform; 
 double RESOLUTION;
 double Crr2var, Cww2var, ResolutionAtom, ResolutionGrid; 
 float  ThreValZero, ThreSDZero; 
 int  ReducedSizeScale,MaxVoxelSize,MaxAtomForAll;
 char omapByteOrder; 
 double gorig[3],gref[3],Rmat[3][3],rmsd; 
 int   Nunit_helix; 
 char ABCstr[256];  
 struct MATRIX AtomMember;
 char EMAlgoType,DelZeroWgdf,DelIdenGdf,VarAtomGdfType,VarGridGdfType;
 /** Variables for MarchingCube */
 struct MC_VERTEX MC_Vhead;
 struct MC_FACE   MC_Fhead;
 struct MC_EDGE   MC_Ehead;
 float  MC_SDthre,MC_VOLthre,MC_threshold; 
 int    MC_NVOXthre;
 char   MC_RGBTstr[MAX_FILENAME],MC_SurfaceWireType,MC_ChainID;
 float  MC_RGBT[4];
 float  MC_Vol, MC_Area; 
 int    MC_offset_atomnum;
 float MCGcen[3];
 float  MC_Min[3],MC_Max[3],MC_MinA[3],MC_MaxA[3],MC_MinB[3],MC_MaxB[3]; 
 struct MC_VERTEX MC_VheadA,MC_VheadB;
 struct MC_FACE   MC_FheadA,MC_FheadB;

 /** Variables for Ellipsoid */
 char ellip_color_type;
 double ellip_cover_ratio;
 char ellip_RGBTstr[MAX_FILENAME];
 float ellip_RGBT[4];
 
 sprintf(MODE,"-");

 /*** Initialize to the default parameters ***/
 sprintf(ABCstr,"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
 ipdbfile[0] = iciffile[0] = itarpdbfile[0] = imapfile[0] = '\0';  
 iciffile[0] = assembly_id[0] = '\0';
 PAR.HETtype = 'F';
 ChainID = '-';
 EMAlgoType = 'G'; 
 chainA = chainB = '-';
 sprintf(residueA,"xxx"); sprintf(residueB,"xxx");
 PAR.MODEL = 'S';
 HeadResidue.next = NULL; 
 PAR.Vtype = 'F';
 ovox.grid_width = 4.0;
 Ngauss = 1;
 Nrepeat = 1000;  
 Srand = 0; Nrepeat_kmeans = 10;
 igmmfile[0] = igmmfile2[0] = ovfile[0] = errfile[0] = '\0';    
 ogmmfile[0] = ovoxfile[0]    = omapfile[0] = oimapfile[0] = oobjfile[0] = '\0';
 ogmmfile_init[0] =  iobjfileA[0] = iobjfileB[0] = '\0';
 opdbfile[0]   = owrlfile[0] = oradwgtfile[0] = '\0';
 owrlfile_ellipsoid[0] = '\0'; 
 ompdbfile[0] = impdbfile[0] = '\0'; 
 opcgmmfile[0] = '\0';
 GMMlistfile [0] = '\0';
 PDBlistfile [0] = '\0';
 ERRlistfile [0] = '\0';
 omaperrfile [0] = '\0';
 InitType = 'K';
 ConvThreEM = 0.001;
 Atom2VoxType =  'I';
 SigmaAtom2Vox = 5.0;
 RESOLUTION = 10.0;

 ResolutionGrid = 0.0; 
 ResolutionAtom = 0.0; 

 PAR.ologfile[0] = '\0';
 PAR.OutStdLog = 'T';
 PAR.TabExpMinX = 'T';
 ThreValZero  = -1.0;
 ThreSDZero   = 3.0;
 DelZeroWgdf  = 'T'; 
 DelIdenGdf   = 'T';

 VarAtomGdfType = 'A';
 VarGridGdfType = 'G';

 AtomRadWgt = 'C';
 AtomSelect = 'A';

 HeadAtom.next = NULL;

 CalCC_AtmGmm = 'F';
 MC_Vhead.next = NULL;
 MC_Fhead.next = NULL;
 MC_Ehead.next = NULL;
 MC_VheadA.next = NULL;
 MC_FheadA.next = NULL;
 MC_VheadB.next = NULL;
 MC_FheadB.next = NULL;

 MC_threshold = -1000.0;
 MC_SDthre    = 3.0; 
 MC_NVOXthre  = -1;
 MC_VOLthre   = -1.0;
 
 MC_SurfaceWireType = 'W'; 
 MC_ChainID = 'X'; 
 MC_offset_atomnum = 0; 
 ReducedSizeScale = 1; 
 MaxVoxelSize  = -1;
 MaxAtomForAll = -1; 
 omapByteOrder = 'N';
 sprintf(MC_RGBTstr,"0:1:0:0");
 sprintf(ellip_RGBTstr,"0:1:0:0");
 output_comment[0] = '\0';
 Nunit_helix = 0;
 Crr2var = 1.0/5.0;
 Cww2var = 1.0/12.0;
 RadiusType = 'V';
 ellip_color_type = 'X';
 ellip_cover_ratio  = 0.5;
 RadiusUniform = 1.9;
 JustCount = 'F';

 if (argc<3){
   printf("gmconvert <options>\n");
   printf(" >> for converting atoms/density map into Gaussian Mixture Model <<\n"); 
   printf(" coded by T.Kawabata. LastModified:%s\n",LAST_MOD_DATE);
   printf("<simple usages>\n");
   printf("*[PDB atom] --> [GMM]\n");
   printf("  gmconvert -ipdb [PDBfile] -ogmm [GMMfile] -ng [Number of gdf]\n");
   printf("*[3D density map] --> [GMM]\n");
   printf("  gmconvert -imap [CCP4/MRC map file] -ogmm [GMMfile] -ng [Number of gdf]\n");
   printf("*[PDB atom] --> [3D density map]\n");
   printf("  gmconvert -ipdb [PDBfile] -omap [CCP4 map file] -reso [resolution]\n");
   printf("*[GMM] --> [Wireframe surface in PDB]\n");
   printf("  gmconvert -igmm [GMM] -opdb [PDB] -gw [grid_width]\n");
   printf("*[GMM] --> [3D density map]\n");
   printf("  gmconvert -igmm [GMM] -omap [PDB] -gw [grid_width]\n");
   printf("*[3D density map] --> [Wireframe surface in PDB]\n");
   printf("  gmconvert -imap [CCP4/MRC mapfile] -opdb [PDB] -gw [grid_width]\n");
   printf("*Comparison between [GMM] and [3D density map]\n");
   printf("  gmconvert -igmm [GMMfile] -imap [CCP4/MRC map file]\n");
   printf("*Comparison between [GMM] and another [GMM]\n");
   printf("  gmconvert -igmm [GMMfile] -igmm2 [GMMfile2]\n");
   printf("*Transform [GMM] to target [GMM] by comparison bwn [PDB] and target [PDB]\n");
   printf("  gmconvert -ipdb [PDB] -itpdb [targetPDB] -igmm [GMM] -ogmm [target GMM]\n");
   printf("*Make helical unit structure by comparison bwn [PDB] and target [PDB]\n");
   printf("  gmconvert -ipdb [PDB] -itpdb [targetPDB] -nhelix [Nunit] -opdb [helical unit PDB] -igmm [GMM] -ogmm [helical unit GMM]\n");
   printf("<for details of options>\n");
   printf("  gmconvert -help\n");
   if ((argc==2)&&((strcmp(argv[1],"-h")==0)||(strcmp(argv[1],"-help")==0))) {
     printf(" -M    : MODE 'A2G' atom-to-gmm    'V2G' voxel-to-gmm, 'A2V'atom-to-voxel\n");
     printf("       :      'V2S' voxel-to-shape 'G2S' gmm-to-shape\n");
     printf("       :      'GcmpG':gmm-vs-gmm comparison 'GcmpV' gmm-voxel-comparison[%s]\n",MODE);
     printf("       :      'Atar' :atom-target transform[%s] Ahelix:atom-target helical transform\n",MODE);
     printf("<options for input>\n"); 
     printf(" -ipdb   : Input PDB   file for atomic model [%s]\n",ipdbfile);
     printf(" -icif   : Input mmCIF file for atomic model [%s]\n",iciffile);
     printf(" -imap   : Input CCP4/MRC density map file (*.map/*.mrc) [%s]\n",imapfile);
     printf(" -igmm   : Input file for Gaussian Mixture Model [%s]\n",igmmfile);
     printf(" -igmm2  : Input file for 2nd Gaussian Mixture Model [%s]\n",igmmfile2);
     printf(" -iobjA  : Input file for Object File     (surface model)  (*.obj) [%s]\n",iobjfileA);
     printf(" -iobjB  : Input file for 2nd Object File (surface model)  (*.obj) [%s]\n",iobjfileB);
     printf(" -justcnt: just count input number_of_atom /map_size, and quit.  ('T' or 'F') [%c]\n",JustCount);
     printf("<options for output>\n"); 
     printf(" -omap  : Output CCP4 density map file (*.map) [%s]\n",omapfile);
     printf(" -ogmm  : Output Gaussian File (*.gmm) [%s]\n",ogmmfile);
     printf(" -ovox  : Output Voxel File [%s]\n",ovoxfile);
     printf(" -ogmmpc: Output Gaussian File with PC axis (*.gmm) [%s]\n",opcgmmfile);
     printf("<options for output for surface/wire model>\n"); 
     printf(" -opdb  : Output PDB File    (only wire model)             (*.pdb) [%s]\n",opdbfile);
     printf(" -owrl  : Output VRML File   (both surface and wire model) (*.wrl) [%s]\n",owrlfile);
     printf(" -oobj  : Output Object File (only surface model)          (*.obj) [%s]\n",oobjfile);
     printf(" -oewrl : Output Ellipsoid VRML File   (both surface and wire model) (*.wrl) [%s]\n",owrlfile_ellipsoid);
     printf("<options for input atomic model (-ipdb or -icif)>\n");
     printf(" -hetatm : Read HETATM ('T' or 'F') [%c]\n",PAR.HETtype);
     printf(" -ch     : Chain ID. (or 'auth_asym_id' in mmCIF). [%c]\n",ChainID);
     printf(" -atmsel : Atom selection. 'A'll atom except hydrogen, 'R'esidue-based (only ' CA ' and ' P ') [%c]\n",AtomSelect);
     printf(" -maxatm : maximum allowed number of atoms for '-atmsel A'. If over '-minatmA', then change '-atmsel R'.[%d]\n",MaxAtomForAll);
     printf(" -model  : 'S':read only single model (for NMR). 'M':read multiple models (for biological unit) [%c]\n",PAR.MODEL);
     printf(" -atmrw  : Model for radius and weight. 'A':atom model, 'R':residue model, 'U':uniform raidus/weight,'C':decide from content. [%c]\n",AtomRadWgt);
     printf(" -radtype: radius type for '-atmrw A'. V:van der Waals radius, C:covalent radius  [%c]\n",RadiusType);
     printf(" -raduni : radius for uniform model for '-atmrw U'. [%f]\n",RadiusUniform);
     printf(" -varatm : Variance type for atom (for '-emalg G'). 'A': var = rr2var * Rvdw*Rvdw for each atom, 'R': var = (resoatm/2.0)^2.[%c]\n",VarAtomGdfType);
     printf(" -rr2var : Constant for variance = Const * Rvdw*Rvdw for -emalg G -varatm A. Default is 1/5 [%lf].\n",Crr2var); 
     printf(" -resoatm: Resolution for atom for -emalg G -varatm R.  [%lf].\n",ResolutionAtom); 
     printf(" -orw    : output radius-weight file [%s]\n",oradwgtfile);
     printf(" -ccatm  : Calculation Corr Coeff bwn Atoms and GMMs.(It takes times..) (T or F)[%c]\n",CalCC_AtmGmm);
     printf(" -assembly: assembly_id for mmCIF file (-icif) [%s]\n",assembly_id);
     printf("<options for input map (-imap)>\n");
     printf(" -zth     : if density < [-zth], it is regarded as zero density. [%f]\n",ThreValZero);
     printf(" -zsd     : if density < MEAN + [-zsd]*SD, it is regarded as zero density. [%f]\n",ThreSDZero);
     printf(" ** Negative values for [-zth] and [-zsd] are not meaningful.\n");
     printf(" ** If both [-zth] and [-zsd] are non-negative, [-zth] has a priority.\n");
     printf(" -maxsize : maximum allowed voxel size of each axis(if over, reducing size). [%d]\n",MaxVoxelSize);
     printf(" -redsize : reducing size scale (2,3,4,...)  [%d]\n",ReducedSizeScale);
     printf(" -vargrd  : Variance type for grid (for -emalg G). 'G':var = ww2var * grid_width * grid_width, 'R': var = (resolution/2.0)^2.[%c]\n",VarGridGdfType);
     printf(" -ww2var  : Constant for variance = Const * grid_width*grid_width for emalg G. Default is 1/12 [%lf].\n",Cww2var); 
     printf(" -resogrd : Resolution for grid for -emalg G -vargrd R.  [%lf].\n",ResolutionGrid); 
     printf(" -oimap   : Output CCP4 density map file (*.map) for -imap after modification [%s]\n",oimapfile);
     printf("<options for output 3D density map (-omap)>\n");
     printf(" -gw   : grid width (angstrom) [%f]\n",ovox.grid_width);
     printf(" -a2v  : Type for atom-to-voxel gaussian. 'I':isotropic 'R':Rvdw ('-mcth' should be 0.5) [%c]\n",Atom2VoxType);
     printf(" -reso : Resolution for atom-to-vox (sigma = reso/2) for -a2v I[%lf]\n",RESOLUTION);
     printf(" -sigma: Sigma      for atom-to-vox (sigma = reso/2)  for -a2v I[%f]\n",SigmaAtom2Vox);
     printf("<options for EM algorithm>\n");
     printf(" -ng   : Number of Gaussian Distribution Functions(gdf) for GMM [%d]\n",Ngauss);
     printf(" -emalg: type for EM algorithm. 'P'oint-input EM,'G'mm-input EM, 'O':one-to-one_atom/grid [%c]\n",EMAlgoType); 
     printf(" -I    : Initialization of GMM. 'K'-means, 'R'andom 'O':one-to-one_atom/grid [%c]\n",InitType); 
     printf(" -delzw : Delete Zero-weight gdfs from the GMM. ('T' or 'F') [%c]\n",DelZeroWgdf);
     printf(" -delid : Delete identical gdfs in the GMM. ('T' or 'F') [%c]\n",DelIdenGdf);
     printf(" -nr   : Number of repeat for EM [%d]\n",Nrepeat);
     printf(" -nk   : Number of repeat for K-means multi-start [%d]\n",Nrepeat_kmeans);
     printf(" -cv   : Convergence threshold for EM logLike [%lf]\n",ConvThreEM);
     printf(" -srand: Seed of rand [%d]\n",Srand);
     printf(" -olog : Output logfile [%s]\n",PAR.ologfile);
     printf(" -stdlog: output convergence log as stdout ('T' or 'F') [%c]\n",PAR.OutStdLog);
     printf(" -ogmminit : Output Initial Gaussian File before EM algorithm (*.gmm) [%s]\n",ogmmfile_init);
     printf(" -tex  : Use table for exp(-x) calculation 'T' or 'F' [%c]\n",PAR.TabExpMinX);
     printf("<options for transforming models>\n"); 
     printf(" -itpdb  : Input target PDB file for transforming GMM with symmetric configurations. [%s]\n",itarpdbfile);
     printf("  ( transforming [-igmm] GMM file into [-ogmm] by the transformation of [-ipdb] to the pose of [-itpdb].) \n");
     printf(" -nhelix : number of generating subunits in helical symmetry [%d]\n",Nunit_helix);
     printf("<options for making solid model(Marching Cube algorithm)>\n");
     printf(" -mcSW   : model type. 'S'urface, 'W'ireframe [%c]\n",MC_SurfaceWireType); 
     printf(" -mcsd   : SD value          for threshold density [%lf]\n",MC_SDthre); 
     printf(" -mcth   : raw density value for threshold density [%lf]\n",MC_threshold); 
     printf(" -mcnv   : number for voxel  for threshold density [%d]\n",MC_NVOXthre); 
     printf(" -mcvo   : volume (A^3)      for threshold density [%lf]\n",MC_VOLthre); 
     printf(" -mcRGBT : RGBT string (red:green:blue:transparency)  [%s]\n",MC_RGBTstr); 
     printf(" -mcch     ChainID for output [%c]\n",MC_ChainID); 
     printf(" -mcofan : offset atom number [%d]\n",MC_offset_atomnum); 
     printf("<options for ellipsoidal model>\n"); 
     printf(" -oewrl : Output Ellipsoid VRML File   (both surface and wire model) (*.wrl) [%s]\n",owrlfile_ellipsoid);
     printf(" -elcov : cover ratio for ellipsoid magnitude (0..1). (larger value -> large ellipsoid) [%lf]\n",ellip_cover_ratio);
     printf(" -elcol : color scheme. 'W':by gdf weight. 'N':by gdf number, otherwise:by '-elRGBT' option [%c]\n",ellip_color_type);
     printf(" -elRGBT: RGBT string (red:green:blue:transparency)  [%s]\n",ellip_RGBTstr); 
     printf("<options for input/output GMM memberships for each atom>\n"); 
     printf(" -ompdb  : Output PDB File with membership values (*.pdb) [%s]\n",ompdbfile);
     printf(" -impdb  : Input  PDB File with membership values (*.pdb) [%s]\n",impdbfile);
  }
  exit(1);
 }

 /****** READ ARGUMENTS ***********/
 PAR.COMMAND[0] = '\0';
 PAR.START_TIME_SEC = Get_Time_in_Second_by_Double();
 sprintf(PAR.START_DATE,"%s",Get_Date_String());

 for (k=0;k<argc;++k){ 
   if (k>0) strcat(PAR.COMMAND," ");
   strcat(PAR.COMMAND,argv[k]);
 }

 printf("#COMMAND %s\n",PAR.COMMAND);
 k = 1;
 while (k<argc){
    if (argv[k][0]=='-'){
      L = strlen(argv[k]);
           if ((L==5)&&(argv[k][1]=='i')&&(argv[k][2]=='g')&&(argv[k][3]=='m')&&(argv[k][4]=='m')) {++k; sprintf(igmmfile,"%s",argv[k]); }
      else if ((L==6)&&(argv[k][1]=='c')&&(argv[k][2]=='c')&&(argv[k][3]=='a')&&(argv[k][4]=='t')&&(argv[k][5]=='m')) 
          {++k; CalCC_AtmGmm= argv[k][0]; }
      else if ((L==6)&&(argv[k][1]=='i')&&(argv[k][2]=='g')&&(argv[k][3]=='m')&&(argv[k][4]=='m')&&(argv[k][5]=='2')) 
          {++k; sprintf(igmmfile2,"%s",argv[k]); }
      else if ((L==5)&&(argv[k][1]=='i')&&(argv[k][2]=='p')&&(argv[k][3]=='d')&&(argv[k][4]=='b')) {++k; sprintf(ipdbfile,"%s",argv[k]); }
      else if ((L==5)&&(argv[k][1]=='i')&&(argv[k][2]=='c')&&(argv[k][3]=='i')&&(argv[k][4]=='f')) {++k; sprintf(iciffile,"%s",argv[k]); }
      else if ((L==6)&&(argv[k][1]=='i')&&(argv[k][2]=='t')&&(argv[k][3]=='p')&&(argv[k][4]=='d')&&(argv[k][5]=='b')) {++k; sprintf(itarpdbfile,"%s",argv[k]); }
      else if ((L==5)&&(argv[k][1]=='i')&&(argv[k][2]=='m')&&(argv[k][3]=='a')&&(argv[k][4]=='p')) {++k; sprintf(imapfile,"%s",argv[k]); }
      else if ((L==5)&&(argv[k][1]=='o')&&(argv[k][2]=='g')&&(argv[k][3]=='m')&&(argv[k][4]=='m')) {++k; sprintf(ogmmfile,"%s",argv[k]); }
      else if ((L==9)&&(argv[k][1]=='o')&&(argv[k][2]=='g')&&(argv[k][3]=='m')&&(argv[k][4]=='m')&&(argv[k][5]=='i')&&(argv[k][6]=='n')&&(argv[k][7]=='i')&&(argv[k][8]=='t')){
       ++k; sprintf(ogmmfile_init,"%s",argv[k]); }
      else if ((L==7)&&(argv[k][1]=='o')&&(argv[k][2]=='g')&&(argv[k][3]=='m')&&(argv[k][4]=='m')&&(argv[k][5]=='p')&&(argv[k][6]=='c')) 
           {++k; sprintf(opcgmmfile,"%s",argv[k]); }
      else if ((L==5)&&(argv[k][1]=='o')&&(argv[k][2]=='w')&&(argv[k][3]=='r')&&(argv[k][4]=='l')) {++k; sprintf(owrlfile,"%s",argv[k]); }
      else if ((L==6)&&(argv[k][1]=='o')&&(argv[k][2]=='e')&&(argv[k][3]=='w')&&(argv[k][4]=='r')&&(argv[k][5]=='l')) 
         {++k; sprintf(owrlfile_ellipsoid,"%s",argv[k]); }
      else if ((L==3)&&(argv[k][1]=='c')&&(argv[k][2]=='v')) {++k; ConvThreEM = (double)atof(argv[k]); }
      else if ((L==3)&&(argv[k][1]=='c')&&(argv[k][2]=='h')) {++k; ChainID = argv[k][0];}
      else if ((L==3)&&(argv[k][1]=='c')&&(argv[k][2]=='A')) {++k; chainA = argv[k][0];}
      else if ((L==3)&&(argv[k][1]=='c')&&(argv[k][2]=='B')) {++k; chainB = argv[k][0];}
      else if ((L==3)&&(argv[k][1]=='r')&&(argv[k][2]=='A')) {++k; sprintf(residueA,"%s",argv[k]) ;}
      else if ((L==3)&&(argv[k][1]=='r')&&(argv[k][2]=='B')) {++k; sprintf(residueB,"%s",argv[k]) ;}
      else if ((L==3)&&(argv[k][1]=='n')&&(argv[k][2]=='g')) {++k; Ngauss = atoi(argv[k]);}
      else if ((L==3)&&(argv[k][1]=='o')&&(argv[k][2]=='v')) {++k; sprintf(ovfile,"%s",argv[k]);}
      else if ((L==4)&&(argv[k][1]=='e')&&(argv[k][2]=='r')&&(argv[k][3]=='r')) {++k; sprintf(errfile,"%s",argv[k]);}
      else if ((L==5)&&(argv[k][1]=='o')&&(argv[k][2]=='v')&&(argv[k][3]=='o')&&(argv[k][4]=='x')) {++k; sprintf(ovoxfile,"%s",argv[k]);}
      else if ((L==5)&&(argv[k][1]=='g')&&(argv[k][2]=='m')&&(argv[k][3]=='m')&&(argv[k][4]=='l')) {++k; sprintf(GMMlistfile,"%s",argv[k]);}
      else if ((L==5)&&(argv[k][1]=='p')&&(argv[k][2]=='d')&&(argv[k][3]=='b')&&(argv[k][4]=='l')) {++k; sprintf(PDBlistfile,"%s",argv[k]);}
      else if ((L==5)&&(argv[k][1]=='e')&&(argv[k][2]=='r')&&(argv[k][3]=='r')&&(argv[k][4]=='l')) {++k; sprintf(ERRlistfile,"%s",argv[k]);}
      else if ((L==5)&&(argv[k][1]=='e')&&(argv[k][2]=='m')&&(argv[k][3]=='a')&&(argv[k][4]=='p')) {++k; sprintf(omaperrfile,"%s",argv[k]);}
      else if ((L==6)&&(argv[k][1]=='m')&&(argv[k][2]=='o')&&(argv[k][3]=='d')&&(argv[k][4]=='e')&&(argv[k][5]=='l')) {++k;PAR.MODEL = argv[k][0];}
      else if ((L==5)&&(argv[k][1]=='o')&&(argv[k][2]=='m')&&(argv[k][3]=='a')&&(argv[k][4]=='p')) {++k; sprintf(omapfile,"%s",argv[k]);}
      else if ((L==6)&&(argv[k][1]=='o')&&(argv[k][2]=='i')&&(argv[k][3]=='m')&&(argv[k][4]=='a')&&(argv[k][5]=='p')) {++k; sprintf(oimapfile,"%s",argv[k]);}
      else if ((L==5)&&(argv[k][1]=='o')&&(argv[k][2]=='o')&&(argv[k][3]=='b')&&(argv[k][4]=='j')) {++k; sprintf(oobjfile,"%s",argv[k]);}
      else if ((L==6)&&(argv[k][1]=='i')&&(argv[k][2]=='o')&&(argv[k][3]=='b')&&(argv[k][4]=='j')&&(argv[k][5]=='A')) {++k; sprintf(iobjfileA,"%s",argv[k]);}
      else if ((L==6)&&(argv[k][1]=='i')&&(argv[k][2]=='o')&&(argv[k][3]=='b')&&(argv[k][4]=='j')&&(argv[k][5]=='B')) {++k; sprintf(iobjfileB,"%s",argv[k]);}
      else if ((L==5)&&(argv[k][1]=='o')&&(argv[k][2]=='p')&&(argv[k][3]=='d')&&(argv[k][4]=='b')) {++k; sprintf(opdbfile,"%s",argv[k]);}
      else if ((L==6)&&(argv[k][1]=='o')&&(argv[k][2]=='m')&&(argv[k][3]=='p')&&(argv[k][4]=='d')&&(argv[k][5]=='b')) {++k; sprintf(ompdbfile,"%s",argv[k]);}
      else if ((L==6)&&(argv[k][1]=='i')&&(argv[k][2]=='m')&&(argv[k][3]=='p')&&(argv[k][4]=='d')&&(argv[k][5]=='b')) {++k; sprintf(impdbfile,"%s",argv[k]);}
      else if ((L==5)&&(argv[k][1]=='m')&&(argv[k][2]=='c')&&(argv[k][3]=='t')&&(argv[k][4]=='h')) {++k; MC_threshold = atof(argv[k]);}
      else if ((L==5)&&(argv[k][1]=='m')&&(argv[k][2]=='c')&&(argv[k][3]=='s')&&(argv[k][4]=='d')) {++k; MC_SDthre    = atof(argv[k]);}
      else if ((L==5)&&(argv[k][1]=='m')&&(argv[k][2]=='c')&&(argv[k][3]=='n')&&(argv[k][4]=='v')) {++k; MC_NVOXthre  = atoi(argv[k]);}
      else if ((L==5)&&(argv[k][1]=='m')&&(argv[k][2]=='c')&&(argv[k][3]=='v')&&(argv[k][4]=='o')) {++k; MC_VOLthre   = atof(argv[k]);}
      else if ((L==5)&&(argv[k][1]=='m')&&(argv[k][2]=='c')&&(argv[k][3]=='c')&&(argv[k][4]=='h')) {++k; MC_ChainID   = argv[k][0];}
      else if ((L==5)&&(argv[k][1]=='m')&&(argv[k][2]=='c')&&(argv[k][3]=='S')&&(argv[k][4]=='W')) {++k; MC_SurfaceWireType = argv[k][0];}
      else if ((L==6)&&(argv[k][1]=='e')&&(argv[k][2]=='m')&&(argv[k][3]=='a')&&(argv[k][4]=='l')&&(argv[k][5]=='g')){
        ++k; EMAlgoType = argv[k][0];
      }
      else if ((L==8)&&(argv[k][1]=='r')&&(argv[k][2]=='a')&&(argv[k][3]=='d')&&(argv[k][4]=='t')&&(argv[k][5]=='y')&&(argv[k][6]=='p')&&(argv[k][7]=='e')){
        ++k; RadiusType = argv[k][0];
      }
      else if ((L==7)&&(argv[k][1]=='r')&&(argv[k][2]=='a')&&(argv[k][3]=='d')&&(argv[k][4]=='u')&&(argv[k][5]=='n')&&(argv[k][6]=='i')){
        ++k; RadiusUniform = atof(argv[k]);
      }
      else if ((L==6)&&(argv[k][1]=='a')&&(argv[k][2]=='t')&&(argv[k][3]=='m')&&(argv[k][4]=='r')&&(argv[k][5]=='w')){
        ++k; AtomRadWgt = argv[k][0];
      }  
      else if ((L==7)&&(argv[k][1]=='a')&&(argv[k][2]=='t')&&(argv[k][3]=='m')&&(argv[k][4]=='s')&&(argv[k][5]=='e')&&(argv[k][6]=='l')){
        ++k; AtomSelect = argv[k][0];
      }  
      else if ((L==6)&&(argv[k][1]=='d')&&(argv[k][2]=='e')&&(argv[k][3]=='l')&&(argv[k][4]=='z')&&(argv[k][5]=='w')){
        ++k; DelZeroWgdf = argv[k][0];
      }
      else if ((L==6)&&(argv[k][1]=='d')&&(argv[k][2]=='e')&&(argv[k][3]=='l')&&(argv[k][4]=='i')&&(argv[k][5]=='d')){
        ++k; DelIdenGdf = argv[k][0];
      }
      else if ((L==7)&&(argv[k][1]=='m')&&(argv[k][2]=='c')&&(argv[k][3]=='R')&&(argv[k][4]=='G')&&(argv[k][5]=='B')&&(argv[k][6]=='T')){
        ++k; sprintf(MC_RGBTstr,"%s",argv[k]);
      }
      else if ((L==7)&&(argv[k][1]=='m')&&(argv[k][2]=='c')&&(argv[k][3]=='o')&&(argv[k][4]=='f')&&(argv[k][5]=='a')&&(argv[k][6]=='n')){
        ++k; MC_offset_atomnum = atoi(argv[k]);
      }
      else if ((L==3)&&(argv[k][1]=='g')&&(argv[k][2]=='w')) {++k; ovox.grid_width = atof(argv[k]);}
      else if ((L==3)&&(argv[k][1]=='n')&&(argv[k][2]=='r')) {++k; Nrepeat     = atoi(argv[k]);}
      else if ((L==3)&&(argv[k][1]=='n')&&(argv[k][2]=='k')) {++k; Nrepeat_kmeans = atoi(argv[k]);}
      else if ((L==6)&&(argv[k][1]=='s')&&(argv[k][2]=='r')&&(argv[k][3]=='a')&&(argv[k][4]=='n')&&(argv[k][5]=='d')) {++k; Srand = atoi(argv[k]);}
      else if ((L==4)&&(argv[k][1]=='a')&&(argv[k][2]=='2')&&(argv[k][3]=='v')) 
         {++k; Atom2VoxType = argv[k][0];}
      else if ((L==4)&&(argv[k][1]=='o')&&(argv[k][2]=='r')&&(argv[k][3]=='w')) 
         {++k; sprintf(oradwgtfile,"%s",argv[k]);}
      else if ((L==6)&&(argv[k][1]=='s')&&(argv[k][2]=='i')&&(argv[k][3]=='g')&&(argv[k][4]=='m')&&(argv[k][5]=='a')) 
         {++k; SigmaAtom2Vox = atof(argv[k]);}
      else if ((L==5)&&(argv[k][1]=='r')&&(argv[k][2]=='e')&&(argv[k][3]=='s')&&(argv[k][4]=='o')) 
         {++k; RESOLUTION = atof(argv[k]); SigmaAtom2Vox = RESOLUTION/2.0; }
      else if ((L==5)&&(argv[k][1]=='o')&&(argv[k][2]=='l')&&(argv[k][3]=='o')&&(argv[k][4]=='g')) {++k; sprintf(PAR.ologfile,"%s",argv[k]);}
      else if ((L==7)&&(argv[k][1]=='s')&&(argv[k][2]=='t')&&(argv[k][3]=='d')&&(argv[k][4]=='l')&&(argv[k][5]=='o')&&(argv[k][6]=='g'))
      {++k; PAR.OutStdLog = argv[k][0];}
      else if ((L==7)&&(argv[k][1]=='h')&&(argv[k][2]=='e')&&(argv[k][3]=='t')&&(argv[k][4]=='a')&&(argv[k][5]=='t')&&(argv[k][6]=='m')) {++k; PAR.HETtype = argv[k][0]; }
      else if ((L==2)&&(argv[k][1]=='I')) {++k; InitType = argv[k][0]; }
      else if ((L==2)&&(argv[k][1]=='M')) {++k; sprintf(MODE,"%s",argv[k]); }
      else if ((L==4)&&(argv[k][1]=='z')&&(argv[k][2]=='t')&&(argv[k][3]=='h')) {++k; ThreValZero = atof(argv[k]);  }
      else if ((L==4)&&(argv[k][1]=='z')&&(argv[k][2]=='s')&&(argv[k][3]=='d')) {++k; ThreSDZero = atof(argv[k]);  }
      else if ((L==4)&&(argv[k][1]=='s')&&(argv[k][2]=='u')&&(argv[k][3]=='b')) 
           {++k; sprintf(PAR.SubType,"%s",argv[k]);}
      else if ((L==4)&&(argv[k][1]=='t')&&(argv[k][2]=='e')&&(argv[k][3]=='x')) 
           {++k; PAR.TabExpMinX = argv[k][0];}
      else if ((L==8)&&(argv[k][1]=='r')&&(argv[k][2]=='e')&&(argv[k][3]=='d')&&(argv[k][4]=='s')&&(argv[k][5]=='i')&&(argv[k][6]=='z')&&(argv[k][7]=='e')) 
           {++k; ReducedSizeScale = atoi(argv[k]);}
      else if ((L==8)&&(argv[k][1]=='m')&&(argv[k][2]=='a')&&(argv[k][3]=='x')&&(argv[k][4]=='s')&&(argv[k][5]=='i')&&(argv[k][6]=='z')&&(argv[k][7]=='e')) 
           {++k; MaxVoxelSize = atoi(argv[k]);}
      else if ((L==7)&&(argv[k][1]=='m')&&(argv[k][2]=='a')&&(argv[k][3]=='x')&&(argv[k][4]=='a')&&(argv[k][5]=='t')&&(argv[k][6]=='m')) 
           {++k; MaxAtomForAll = atoi(argv[k]);}
      else if ((L==7)&&(argv[k][1]=='n')&&(argv[k][2]=='h')&&(argv[k][3]=='e')&&(argv[k][4]=='l')&&(argv[k][5]=='i')&&(argv[k][6]=='x')){
        ++k; Nunit_helix = atoi(argv[k]);
      }
      else if ((L==7)&&(argv[k][1]=='r')&&(argv[k][2]=='r')&&(argv[k][3]=='2')&&(argv[k][4]=='v')&&(argv[k][5]=='a')&&(argv[k][6]=='r')){
        ++k; Crr2var = atof(argv[k]);
      }
      else if ((L==7)&&(argv[k][1]=='w')&&(argv[k][2]=='w')&&(argv[k][3]=='2')&&(argv[k][4]=='v')&&(argv[k][5]=='a')&&(argv[k][6]=='r')){
        ++k; Cww2var = atof(argv[k]);
      }
      else if ((L==8)&&(argv[k][1]=='r')&&(argv[k][2]=='e')&&(argv[k][3]=='s')&&(argv[k][4]=='o')&&(argv[k][5]=='a')&&(argv[k][6]=='t')&&(argv[k][7]=='m')){
        ++k; ResolutionAtom = atof(argv[k]);
      }
      else if ((L==8)&&(argv[k][1]=='r')&&(argv[k][2]=='e')&&(argv[k][3]=='s')&&(argv[k][4]=='o')&&(argv[k][5]=='g')&&(argv[k][6]=='r')&&(argv[k][7]=='d')){
        ++k; ResolutionGrid = atof(argv[k]);
      }
      else if ((L==7)&&(argv[k][1]=='v')&&(argv[k][2]=='a')&&(argv[k][3]=='r')&&(argv[k][4]=='g')&&(argv[k][5]=='r')&&(argv[k][6]=='d')){
        ++k; VarGridGdfType = argv[k][0];
      }
      else if ((L==7)&&(argv[k][1]=='v')&&(argv[k][2]=='a')&&(argv[k][3]=='r')&&(argv[k][4]=='a')&&(argv[k][5]=='t')&&(argv[k][6]=='m')){
        ++k; VarAtomGdfType = argv[k][0];
      }
      else if ((L==9)&&(argv[k][1]=='a')&&(argv[k][2]=='s')&&(argv[k][3]=='s')&&(argv[k][4]=='e')&&(argv[k][5]=='m')&&(argv[k][6]=='b')&&(argv[k][7]=='l')&&(argv[k][8]=='y')){
        ++k; sprintf(assembly_id,"%s",argv[k]);
      }
      else if ((L==6)&&(argv[k][1]=='e')&&(argv[k][2]=='l')&&(argv[k][3]=='c')&&(argv[k][4]=='o')&&(argv[k][5]=='v')){
        ++k; ellip_cover_ratio = atof(argv[k]);
      }

      else if ((L==6)&&(argv[k][1]=='e')&&(argv[k][2]=='l')&&(argv[k][3]=='c')&&(argv[k][4]=='o')&&(argv[k][5]=='l')){
        ++k; ellip_color_type = argv[k][0];
      }
      else if ((L==7)&&(argv[k][1]=='e')&&(argv[k][2]=='l')&&(argv[k][3]=='R')&&(argv[k][4]=='G')&&(argv[k][5]=='B')&&(argv[k][6]=='T')){
        ++k; sprintf(ellip_RGBTstr,"%s",argv[k]);
      }
      else if ((L==8)&&(argv[k][1]=='j')&&(argv[k][2]=='u')&&(argv[k][3]=='s')&&(argv[k][4]=='t')&&(argv[k][5]=='c')&&(argv[k][6]=='n')&&(argv[k][7]=='t')){
        ++k; JustCount = argv[k][0];
      }
      else { printf("#ERROR:Can't understand option '%s'.\n",argv[k]); exit(1);}

    } /* if '-' */

    ++k;

  } /* while k */

  /*****************************/
  /*** Guess MODE from input ***/
  /*****************************/
  
  if (strcmp(MODE,"-")==0){
        if ((ipdbfile[0] != '\0')&&(omapfile[0]!='\0')) sprintf(MODE,"A2V");
   else if ((ipdbfile[0] != '\0') && (itarpdbfile[0] != '\0') && (Nunit_helix>1)) sprintf(MODE,"Ahel");
   else if ((ipdbfile[0] != '\0') && (itarpdbfile[0] != '\0')) sprintf(MODE,"Atar");
   else if (ipdbfile[0] != '\0')  sprintf(MODE,"A2G");
   else if (iciffile[0] != '\0')  sprintf(MODE,"A2G");
   else if (impdbfile[0] != '\0') sprintf(MODE,"AM2G");
   else if ((igmmfile[0] != '\0') && (imapfile[0] != '\0')) sprintf(MODE,"GcmpV");
   else if ((igmmfile[0] != '\0') && (igmmfile2[0] != '\0')) sprintf(MODE,"GcmpG");
   else if ((igmmfile[0] != '\0') && (opdbfile[0]  != '\0')) sprintf(MODE,"G2S");
   else if ((igmmfile[0] != '\0') && (owrlfile[0] != '\0')) sprintf(MODE,"G2S");
   else if ((igmmfile[0] != '\0') && (omapfile[0] != '\0')) sprintf(MODE,"G2S");
   else if ((igmmfile[0] != '\0') && (oobjfile[0] != '\0')) sprintf(MODE,"G2S");
   else if ((igmmfile[0] != '\0') && (owrlfile_ellipsoid[0] != '\0')) sprintf(MODE,"G2E");
   else if ((imapfile[0] != '\0') && (opdbfile[0]  != '\0')) sprintf(MODE,"V2S");
   else if ((imapfile[0] != '\0') && (owrlfile[0] != '\0')) sprintf(MODE,"V2S");
   else if ((imapfile[0] != '\0') && (oobjfile[0] != '\0')) sprintf(MODE,"V2S");
   else if ((iobjfileA[0] != '\0') && (iobjfileB[0] != '\0')) sprintf(MODE,"ScmpS");
   else if (imapfile[0] != '\0') sprintf(MODE,"V2G");
   else{
     printf("#ERROR:required to input  PDB file (-ip) or CCP4 map file (-im).\n");
     exit(1);   
   }
  }

  printf("#MODE:%s\n",MODE);


  if ((strcmp(MODE,"A2G")==0)||(strcmp(MODE,"A2V")==0)||(strcmp(MODE,"AM2G")==0)){
    Set_Default_Radius(&HeadRadius,RadiusType); 
  }

  /** Set Logfile **/
  if (PAR.ologfile[0]!='\0'){ 
   PAR.fp_olog = fopen(PAR.ologfile,"w");
   fprintf(PAR.fp_olog,"#COMMAND %s\n",PAR.COMMAND);
   fprintf(PAR.fp_olog,"#DATE    %s\n",Get_Date_String());
  }


  /** Set Table for exp(-x) **/
  if (PAR.TabExpMinX == 'T'){ 
    Make_Exp_minus_X_Table(&TABEXP,20.0, 0.01);
  } 

   RGBTstr_to_floatRGBT(MC_RGBTstr,MC_RGBT);
   RGBTstr_to_floatRGBT(ellip_RGBTstr,ellip_RGBT);
  /*********************************/
  /*** MODE 'A2G': PDB -> Gauss  ***/
  /*********************************/
 
  if (strcmp(MODE,"A2G")==0){ 
    /** Managing File name **/ 
    Find_Filename_Core(corename,ipdbfile);
    if (ChainID!='-') 
     { sprintf(buff,"%s",corename); 
       sprintf(corename,"%s%c",buff,ChainID); }
  
    if (ogmmfile[0]=='\0') sprintf(ogmmfile,"%s_%d.gmm",corename,Ngauss);
 
    if (ipdbfile[0] != '\0'){ 
      Read_PDB_File(ipdbfile,&HeadAtom,AtomSelect,ChainID);
    }
    else if (iciffile[0] != '\0'){ 
      Read_mmCIF_File(iciffile,&HeadAtom,AtomSelect,ChainID,assembly_id,MaxAtomForAll);
    }

    Natom = Number_Of_Atom(&HeadAtom);
    printf("NATOM %d\n",Natom); 
    printf("NCHAIN %d\n",Number_Of_Chain(&HeadAtom));
    if (JustCount=='T'){ exit(1);}

    if ((ipdbfile[0]!='\0')&&(AtomSelect=='A') && (MaxAtomForAll > 0) && (Natom > MaxAtomForAll)){
      Free_ATOMs(&HeadAtom);
      Read_PDB_File(ipdbfile,&HeadAtom,'R',ChainID);
      printf("#Natom %d -->(by change residue model)-->",Natom); 
      Natom = Number_Of_Atom(&HeadAtom);
      printf(" %d\n",Natom); 
    } 


    if ((AtomRadWgt != 'A') && (AtomRadWgt != 'R') && (AtomRadWgt != 'U')){ 
       AtomRadWgt = Judge_Atom_or_Residue_Model(&HeadAtom);
    }

    if (AtomRadWgt=='A'){
      Assign_Radius(&HeadAtom,&HeadRadius);
      Assign_Weight(&HeadAtom);
    }
    else if (AtomRadWgt=='R'){
      Assign_Radius_by_Residue(&HeadAtom);
      Assign_Weight_by_Residue(&HeadAtom);
    }
    else if (AtomRadWgt=='U'){
      Assign_Uniform_Radius(&HeadAtom,RadiusUniform);
      Assign_Uniform_Weight(&HeadAtom);
    }


    if (oradwgtfile[0] != '\0'){ write_radius_weight_of_atoms(oradwgtfile,&HeadAtom);}

   if (opdbfile[0] != '\0'){
     /* sprintf(output_comment,"",1); */
     Write_PDB_File(opdbfile,&HeadAtom,'w','-',"",PAR.COMMAND);
   }



    Set_Variance_of_Atoms(&HeadAtom, VarAtomGdfType, Crr2var, ResolutionAtom);

    if (Natom==0){ 
      printf("#ERROR:PDBfile \"%s\" does not contain any atom.\n",ipdbfile);
      exit(1);
    }
 /*
    if (Natom < Ngauss){
      printf("#WARNING:Natom %d is less than Ngauss %d. Ngauss is set to %d.\n",Ngauss,Natom,Natom);
      Ngauss = Natom;
    }
 */
   /* Make_Residue_List(&HeadAtom,&HeadResidue); */
  
   /* Cal Mean and CovarMat for Protein Atoms */
    Cal_Mean_and_CovarMat_from_Atoms(&HeadAtom,Mean_Mol,CovM_Mol,Min_Mol,Max_Mol);
    Cal_Inverse_Matrix3D_by_Cramer_Rule(iCovM_Mol, CovM_Mol,&Det_Mol);
    printf("#Det_Mol %lf\n",Det_Mol);
   
    GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss); 

    if (InitType=='K'){
     K_Means_Clustering_For_Atoms_Multiple_Start(Ngauss,GaussArray,&HeadAtom,Crr2var,&logLike_fin,Nrepeat_kmeans); 
     for (g=0;g<Ngauss;++g){
       printf("#g %d W %lf M %lf %lf %lf\n",g,GaussArray[g].Weight, GaussArray[g].M[0], GaussArray[g].M[1], GaussArray[g].M[2]);
     }
    }
    else if (InitType=='O'){
      free(GaussArray);
      Ngauss = Number_Of_Atom(&HeadAtom);
      GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss); 
      GaussMix_by_One_Atom_One_Gdf_Assignment(Ngauss,GaussArray,&HeadAtom);
    }
    else{
      srand(Srand);
      Initialize_Gauss_From_Randomly_Chosen_Atoms(&HeadAtom,Ngauss,GaussArray);
    }

    if (ogmmfile_init[0] != '\0'){
      Write_Gaussian3D_File(ogmmfile_init,'w',Ngauss,GaussArray,Chain_Of_Atoms(&HeadAtom),""); 
    }

    if (EMAlgoType=='G'){  
      EM_optimize_GaussMix_For_Gaussian_Atoms(Nrepeat,Ngauss,GaussArray,&HeadAtom,&logLike_fin,ConvThreEM);
    }
    else if (EMAlgoType=='O'){  
      free(GaussArray);
      Ngauss = Number_Of_Atom(&HeadAtom);
      GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss); 
      GaussMix_by_One_Atom_One_Gdf_Assignment(Ngauss,GaussArray,&HeadAtom);
      Write_Gaussian3D_File("1vs1.gmm",'w',Ngauss,GaussArray,Chain_Of_Atoms(&HeadAtom),""); 
    } 
    else{
      EM_optimize_GaussMix_For_Atoms(Nrepeat,Ngauss,GaussArray,&HeadAtom,&logLike_fin,ConvThreEM);
    }

    if (DelZeroWgdf == 'T') { Delete_ZeroWeight_gdfs_of_GMM(&Ngauss,GaussArray);}
    if (DelIdenGdf  == 'T') { Delete_Identical_gdfs_of_GMM(&Ngauss,GaussArray);}

    /* printf("#Ngauss %d logLike_fin %f\n",Ngauss,logLike_fin); */
  
    /*** output gaussians ***/ 
    if (ogmmfile[0]!='\0'){ Write_Gaussian3D_File(ogmmfile,'w',Ngauss,GaussArray,Chain_Of_Atoms(&HeadAtom),""); } 

    if (opcgmmfile[0]!='\0'){
      for (g=0;g<Ngauss;++g){ Cal_EigenVectors_For_Matrix3D(GaussArray[g].CovM,GaussArray[g].PCvar,GaussArray[g].PCaxis); }
      Write_Gaussian3D_File_with_PCaxis(opcgmmfile,'w',Ngauss,GaussArray,Chain_Of_Atoms(&HeadAtom),"");
    }

    CC = 0.0;
    if (CalCC_AtmGmm=='T'){
      CC = CorrCoeff_Bwn_Atoms_and_GMM(&HeadAtom,Ngauss, GaussArray);
      printf("Sigma %lf Reso %lf Ngauss %d CC %f\n",SigmaAtom2Vox,2.0*SigmaAtom2Vox,Ngauss,CC);
    }

    printf("NGAUSS %d CorrCoeff %lf logLike %e Sigma %lf\n",Ngauss,CC,logLike_fin,SigmaAtom2Vox);
 
   /*** output voxels ***/
   if ((ovoxfile[0]!= '\0')||(omapfile[0]!='\0')||(owrlfile[0]!='\0')){
     Malloc_Voxel_From_Gaussians(&ovox,Ngauss,GaussArray);
     Set_Voxel_Value_By_Gaussians(&ovox,Ngauss,GaussArray);
     if (ovoxfile[0]!='\0'){Write_Voxel_File(&ovox,ovoxfile);} 
     if (omapfile[0]!='\0'){Write_MapCCP4_File(&ovox,omapfile);} 
   }

   if (ompdbfile[0] != '\0'){
    Malloc_MATRIX(&AtomMember,Natom,Ngauss);
    Cal_Memberships_for_Atoms(&AtomMember,&HeadAtom,Ngauss,GaussArray);
    Write_PDB_File_With_GMM_and_Memberships(ompdbfile,&HeadAtom,Ngauss,GaussArray,&AtomMember,"",PAR.COMMAND);
    Estimate_M_and_CovM_from_Memberships(Ngauss,GaussArray,&HeadAtom,&AtomMember);
    /*
    Write_Gaussian3D_File("reest.gmm",'w',Ngauss,GaussArray,Chain_Of_Atoms(&HeadAtom),"");
    */
    Free_MATRIX(&AtomMember);
   }




 }

  /***************************************/
  /*** MODE 'V2G':  CCP4_MAP -> Gauss  ***/
  /***************************************/
  
  if (strcmp(MODE,"V2G")==0){
   Read_Density_File(imapfile,&ivox,'T',-1.0);
   /* printf("#grid_width %lf N %d %d %d\n",ivox.grid_width,ivox.N[0],ivox.N[1],ivox.N[2]);  */
   Find_Filename_Core(corename,imapfile);
   
   if ((ivox.N[0]<=2)||(ivox.N[1]<=2)||(ivox.N[2]<=2)){
     printf("#ERROR(vox_length): density map ('%s') is too small vox length(%d %d %d).\n",
       imapfile,ivox.N[0],ivox.N[1],ivox.N[2]);
     exit(1); 
   } 

   if (ogmmfile[0]=='\0') sprintf(ogmmfile,"%s_%d.gmm",corename,Ngauss);

   if (ThreValZero>=0.0){
     Set_Less_Than_Threshold_Voxels_to_Zero(&ivox,ThreValZero);  
     printf("THRESHOLD_VALUE %f\n",ThreValZero);
   }
   else if (ThreSDZero>0.0){ 
     Voxel_Statistics(&ivox,&(ivox.min),&(ivox.max),&(ivox.ave),&(ivox.sd));
     Set_Less_Than_Threshold_Voxels_to_Zero(&ivox,ivox.ave + ThreSDZero * ivox.sd);  
     printf("DENSITY_AVERAGE %f\n",ivox.ave);
     printf("DENSITY_SD      %f\n",ivox.sd);  
     printf("THRESHOLD_VALUE %f\n",ivox.ave + ThreSDZero * ivox.sd);  
 
   }

   if (JustCount=='T'){ exit(1);}
   
   if (MaxVoxelSize>0){
     VoxelVol = (double)(ivox.N[0]*ivox.N[1]*ivox.N[2]);
     if (VoxelVol>(double)(MaxVoxelSize*MaxVoxelSize*MaxVoxelSize)){
        ReducedSizeScale = (int)ceil(pow(VoxelVol,1.0/3.0)/(double)MaxVoxelSize); 
     }
     printf("#ReducedSizeByMaxVoxelSize ivox size %d %d %d MaxVoxelSize %d --> ReducedSizeScale %d\n",ivox.N[0], ivox.N[1], ivox.N[2], MaxVoxelSize, ReducedSizeScale); 
   }

   if (ReducedSizeScale>=2){
     for (k=0;k<3;++k){ovox.N[k] = ivox.N[k]/ReducedSizeScale;}
     Malloc_Voxel(&ovox,ovox.N[0],ovox.N[1],ovox.N[2]);
     Make_Reduced_Size_Voxel(&ovox,&ivox,ReducedSizeScale);
     Free_Voxel(&ivox);
     for (k=0;k<3;++k){ivox.N[k] = ovox.N[k];}
     Malloc_Voxel(&ivox,ivox.N[0],ivox.N[1],ivox.N[2]);
     Copy_Voxel(&ivox,&ovox);
     Free_Voxel(&ovox);
   }
       
   // now read file with list of GMMs to rescale ivox
   if (GMMlistfile[0]!='\0'){
      // prepare vox files
      // allocate space for output ovox
      Malloc_Voxel(&ovox,ivox.N[0],ivox.N[1],ivox.N[2]); 
      // set parameters equal to input map
      ovox.grid_width = ivox.grid_width;
      ovox.OrigPos[0] = ivox.OrigPos[0];
      ovox.OrigPos[1] = ivox.OrigPos[1];
      ovox.OrigPos[2] = ivox.OrigPos[2];
      // allocate space for output map of i-th GMM component
      Malloc_Voxel(&ovox_GMM,ivox.N[0],ivox.N[1],ivox.N[2]);      
      // set parameters equal to input map
      ovox_GMM.grid_width = ivox.grid_width;
      ovox_GMM.OrigPos[0] = ivox.OrigPos[0];
      ovox_GMM.OrigPos[1] = ivox.OrigPos[1];
      ovox_GMM.OrigPos[2] = ivox.OrigPos[2];     
      // read  GMMlistfile line by line
      ptr_file = fopen(GMMlistfile,"r");
      while (fgets(GMMbuf,1000, ptr_file)!=NULL){
       token = strtok(GMMbuf, " ");
       sprintf(igmmfile2, "%s", token); 
       token = strtok(NULL, " \t");
       x = atoi(token);
       printf("Reading GMM %s Component %d\n",igmmfile2,x); 
       // read number of Gaussians
       Ngauss2 = Number_Of_Gaussian3D_in_File(igmmfile2);
       // allocate array
       GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss2); 
       // read file
       chain = Read_Gaussian3D_File(igmmfile2,Ngauss2,&Ngauss2,GaussArray);
       // put the value of GMM into ovox
       Set_Voxel_Value_By_Gaussians(&ovox,Ngauss2,GaussArray);
       // put the value of GMM into ovox_GMM
       Set_Voxel_Value_By_Gaussian(&ovox_GMM,&GaussArray[x]);
       // rescale ivox
       for (x=0;x< ivox.N[0];++x){
        for (y=0;y< ivox.N[1];++y){
         for (z=0;z< ivox.N[2];++z){
         // rescale value
          if( ivox.dat[x][y][z] == 0.0 || ovox.dat[x][y][z]  == 0.0){
            ivox.dat[x][y][z] = 0.0;
          } else {         
            ivox.dat[x][y][z] = ivox.dat[x][y][z] * ovox_GMM.dat[x][y][z] / ovox.dat[x][y][z];
          }
         }
        }
       }   
       // deallocate array
       free(GaussArray);
      }
      fclose(ptr_file);
      // deallocate voxels
      Free_Voxel(&ovox); 
      Free_Voxel(&ovox_GMM);
   }
 
   // printout the input map after thresolding and rescaling
   if (oimapfile[0]!='\0') {Write_MapCCP4_File(&ivox,oimapfile);}
   
   Cal_Mean_and_CovarMat_from_Voxel(&ivox,Mean_Mol,CovM_Mol,Min_Mol,Max_Mol);
   Cal_Inverse_Matrix3D_by_Cramer_Rule(iCovM_Mol, CovM_Mol,&Det_Mol);
   Cal_EigenVectors_For_Matrix3D(CovM_Mol,PCvar_Mol,PCaxis_Mol);

   for (i=0;i<3;++i){printf("#%d Min %lf Max %lf\n",i,Min_Mol[i],Max_Mol[i]);}
   printf("MOL_M      %15.10lf %15.10lf %15.10lf\n", Mean_Mol[0], Mean_Mol[1], Mean_Mol[2]);
   printf("MOL_CovM_X %15.10lf %15.10lf %15.10lf\n", CovM_Mol[0][0], CovM_Mol[0][1], CovM_Mol[0][2]);
   printf("MOL_CovM_Y %15.10lf %15.10lf %15.10lf\n", CovM_Mol[1][0], CovM_Mol[1][1], CovM_Mol[1][2]);
   printf("MOL_CovM_Z %15.10lf %15.10lf %15.10lf\n", CovM_Mol[2][0], CovM_Mol[2][1], CovM_Mol[2][2]);
   printf("MOL_PCvar  %15.10lf %15.10lf %15.10lf\n", PCvar_Mol[0], PCvar_Mol[1], PCvar_Mol[2]);

   if (Ngauss==0){
     printf("#ERROR(-ng 0):Cannot estimate GMM with Ngauss=0.\n");
     exit(1);
   }

   GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss); 
   srand(Srand);

   Initialize_Gauss_From_M_and_CovM(Ngauss,GaussArray,Mean_Mol,CovM_Mol);
   if (InitType=='K'){
    K_Means_Clustering_for_3D_Map_Multiple_Start(Ngauss,GaussArray,&ivox,&logLike_fin,Nrepeat_kmeans,VarGridGdfType,Cww2var,ResolutionGrid);
   }
   else if (InitType=='O'){
    free(GaussArray);
    cal_sumDensity_and_Nvox_posi_dens(&ivox,&sumDensity,&Ngauss);
    GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss); 
    GaussMix_by_One_Voxel_One_Gdf_Assignment(&Ngauss,GaussArray,&ivox,VarGridGdfType,Cww2var,ResolutionGrid);
  }
  else {
    Initialize_Gauss_From_M_and_CovM(Ngauss,GaussArray,Mean_Mol,CovM_Mol);
    /*
    Initialize_Gauss_From_Randomly_Chosen_Voxels(&ivox,Ngauss,GaussArray);
    */
  } 

  printf("#DelZeroWgdf %c DelIdenGdf %c\n",DelZeroWgdf,DelIdenGdf); 
  if (DelZeroWgdf == 'T'){ Delete_ZeroWeight_gdfs_of_GMM(&Ngauss,GaussArray);}
  if (DelIdenGdf == 'T') { Delete_Identical_gdfs_of_GMM(&Ngauss,GaussArray);}

  CCingmm = Corr_Coeff_Bwn_VoxelGMM_and_GMM(&ivox,Ngauss,GaussArray,VarGridGdfType,Cww2var,ResolutionGrid);
  printf("#CCingmm_initial %lf\n",CCingmm);
  if (ogmmfile_init[0] != '\0'){
    Write_Gaussian3D_File(ogmmfile_init,'w',Ngauss,GaussArray,'I',""); 
  }

   emok = 0;
 
   if (EMAlgoType=='s'){
     emok = EM_optimize_GaussMix_for_3D_Map_SmallMemory(
           Nrepeat,Ngauss,GaussArray,&ivox,&logLike_fin,ConvThreEM);
   }
   else if (EMAlgoType=='G'){ 
     emok = EM_optimize_GaussMix_for_Gaussian_3D_Map(
           Nrepeat,Ngauss,GaussArray,&ivox,VarGridGdfType,Cww2var,ResolutionGrid,&logLike_fin,ConvThreEM);
   }
   else if (EMAlgoType=='O'){ 
     free(GaussArray);
     cal_sumDensity_and_Nvox_posi_dens(&ivox,&sumDensity,&Ngauss);
     GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss); 
     emok = GaussMix_by_One_Voxel_One_Gdf_Assignment(&Ngauss,GaussArray,&ivox,VarGridGdfType,Cww2var,ResolutionGrid);
   } 
   else { 
     emok = EM_optimize_GaussMix_for_3D_Map(
           Nrepeat,Ngauss,GaussArray,&ivox,&logLike_fin,ConvThreEM);
   } 

   if (emok==0){
     printf("#ERROR(EM_3Dmap): density map ('%s') cannot be changed into GMM due to non-normal parameters.\n",imapfile);
     exit(1); 
   }

  printf("#DelZeroWgdf %c DelIdenGdf %c\n",DelZeroWgdf,DelIdenGdf); 
  if (DelZeroWgdf == 'T'){ Delete_ZeroWeight_gdfs_of_GMM(&Ngauss,GaussArray);}
  if (DelIdenGdf == 'T') { Delete_Identical_gdfs_of_GMM(&Ngauss,GaussArray);}

  /*** output gaussians ***/ 
  /* logLike_fin = Log_Likelihood_Of_GaussMix_For_3Dmap(Ngauss,GaussArray,&ivox); */
  CC      = Corr_Coeff_Bwn_Voxel_and_GMM(&ivox,Ngauss,GaussArray);
  CCingmm = Corr_Coeff_Bwn_VoxelGMM_and_GMM(&ivox,Ngauss,GaussArray,VarGridGdfType,Cww2var,ResolutionGrid);
  sprintf(output_comment,"Corr.Coeff. %lf",CC); 
  printf("NGAUSS %d CorrCoeff %lf CorrCoeffInGMM %lf logLike %e\n",Ngauss,CC,CCingmm,logLike_fin);

  // allocate space for output ovox
  Malloc_Voxel(&ovox,ivox.N[0],ivox.N[1],ivox.N[2]);
  // set parameters equal to input map
  ovox.grid_width = ivox.grid_width;
  ovox.OrigPos[0] = ivox.OrigPos[0];
  ovox.OrigPos[1] = ivox.OrigPos[1];
  ovox.OrigPos[2] = ivox.OrigPos[2];
  // put the GMM into ovox
  Set_Voxel_Value_By_Gaussians(&ovox,Ngauss,GaussArray);

  // calculate scaling factor 
  scale_n = 0.0; scale_d = 0.0;
  // cycle on x - y - z 
  for (x=0;x<ivox.N[0];++x){
    for (y=0;y<ivox.N[1];++y){
      for (z=0;z<ivox.N[2];++z){
       scale_n += ivox.dat[x][y][z] * ovox.dat[x][y][z];
       scale_d += ovox.dat[x][y][z] * ovox.dat[x][y][z];
      }
    }
  }
  
  // multiply GMM weights by scaling factor 
  for (g=0;g<Ngauss;++g) GaussArray[g].Weight *= scale_n / scale_d;
    
  // printout Gaussian file
  if (ogmmfile[0]!='\0') Write_Gaussian3D_File(ogmmfile,'w',Ngauss,GaussArray,'I',output_comment);
  
 }

  
 /*********************************************/
 /*** MODE 'GcmpG': Gauss-Gauss comparison  ***/
 /*********************************************/
 if (strcmp(MODE,"GcmpG")==0){

  Ngauss_malloc = Number_Of_Gaussian3D_in_File(igmmfile);
  GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss_malloc); 
  chain = Read_Gaussian3D_File(igmmfile,Ngauss_malloc,&Ngauss,GaussArray);

  Ngauss_malloc2 = Number_Of_Gaussian3D_in_File(igmmfile2);
  GaussArray2 = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss_malloc2); 
  chain = Read_Gaussian3D_File(igmmfile2,Ngauss_malloc2,&Ngauss2,GaussArray2);

  printf("IGAUSSFILE  %s\n", igmmfile);
  printf("IGAUSSFILE2 %s\n", igmmfile2);
  printf("CC %lf\n",Corr_Coeff_Bwn_Two_GAUSS3D_Arrays(Ngauss,GaussArray,Ngauss2,GaussArray2));
 }

 
 /*****************************************/
 /*** MODE 'GcmpV': Gauss-Voxel comparison  ***/
 /*****************************************/
 if (strcmp(MODE,"GcmpV")==0){
   Ngauss_malloc = Number_Of_Gaussian3D_in_File(igmmfile);
   GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss_malloc); 
   chain = Read_Gaussian3D_File(igmmfile,Ngauss_malloc,&Ngauss,GaussArray);

   Read_Density_File(imapfile,&ivox,'T',-1.0);
   /*
   printf("#grid_width %lf N %d %d %d\n",ivox.grid_width,ivox.N[0],ivox.N[1],ivox.N[2]); 
   */

   for (i=0;i<3;++i){
     ovox.N[i] = ivox.N[i];
     ovox.OrigPos[i] = ivox.OrigPos[i];
   }
   ovox.grid_width = ivox.grid_width; 
   Malloc_Voxel(&ovox,ovox.N[0],ovox.N[1],ovox.N[2]);
   Set_GaussMix_Density_To_3Dmap(Ngauss,GaussArray,&ovox);

   // write ovox to file 
   if (omapfile[0]!='\0') Write_MapCCP4_File(&ovox,omapfile); 

   // calculate integral of ivox and ovox
   isumDensity = 0.0;
   osumDensity = 0.0;
   // cycle on x
   for (x=0;x<ivox.N[0];++x){
     // set bin size
     double dx = ivox.grid_width;
     if(x==0 || x==ivox.N[0]-1) dx /= 2.0;
     // cycle on y
     for (y=0;y<ivox.N[1];++y){
       // set bin size
       double dy = ivox.grid_width;
       if(y==0 || y==ivox.N[1]-1) dy /= 2.0;
       // cycle on z
       for (z=0;z<ivox.N[2];++z){
        // set bin size
        double dz = ivox.grid_width;
        if(z==0 || z==ivox.N[2]-1) dz /= 2.0;
        isumDensity += ivox.dat[x][y][z] * dx * dy * dz;
        osumDensity += ovox.dat[x][y][z] * dx * dy * dz;
       }
     }
   }

   printf("IGAUSSFILE %s\n", igmmfile);
   printf("IMAPFILE   %s\n", imapfile);
   printf("CC %lf\n",Corr_Coeff_Bwn_Two_Voxels(&ivox,&ovox));
   printf("INTEGRAL IVOX (grid) %lf\n", isumDensity);
   printf("INTEGRAL OVOX (grid) %lf\n", osumDensity);

   // calculate integral of GMM as sum of weights
   osumDensity = 0.0;
   for (g=0;g<Ngauss;++g) osumDensity += GaussArray[g].Weight;
   printf("INTEGRAL OVOX (GMM)  %lf\n", osumDensity);

   // calculate overlaps
   if(ovfile[0] != '\0'){
     printf("Writing overlap file with GMM error %s\n", ovfile);
     fp_ovfile = fopen(ovfile,"w");
     // allocate grid for individual components
     Malloc_Voxel(&ovox_GMM,ivox.N[0],ivox.N[1],ivox.N[2]);
     // set parameters equal to input map
     ovox_GMM.grid_width = ivox.grid_width;
     ovox_GMM.OrigPos[0] = ivox.OrigPos[0];
     ovox_GMM.OrigPos[1] = ivox.OrigPos[1];
     ovox_GMM.OrigPos[2] = ivox.OrigPos[2];
     // cycle on components
     for(g=0; g<Ngauss; ++g){
      // put the value of GMM into ovox_GMM
      Set_Voxel_Value_By_Gaussian(&ovox_GMM,&GaussArray[g]);
      // calculate overlap
      isumDensity = 0.0;
      osumDensity = 0.0;
      // cycle on x
      for (x=0;x<ivox.N[0];++x){
        // set bin size
        double dx = ivox.grid_width;
        if(x==0 || x==ivox.N[0]-1) dx /= 2.0;
        // cycle on y
        for (y=0;y<ivox.N[1];++y){
          // set bin size
          double dy = ivox.grid_width;
          if(y==0 || y==ivox.N[1]-1) dy /= 2.0;
          // cycle on z
          for (z=0;z<ivox.N[2];++z){
           // set bin size
           double dz = ivox.grid_width;
           if(z==0 || z==ivox.N[2]-1) dz /= 2.0;
           isumDensity += ovox_GMM.dat[x][y][z] * ivox.dat[x][y][z] * dx * dy * dz;
           osumDensity += ovox_GMM.dat[x][y][z] * ovox.dat[x][y][z] * dx * dy * dz;
          }
        }
      }
      // print to file
      fprintf(fp_ovfile,"%e %e %e\n", isumDensity, osumDensity, fabs(isumDensity-osumDensity));
     }
     fclose(fp_ovfile);
     // free stuff
     Free_Voxel(&ovox_GMM);
   }

   // calculate overlap errors from variance map
   if(errfile[0] != '\0'){
     printf("Writing overlap errors to file %s\n", errfile);
     fp_errfile = fopen(errfile,"w");
     // allocate grid for individual components
     Malloc_Voxel(&ovox_GMM,ivox.N[0],ivox.N[1],ivox.N[2]);
     // set parameters equal to input map
     ovox_GMM.grid_width = ivox.grid_width;
     ovox_GMM.OrigPos[0] = ivox.OrigPos[0];
     ovox_GMM.OrigPos[1] = ivox.OrigPos[1];
     ovox_GMM.OrigPos[2] = ivox.OrigPos[2];
     // cycle on components
     for(g=0; g<Ngauss; ++g){
      // put the value of GMM into ovox_GMM
      Set_Voxel_Value_By_Gaussian(&ovox_GMM,&GaussArray[g]);
      // calculate error with propagation from overlap formula
      isumDensity = 0.0;
      // cycle on x
      for (x=0;x<ivox.N[0];++x){
        // set bin size
        double dx = ivox.grid_width;
        if(x==0 || x==ivox.N[0]-1) dx /= 2.0;
        // cycle on y
        for (y=0;y<ivox.N[1];++y){
          // set bin size
          double dy = ivox.grid_width;
          if(y==0 || y==ivox.N[1]-1) dy /= 2.0;
          // cycle on z
          for (z=0;z<ivox.N[2];++z){
           // set bin size
           double dz = ivox.grid_width;
           if(z==0 || z==ivox.N[2]-1) dz /= 2.0;
           double der = ovox_GMM.dat[x][y][z] * dx * dy * dz;
           isumDensity += der * der * ivox.dat[x][y][z];
          }
        }
      }
      // print to file
      fprintf(fp_errfile,"%e\n", sqrt(isumDensity));
     }
     fclose(fp_errfile);
     // free stuff
     Free_Voxel(&ovox_GMM);
   }

}

 /**********************************/
 /*** MODE 'G2S': GMM -> Shape  ***/
 /*********************************/

 if (strcmp(MODE,"G2S")==0){ 
   Ngauss_malloc = Number_Of_Gaussian3D_in_File(igmmfile);
   GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss_malloc); 
   chain = Read_Gaussian3D_File(igmmfile,Ngauss_malloc,&Ngauss,GaussArray);
   Malloc_Voxel_From_Gaussians(&ovox,Ngauss,GaussArray,5.0);
   /*
   Malloc_Voxel_From_Gaussians(&ovox,Ngauss,GaussArray,8.0);
   */  
   Set_Voxel_Value_By_Gaussians(&ovox,Ngauss,GaussArray);
   if (omapfile[0] != '\0') Write_MapCCP4_File(&ovox,omapfile,omapByteOrder);
   
   // read  error file
   if (ERRlistfile[0]!='\0'){
     // open error file
     ptr_file = fopen(ERRlistfile,"r");
     printf("Reading error file %s\n", ERRlistfile);
     // cycle on file
     g=0;
     while (fgets(GMMbuf,1000, ptr_file)!=NULL){
       // convert to double
       PDBw = atof(GMMbuf);
       // printout
       //printf("Error of component %d : %f\n", g, PDBw);
       // rescale gaussian
       GaussArray[g].Weight *= PDBw;
       // increment counter
       ++g;
     }
     // close file
     fclose(ptr_file);

     // prepare new map
     Malloc_Voxel(&ovox_GMM,ovox.N[0],ovox.N[1],ovox.N[2]);      
     // set parameters equal to omap
     ovox_GMM.grid_width = ovox.grid_width;
     ovox_GMM.OrigPos[0] = ovox.OrigPos[0];
     ovox_GMM.OrigPos[1] = ovox.OrigPos[1];
     ovox_GMM.OrigPos[2] = ovox.OrigPos[2];
     // add Gaussians to map
     Set_Voxel_Value_By_Gaussians(&ovox_GMM,Ngauss,GaussArray);
     // divide by original map
     for (x=0;x< ovox.N[0];++x){
        for (y=0;y< ovox.N[1];++y){
         for (z=0;z< ovox.N[2];++z){
          // rescale value
          if( ovox.dat[x][y][z] != 0.0 ){
           ovox_GMM.dat[x][y][z] /= ovox.dat[x][y][z];
          } else {
           ovox_GMM.dat[x][y][z] = 0.0;
          }
         }
        }
     }   
     // write to file
     if (omaperrfile[0] != '\0') Write_MapCCP4_File(&ovox_GMM,omaperrfile,omapByteOrder);
        
   }

   Voxel_Statistics(&ovox,&(ovox.min),&(ovox.max),&(ovox.ave),&(ovox.sd));
   if (MC_threshold<-100.0){
     if (MC_VOLthre>0){
       MC_NVOXthre = (int)(MC_VOLthre/(ovox.grid_width*ovox.grid_width*ovox.grid_width));
       MC_threshold = Find_Threshold_Density_For_NumVoxel(&ovox,MC_NVOXthre);
       printf("#MC_VOLthre %.1f -> MC_NVOXthre %d -> MC_threshold %e\n",MC_VOLthre,MC_NVOXthre,MC_threshold);
     }
     else if (MC_NVOXthre>0){ MC_threshold = Find_Threshold_Density_For_NumVoxel(&ovox,MC_NVOXthre); }
     else if (MC_SDthre>0.0){MC_threshold = ovox.ave + ovox.sd * MC_SDthre;}
   }
   Marching_Cube_Tetrahedral(&ovox,&MC_Vhead,&MC_Fhead,&MC_Ehead,MC_threshold,'C');
  
   MC_Vol = Volume_Inside_Faces(ovox.grid_width,&MC_Vhead,&MC_Fhead);
   MC_Area = Area_Faces(ovox.grid_width,&MC_Fhead);
   sprintf(buff,"Surface_Model: Volume %.2f (A^3) Area %.2f (A^2)",MC_Vol,MC_Area);
   printf("#%s\n",buff); strcat(output_comment,buff);
 
   if (oobjfile[0] != '\0'){ 
     Output_MC_Faces_Object(oobjfile,'w',&MC_Vhead,&MC_Fhead,&ovox,output_comment);
   }
   if (owrlfile[0] != '\0'){
     if (MC_SurfaceWireType=='W'){
       Output_MC_Edges_VRML(owrlfile,'w',&MC_Vhead,&MC_Ehead,&ovox,MC_RGBT,output_comment);
     }
     else{
       Output_MC_Faces_VRML(owrlfile,'w',&MC_Vhead,&MC_Fhead,&ovox,MC_RGBT,output_comment);
     }
   }
   if (opdbfile[0] != '\0'){
     n = Output_MC_Edges_PDB(opdbfile,'w',&MC_Vhead,&MC_Ehead,&ovox,1.0,output_comment,MC_ChainID,MC_offset_atomnum);
     Write_Gaussian3D_Center_Points(opdbfile,'a',Ngauss,GaussArray,MC_ChainID,MC_offset_atomnum+n,"");
   }




}

 /*************************************/
 /*** MODE 'G2E': GMM -> Ellipsoids ***/
 /*************************************/
 if (strcmp(MODE,"G2E")==0){ 
   Ngauss_malloc = Number_Of_Gaussian3D_in_File(igmmfile);
   GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss_malloc); 
   chain = Read_Gaussian3D_File(igmmfile,Ngauss_malloc,&Ngauss,GaussArray);
   write__GAUSS3Ds_in_ellipsoid_VRML(owrlfile_ellipsoid,Ngauss, GaussArray,ellip_cover_ratio,ellip_color_type,ellip_RGBT);
 }


 /**********************************/
 /*** MODE 'V2S': Voxel -> Shape  ***/
 /*********************************/

 if (strcmp(MODE,"V2S")==0){ 
   if (imapfile[0]!='\0') {Read_Density_File(imapfile,&ivox,'T',-1.0);}
   /*
   printf("#grid_width %lf N %d %d %d\n",ivox.grid_width,ivox.N[0],ivox.N[1],ivox.N[2]); 
   */

   if (MaxVoxelSize>0){
     if ((ivox.N[0]>MaxVoxelSize)|| (ivox.N[1]>MaxVoxelSize)|| (ivox.N[2]>MaxVoxelSize)){
            if ((ivox.N[0]>=ivox.N[1])&&(ivox.N[0]>=ivox.N[2])){ReducedSizeScale = (int)ceil((double)ivox.N[0]/(double)MaxVoxelSize);} 
       else if ((ivox.N[1]>=ivox.N[0])&&(ivox.N[1]>=ivox.N[2])){ReducedSizeScale = (int)ceil((double)ivox.N[1]/(double)MaxVoxelSize);} 
       else if ((ivox.N[2]>=ivox.N[0])&&(ivox.N[2]>=ivox.N[1])){ReducedSizeScale = (int)ceil((double)ivox.N[2]/(double)MaxVoxelSize);} 
     }
     printf("#ReducedSizeByMaxVoxelSize ivox size %d %d %d MaxVoxelSize %d --> ReducedSizeScale %d\n",ivox.N[0], ivox.N[1], ivox.N[2], MaxVoxelSize, ReducedSizeScale); 
   }

   if (ReducedSizeScale>=2){
     for (k=0;k<3;++k){ovox.N[k] = ivox.N[k]/ReducedSizeScale;}
     Malloc_Voxel(&ovox,ovox.N[0],ovox.N[1],ovox.N[2]);
     Make_Reduced_Size_Voxel(&ovox,&ivox,ReducedSizeScale);
     Free_Voxel(&ivox);
     for (k=0;k<3;++k){ivox.N[k] = ovox.N[k];}
     Malloc_Voxel(&ivox,ivox.N[0],ivox.N[1],ivox.N[2]);
     Copy_Voxel(&ivox,&ovox);
     Free_Voxel(&ovox);
     if (omapfile[0]!='\0') {Write_MapCCP4_File(&ivox,omapfile);}
   }

   if (omapfile[0] != '\0') Write_MapCCP4_File(&ivox,omapfile,omapByteOrder);

   Voxel_Statistics(&ivox,&(ivox.min),&(ivox.max),&(ivox.ave),&(ivox.sd));
   if (MC_threshold < -100.0){
     if (MC_VOLthre>0.0){
       MC_NVOXthre = (int)(MC_VOLthre/(ivox.grid_width*ivox.grid_width*ivox.grid_width));
       MC_threshold = Find_Threshold_Density_For_NumVoxel(&ovox,MC_NVOXthre);
       printf("#MC_VOLthre %.1f -> MC_NVOXthre %d -> MC_threshold %e\n",MC_VOLthre,MC_NVOXthre,MC_threshold);
     }
     else if (MC_NVOXthre>0){ MC_threshold = Find_Threshold_Density_For_NumVoxel(&ivox,MC_NVOXthre); }
     else if (MC_SDthre>0.0){ MC_threshold = ivox.ave + ivox.sd * MC_SDthre;}
   }

   Marching_Cube_Tetrahedral(&ivox,&MC_Vhead,&MC_Fhead,&MC_Ehead,MC_threshold,'C');

   MC_Vol = Volume_Inside_Faces(ivox.grid_width,&MC_Vhead,&MC_Fhead);
   MC_Area = Area_Faces(ivox.grid_width,&MC_Fhead);
   Cal_Gcenter_Of_Vertexes(MCGcen,&MC_Vhead);
   sprintf(buff,"Surface_Model: Volume %.2f (A^3) Area %.2f (A^2)",MC_Vol,MC_Area);
   printf("#%s\n",buff); strcat(output_comment,buff);
   printf("#MCGcenter %f %f %f\n", ivox.grid_width*MCGcen[0], ivox.grid_width*MCGcen[1], ivox.grid_width*MCGcen[2]);

   if (oobjfile[0] != '\0'){ 
     Output_MC_Faces_Object(oobjfile,'w',&MC_Vhead,&MC_Fhead,&ivox,output_comment);
   }

   if (owrlfile[0] != '\0'){
     if (MC_SurfaceWireType=='W'){
       Output_MC_Edges_VRML(owrlfile,'w',&MC_Vhead,&MC_Ehead,&ivox,MC_RGBT,"");
     }
     else{
       Output_MC_Faces_VRML(owrlfile,'w',&MC_Vhead,&MC_Fhead,&ivox,MC_RGBT,"");
     }
   }
   if (opdbfile[0] != '\0'){
     Output_MC_Edges_PDB(opdbfile,'w',&MC_Vhead,&MC_Ehead,&ivox,1.0,"",MC_ChainID,MC_offset_atomnum);
   }
}


 /**********************************/
 /*** MODE 'A2V': PDB -> Voxel  ***/
 /*********************************/
 if (strcmp(MODE,"A2V")==0){
  //  read file with list of pdbs 
   if (PDBlistfile[0]!='\0'){
      // reset Min and Max
      for(i=0; i<3; ++i) Min[i]=1.0e+9;
      for(i=0; i<3; ++i) Max[i]=-1.0e+9;     
      // read  PDBlistfile line by line
      ptr_file = fopen(PDBlistfile,"r");
      while (fgets(GMMbuf,1000, ptr_file)!=NULL){
       token = strtok(GMMbuf, " ");
       sprintf(ipdbfile, "%s", token); 
       token = strtok(NULL, " \t");
       PDBw = atof(token);
       printf("Reading PDB %s Weight %lf\n",ipdbfile,PDBw);
       Find_Filename_Core(corename,ipdbfile);
       if (ChainID!='-') 
        { sprintf(buff,"%s",corename); sprintf(corename,"%s%c",buff,ChainID); }
       Read_PDB_File(ipdbfile,&HeadAtom,AtomSelect,ChainID);
       Assign_Radius(&HeadAtom,&HeadRadius);
       Natom = Number_Of_Atom(&HeadAtom);
       printf("Natom %d\n",Natom); 
  
       if (Natom==0)
       { printf("#ERROR:PDBfile \"%s\" does not contain any atom.\n",ipdbfile);
        exit(1); }
     
       // find minimum and maximum in the current pdb
       Get_Min_Max_From_Atoms(&HeadAtom,&min,&max);
       // check for global min and max
       for(i=0; i<3; ++i) if(min[i]<Min[i]) Min[i] = min[i];
       for(i=0; i<3; ++i) if(max[i]>Max[i]) Max[i] = max[i];
          
       // Free Memory
       Free_ATOMs(&HeadAtom);
      }
      // printout
      printf("[0]Min %f Max %f [1]Min %f Max %f [2]Min %f Max %f\n",
      Min[0], Max[0], Min[1], Max[1], Min[2], Max[2]);     
      // close list file
      fclose(ptr_file);
      // allocate stuff
      Malloc_Voxel_From_Min_Max(&ovox,Min,Max,3.0*SigmaAtom2Vox);
      
      // reopen list
      ptr_file = fopen(PDBlistfile,"r");
      // counter
      i=0;
      while (fgets(GMMbuf,1000, ptr_file)!=NULL){
       token = strtok(GMMbuf, " ");
       sprintf(ipdbfile, "%s", token); 
       token = strtok(NULL, " \t");
       PDBw = atof(token);
       Find_Filename_Core(corename,ipdbfile);
       if (ChainID!='-') 
        { sprintf(buff,"%s",corename); sprintf(corename,"%s%c",buff,ChainID); }
       Read_PDB_File(ipdbfile,&HeadAtom,AtomSelect,ChainID);
       Assign_Radius(&HeadAtom,&HeadRadius);
       Natom = Number_Of_Atom(&HeadAtom);
       
       // add stuff to ovox
       if(i==0) Set_Voxel_Value_By_Atoms(&ovox,&HeadAtom,Atom2VoxType,SigmaAtom2Vox,true);
       else     Set_Voxel_Value_By_Atoms(&ovox,&HeadAtom,Atom2VoxType,SigmaAtom2Vox,false);
       
       // increase counter
       ++i;
       
       // Free Memory
       Free_ATOMs(&HeadAtom);
      }
      
      // divide map by the number of frames
      for (x=0;x< ovox.N[0];++x){
       for (y=0;y< ovox.N[1];++y){
        for (z=0;z< ovox.N[2];++z){
           ovox.dat[x][y][z] /= (double) i;
        }
       }
      }
      
      // print map
      if (ovoxfile[0]!='\0') Write_Voxel_File(&ovox,ovoxfile); 
      if (omapfile[0]!='\0') Write_MapCCP4_File(&ovox,omapfile);
    
   } else {
    // traditional single PDB mode      
    Find_Filename_Core(corename,ipdbfile);
    if (ChainID!='-') 
    { sprintf(buff,"%s",corename); sprintf(corename,"%s%c",buff,ChainID); }
    Read_PDB_File(ipdbfile,&HeadAtom,AtomSelect,ChainID);
    Assign_Radius(&HeadAtom,&HeadRadius);
    Natom = Number_Of_Atom(&HeadAtom);
    printf("Natom %d\n",Natom); 
  
    if (Natom==0)
    { printf("#ERROR:PDBfile \"%s\" does not contain any atom.\n",ipdbfile);
     exit(1); }
  
    Malloc_Voxel_From_Atoms(&ovox,&HeadAtom,3.0*SigmaAtom2Vox);
    Set_Voxel_Value_By_Atoms(&ovox,&HeadAtom,Atom2VoxType,SigmaAtom2Vox,true);

    if (ovoxfile[0]!='\0') Write_Voxel_File(&ovox,ovoxfile); 
    if (omapfile[0]!='\0') Write_MapCCP4_File(&ovox,omapfile);

  }
 }

 /****************************************************/
 /*** MODE 'AM2G': Atoms with memberships  -> GMM  ***/
 /****************************************************/
 if (strcmp(MODE,"AM2G")==0){

   Read_PDB_File(impdbfile,&HeadAtom,AtomSelect,ChainID);
   Natom = Number_Of_Atom(&HeadAtom);
   printf("Natom %d\n",Natom); 
   if (Natom==0){ 
     printf("#ERROR:PDBfile \"%s\" does not contain any atom.\n",ipdbfile);
     exit(1); }
  
   Ngauss     = Number_Of_Gaussian3D_in_File(impdbfile);
   GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss); 
   printf("#Ngauss %d\n",Ngauss); 
   Malloc_MATRIX(&AtomMember,Natom,Ngauss);
   Read_PDB_File_With_GMM_and_Memberships(impdbfile,Ngauss,GaussArray,&AtomMember);
   Estimate_M_and_CovM_from_Memberships(Ngauss,GaussArray,&HeadAtom,&AtomMember);
   if (ogmmfile[0]!='\0'){ Write_Gaussian3D_File(ogmmfile,'w',Ngauss,GaussArray,Chain_Of_Atoms(&HeadAtom),""); } 
   if (opcgmmfile[0]!='\0'){
     for (g=0;g<Ngauss;++g){ Cal_EigenVectors_For_Matrix3D(GaussArray[g].CovM,GaussArray[g].PCvar,GaussArray[g].PCaxis); }
     Write_Gaussian3D_File_with_PCaxis(opcgmmfile,'w',Ngauss,GaussArray,Chain_Of_Atoms(&HeadAtom),"");
   }
   Free_MATRIX(&AtomMember);
 }



 /*********************************************/
 /*** MODE 'Atar': Atom reference transform ***/
 /*********************************************/
 if (strcmp(MODE,"Atar")==0){ 
   Find_Filename_Core(corename,ipdbfile);
   if (ChainID!='-') 
    { sprintf(buff,"%s",corename); sprintf(corename,"%s%c",buff,ChainID); }
  
   Read_PDB_File(ipdbfile,&HeadAtom,AtomSelect,ChainID);
   /* Assign_Radius(&HeadAtom,&HeadRadius); */
   Natom = Number_Of_Atom(&HeadAtom);
   printf("#Natom %d\n",Natom); 

   Read_PDB_File(itarpdbfile,&HeadAtomTar,AtomSelect,ChainID);
   /* Assign_Radius(&HeadAtomTar,&HeadRadius); */
   Natomref = Number_Of_Atom(&HeadAtomTar);
   printf("#Natomref %d\n",Natomref); 

  /*
   rmsd = Calculate_CRMS_Bwn_Two_ATOMs(&HeadAtom,&HeadAtomTar,gorig,gref,Rmat);
   */
   rmsd = Calculate_CRMS_Bwn_Two_ATOMs_Using_CA_With_Same_Rnum(&HeadAtom,&HeadAtomTar,gorig,gref,Rmat);
 
   printf("#gorig %lf %lf %lf gref %lf %lf %lf\n",gorig[0],gorig[1],gorig[2],gref[0],gref[1],gref[2]);
   printf("#rmsd %lf\n",rmsd);
  
   if (igmmfile[0] != '\0'){ 
     Ngauss_malloc = Number_Of_Gaussian3D_in_File(igmmfile);
     GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss_malloc); 
     chain = Read_Gaussian3D_File(igmmfile,Ngauss_malloc,&Ngauss,GaussArray);
     Transform_Gaussians_By_Rmat_Gorig_Gnew(Ngauss, GaussArray, Rmat,gorig, gref);
     if (ogmmfile[0]!='\0'){ 
       Write_Gaussian3D_File(ogmmfile,'w',Ngauss,GaussArray,Chain_Of_Atoms(&HeadAtomTar),"");
     } 
   }
 }
 /****************************************************/
 /*** MODE 'ScmpS': Comparison between two shapes  ***/
 /****************************************************/
 if (strcmp(MODE,"ScmpS")==0){ 
    Read_MC_Faces_Object(iobjfileA,&MC_VheadA,&MC_FheadA);
    /* Output_MC_Faces_Object("outA.obj",'w',&MC_VheadA,&MC_FheadA,&ovox,output_comment); */
    Read_MC_Faces_Object(iobjfileB,&MC_VheadB,&MC_FheadB);
    printf("#A Nver %d Nface %d\n", Number_Of_MC_VERTEX(&MC_VheadA), Number_Of_MC_FACE(&MC_FheadA)); 
    printf("#B Nver %d Nface %d\n", Number_Of_MC_VERTEX(&MC_VheadB), Number_Of_MC_FACE(&MC_FheadB)); 
   
    MinMax_XYZ_of_MCVertex(MC_MinA,MC_MaxA,&MC_VheadA); 
    MinMax_XYZ_of_MCVertex(MC_MinB,MC_MaxB,&MC_VheadB); 
 

    for (i=0;i<3;++i){
      if (MC_MinA[i]<MC_MinB[i]) MC_Min[i] = MC_MinA[i];
                            else MC_Min[i] = MC_MinB[i];
      if (MC_MaxA[i]>MC_MaxB[i]) MC_Max[i] = MC_MaxA[i];
                            else MC_Max[i] = MC_MaxB[i];
      ovox.OrigPos[i] = MC_Min[i] - ovox.grid_width;  
      ovox.N[i] = (int)ceil((MC_Max[i]-MC_Min[i] + 2*ovox.grid_width)/ovox.grid_width);
      printf("#[%d] MinA %f MinB %f Min %f MaxA %f MaxB %f Max %f N %d\n",i, MC_MinA[i],MC_MinB[i],MC_Min[i], MC_MaxA[i],MC_MaxB[i],MC_Max[i],ovox.N[i]);
   }
   Malloc_Voxel(&ovox,ovox.N[0],ovox.N[1],ovox.N[2]);
   Initialize_Voxel(&ovox);

   setup_static_inout_variables_for_MCFaces(&MC_FheadA);
   check_inside_MCFace_for_VOXEL(&ovox,&MC_FheadA,1.0);
   free_static_inout_variables_for_MCFaces(&MC_FheadA);

   setup_static_inout_variables_for_MCFaces(&MC_FheadB);
   check_inside_MCFace_for_VOXEL(&ovox,&MC_FheadB,2.0);
   free_static_inout_variables_for_MCFaces(&MC_FheadB);

   Marching_Cube_Tetrahedral(&ovox,&MC_Vhead,&MC_Fhead,&MC_Ehead,2.9999,'C');
   Output_MC_Edges_PDB("AandB.pdb",'w',&MC_Vhead,&MC_Ehead,&ovox,1.0,output_comment,"C",0);
   Output_MC_Faces_VRML("AandB.wrl",'w',&MC_Vhead,&MC_Fhead,&ovox,0.0,1.0,0.0,output_comment);
   Output_MC_Edges_VRML("AandBe.wrl",'w',&MC_Vhead,&MC_Ehead,&ovox,0.0,1.0,0.0,output_comment);
   Free_MC_Face_Stack(&MC_Fhead);
   Free_MC_Edge_Stack(&MC_Ehead);
   Free_MC_Vertex_Stack(&MC_Vhead);
   set_zero_for_higher_threshold_VOXEL(&ovox,2.9999);

   Marching_Cube_Tetrahedral(&ovox,&MC_Vhead,&MC_Fhead,&MC_Ehead,1.9999,'C');
   printf("#MC Nver %d Nface %d\n", Number_Of_MC_VERTEX(&MC_Vhead), Number_Of_MC_FACE(&MC_Fhead)); 
   Output_MC_Faces_VRML("diffB.wrl",'w',&MC_Vhead,&MC_Fhead,&ovox,1.0,0.0,0.0,output_comment);
   Output_MC_Edges_VRML("diffBe.wrl",'w',&MC_Vhead,&MC_Ehead,&ovox,1.0,0.0,0.0,output_comment);
   Free_MC_Face_Stack(&MC_Fhead);
   Free_MC_Edge_Stack(&MC_Ehead);
   Free_MC_Vertex_Stack(&MC_Vhead);
   
   set_zero_for_higher_threshold_VOXEL(&ovox,1.9999);
   Marching_Cube_Tetrahedral(&ovox,&MC_Vhead,&MC_Fhead,&MC_Ehead,0.9999,'C');
   /* Output_MC_Edges_PDB("diffA.pdb",'w',&MC_Vhead,&MC_Ehead,&ovox,1.0,output_comment,"C",0); */
   Output_MC_Faces_VRML("diffA.wrl",'w',&MC_Vhead,&MC_Fhead,&ovox,0.0,0.0,1.0,output_comment);
   Output_MC_Edges_VRML("diffAe.wrl",'w',&MC_Vhead,&MC_Ehead,&ovox,0.0,0.0,1.0,output_comment);

}


 /******************************************************/
 /*** MODE 'Ahel': Atom reference helical transform  ***/
 /******************************************************/
 if (strcmp(MODE,"Ahel")==0){ 
   Find_Filename_Core(corename,ipdbfile);
   if (ChainID!='-') 
    { sprintf(buff,"%s",corename); sprintf(corename,"%s%c",buff,ChainID); }
  
   Read_PDB_File(ipdbfile,&HeadAtom,AtomSelect,ChainID);
   Assign_Radius(&HeadAtom,&HeadRadius);
   Natom = Number_Of_Atom(&HeadAtom);
   printf("Natom %d\n",Natom); 

   Read_PDB_File(itarpdbfile,&HeadAtomTar,AtomSelect,ChainID);
   Assign_Radius(&HeadAtomTar,&HeadRadius);
   Natomref = Number_Of_Atom(&HeadAtomTar);
   printf("Natomref %d\n",Natomref); 
   rmsd = Calculate_CRMS_Bwn_Two_ATOMs_Using_CA_With_Same_Rnum(&HeadAtom,&HeadAtomTar,gorig,gref,Rmat);

   if (opdbfile[0]=='\0') sprintf(opdbfile,"helix.pdb");

   if (opdbfile[0] != '\0'){
     sprintf(output_comment,"HELICAL_SYMMETY_UNIT %d",1);
     Write_PDB_File(opdbfile,&HeadAtom,'w',ABCstr[0],output_comment,PAR.COMMAND);

     for (n=1;n<Nunit_helix;++n){
       Transform_Atoms_By_Rmat_Gorig_Gnew(&HeadAtom, Rmat,gorig, gref);
       sprintf(output_comment,"HELICAL_SYMMETY_UNIT %d",n+1);
       Write_PDB_File(opdbfile,&HeadAtom,'a',ABCstr[n],output_comment,PAR.COMMAND);
     }  
   }

   if ((igmmfile[0] != '\0') && (ogmmfile[0] != '\0')){ 
     Ngauss_malloc = Number_Of_Gaussian3D_in_File(igmmfile);
     GaussArray = (struct GAUSS3D *)malloc(sizeof(struct GAUSS3D)*Ngauss_malloc); 
     chain = Read_Gaussian3D_File(igmmfile,Ngauss_malloc,&Ngauss,GaussArray);
     Write_Gaussian3D_File(ogmmfile,'w',Ngauss,GaussArray,ABCstr[0],"");
     for (n=1;n<Nunit_helix;++n){
       Transform_Gaussians_By_Rmat_Gorig_Gnew(Ngauss, GaussArray, Rmat,gorig, gref);
       Write_Gaussian3D_File(ogmmfile,'a',Ngauss,GaussArray,ABCstr[n],"");
     } 
   }
 
 }

 if (PAR.ologfile[0]!='\0') fclose(PAR.fp_olog);
 
 return(1);

} /* end main() */ 
