/*

 <MapCCP4.c>

 Functions for changing CCP4 map file to a simple text file


>> FILE FORMAT CCP4 MAP DATA <<
 ( taken from http://www.ccp4.ac.uk/html/maplib.html#description )
 The term "word" may be 4 characters in following text.


2) DETAILED DESCRIPTION OF THE MAP FORMAT

The overall layout of the file is as follows:

    File header (256 longwords)
    Symmetry information
    Map, stored as a 3-dimensional array 

The header is organised as 56 words followed by space for ten 80 character text labels as follows:

 
 1      NC              # of Columns    (fastest changing in map)
 2      NR              # of Rows
 3      NS              # of Sections   (slowest changing in map)
 4      MODE            Data type
                          0 = envelope stored as signed bytes (from
                              -128 lowest to 127 highest)
                          1 = Image     stored as Integer*2
                          2 = Image     stored as Reals
                          3 = Transform stored as Complex Integer*2
                          4 = Transform stored as Complex Reals
                          5 == 0	
 
                          Note: Mode 2 is the normal mode used in
                                the CCP4 programs. Other modes than 2 and 0
                                may NOT WORK
 5      NCSTART         Number of first COLUMN  in map
 6      NRSTART         Number of first ROW     in map
 7      NSSTART         Number of first SECTION in map
 8      NX              Number of intervals along X
 9      NY              Number of intervals along Y
10      NZ              Number of intervals along Z
11      X length        Cell Dimensions (Angstroms)
12      Y length                     "
13      Z length                     "
14      Alpha           Cell Angles     (Degrees)
15      Beta                         "
16      Gamma                        "
17      MAPC            Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)
18      MAPR            Which axis corresponds to Rows   (1,2,3 for X,Y,Z)
19      MAPS            Which axis corresponds to Sects. (1,2,3 for X,Y,Z)
20      AMIN            Minimum density value
21      AMAX            Maximum density value
22      AMEAN           Mean    density value    (Average)
23      ISPG            Space group number
24      NSYMBT          Number of bytes used for storing symmetry operators
25      LSKFLG          Flag for skew transformation, =0 none, =1 if foll
26-34   SKWMAT          Skew matrix S (in order S11, S12, S13, S21 etc) if
                        LSKFLG .ne. 0.
35-37   SKWTRN          Skew translation t if LSKFLG .ne. 0.
                        Skew transformation is from standard orthogonal
                        coordinate frame (as used for atoms) to orthogonal
                        map frame, as
 
                                Xo(map) = S * (Xo(atoms) - t)
  
38      future use       (some of these are used by the MSUBSX routines
 .          "              in MAPBRICK, MAPCONT and FRODO)
 .          "   (all set to zero by default)
 .          "
52          "

53	MAP	        Character string 'MAP ' to identify file type
54	MACHST		Machine stamp indicating the machine type
			which wrote file
55      ARMS            Rms deviation of map from mean density
56      NLABL           Number of labels being used
57-256  LABEL(20,10)    10  80 character text labels (ie. A4 format)


Symmetry records follow - if any - stored as text as in International Tables, operators separated by * and grouped into 'lines' of 80 characters (i.e. symmetry operators do not cross the ends of the 80-character 'lines' and the 'lines' do not terminate in a *).

Map data array follows.

Note on the machine stamp: The machine stamp (word 54) is a 32-bit quantity containing a set of four `nibbles' (half-bytes) - only half the space is used. Each nibble is a number specifying the representation of (in C terms) double (d), float (f), int (i) and unsigned char (c) types. Thus each stamp is of the form 0xdfic0000. For little endian hardware the stamp is 0x44, 0x41, 0x00, 0x00 while the big endian stamp is 0x11, 0x11, 0x00, 0x00. 


---------------------------------------------------------------------------------------
                                                                                                            
>> FILE FORMAT "MRC" MAP DATA <<
                                                                                                            
The following information is taken from
  http://splorg.org:8080/~tobin/projects/downing/mrc/mrc_format.html
  http://ami.scripps.edu/prtl_data/mrc_specification.htm
                                                                                                            
* The term "word" may be 4 characters in following text.
* Basically, it is the same format as the CCP4 format until the NSYMBT word (the word 24).


 format for MRC image files
                                                                                                            
*                Map/Image Header Format for imsubs2000                 *
*        Length = 1024 bytes, organized as 56 LONG words followed       *
*       by space for 10 80 byte text labels.                            *
*                                                                       *
*       1       NX      number of columns (fastest changing in map)     *
*       2       NY      number of rows                                  *
*       3       NZ      number of sections (slowest changing in map)    *
*       4       MODE    data type :                                     *
*                       0       image : signed 8-bit bytes range -128   *
*                                       to 127                          *
*                       1       image : 16-bit halfwords                *
*                       2       image : 32-bit reals                    *
*                       3       transform : complex 16-bit integers     *
*                       4       transform : complex 32-bit reals        *
*       5       NXSTART number of first column in map (Default = 0)     *
*       6       NYSTART number of first row in map       "              *
*       7       NZSTART number of first section in map   "              *
*       8       MX      number of intervals along X                     *
*       9       MY      number of intervals along Y                     *
*       10      MZ      number of intervals along Z                     *
*       11-13   CELLA   cell dimensions in angstroms                    *
*       14-16   CELLB   cell angles in degrees                          *
*       17      MAPC    axis corresp to cols (1,2,3 for X,Y,Z)          *
*       18      MAPR    axis corresp to rows (1,2,3 for X,Y,Z)          *
*       19      MAPS    axis corresp to sections (1,2,3 for X,Y,Z)      *
*       20      DMIN    minimum density value                           *
*       21      DMAX    maximum density value                           *
*       22      DMEAN   mean density value                              *
*       23      ISPG    space group number 0 or 1 (default=0)           *
*       24      NSYMBT  number of bytes used for symmetry data (0 or 80)*
*       25-49   EXTRA   extra space used for anything  - 0 by default   *
*       50-52   ORIGIN  origin in X,Y,Z used for transforms             *
*       53      MAP     character string 'MAP ' to identify file type   *
*       54      MACHST  machine stamp                                   *
*       55      RMS     rms deviation of map from mean density          *
*       56      NLABL   number of labels being used                     *
*       57-256  LABEL(20,10) 10 80-character text labels                *
*                                                                       *
*       Symmetry records follow - if any - stored as text as in         *
*       International Tables, operators separated by * and grouped into *
*       'lines' of 80 characters (ie. symmetry operators do not cross   *
*       the ends of the 80-character 'lines' and the 'lines' do not     *
*       terminate in a *).                                              *
*                                                                       *
*       Data records follow.                                            *
---------------------------------------------------------------------------------------

>> Difference between *.map file and *.mrc file.

>>> CCP4 <<
24      NSYMBT          Number of bytes used for storing symmetry operators
25      LSKFLG          Flag for skew transformation, =0 none, =1 if foll
26-34   SKWMAT          Skew matrix S (in order S11, S12, S13, S21 etc) if
                        LSKFLG .ne. 0.
35-37   SKWTRN          Skew translation t if LSKFLG .ne. 0.
38      future use       (some of these are used by the MSUBSX routines
52          "

53	MAP	        Character string 'MAP ' to identify file type

>>> MRC <<
*       24      NSYMBT  number of bytes used for symmetry data (0 or 80)*
*       25-49   EXTRA   extra space used for anything  - 0 by default   *
*       50-52   ORIGIN  origin in X,Y,Z used for transforms             *
*       53      MAP     character string 'MAP ' to identify file type   *
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
#include <arpa/inet.h>
#include "Voxel.h" 

/*** FUNCTIONS (GLOBAL) ***/
void Read_Density_File();
int  Write_MapCCP4_File();
int  Write_MapCCP4_File_HtoN();
int  Write_MapCCP4_File_NtoH();

/*** FUNCTIONS (LOCAL) ***/
static void Read_String();
static int   Read_Int();
static float Read_Float();
static short Read_Short();
static char  Read_Char();
static void Inverse_4bytes();
static void Inverse_2bytes();
static void Show_4bytes();


void Read_Density_File(ifname,vox,stdout_type,maxMbyte)
  char   *ifname; 
  struct VOXEL *vox;
  char   stdout_type;  /* if (stdout_type=='F'),  not stdout output */
  double maxMbyte;    /* maximum memory for the image. if (maxMbyte < 0.0), no limit for memory. */
{
 FILE *fp;
 char Order,buff[128];
 int L,i,j,x,y,z,malloc_ok;
 int NC,NR,NS,MODE,NCSTART,NRSTART,NSSTART,NX,NY,NZ,MAPC,MAPR,MAPS,ISPG,NSYMBT,LSKFLG,NLABL;
 float Xlength, Ylength, Zlength,Alpha,Beta,Gamma,AMIN,AMAX,AMEAN,ARMS;
 float SKWMAT[3][3],SKWTRN[3],ORIGIN[3];
 char MAP[5],MACHST[5],LABEL[10][81],command[128];
 char  FileType; /* 'C':ccp4 (*.map), 'M':mrc (*.mrc)  */
 double byte, Mbyte; 

 L = strlen(ifname);
      if ((L>=5)&&(ifname[L-4]=='.')&&(ifname[L-3]=='m')&&(ifname[L-2]=='a')&&(ifname[L-1]=='p')) {FileType = 'C';}
 else if ((L>=5)&&(ifname[L-4]=='.')&&(ifname[L-3]=='M')&&(ifname[L-2]=='A')&&(ifname[L-1]=='P')) {FileType = 'C';} 
 else if ((L>=5)&&(ifname[L-4]=='.')&&(ifname[L-3]=='m')&&(ifname[L-2]=='r')&&(ifname[L-1]=='c')) {FileType = 'M';} 
 else if ((L>=5)&&(ifname[L-4]=='.')&&(ifname[L-3]=='M')&&(ifname[L-2]=='R')&&(ifname[L-1]=='C')) {FileType = 'M';} 
 else {FileType = 'C';}

 if (stdout_type!='F') {printf("#Read_Density_File(\"%s\" FileType '%c' maxMbyte %lf)\n",ifname,FileType,maxMbyte);}

 /* Checking byte number for int and float */

 if (sizeof(int)!=4)   {printf("#ERROR:Size of int(%d) is not 4\n",(int)sizeof(int));}
 if (sizeof(float)!=4) {printf("#ERROR:Size of float(%d) is not 4\n",(int)sizeof(float));}


 L = strlen(ifname); 
 if ((L>3)&&(ifname[L-3]=='.')&&(ifname[L-2]=='g')&&(ifname[L-1]=='z')){
   fp = fopen(ifname,"r");
   if (fp==NULL) {printf("#ERROR:Can't open mapfile \"%s\"\n",ifname); exit(1);}
   else {fclose(fp);}
   sprintf(command,"zcat %s",ifname);
   fp = popen(command,"r");
 }
 else{
   fp = fopen(ifname,"r");
 }

 if (fp==NULL) {printf("#ERROR:Can't open mapfile \"%s\"\n",ifname); exit(1);}

 /** Deciding "Byte Order" from the value of NC,**/ 
 Order = 'N'; /* First try Normal byte order */
 NC = Read_Int(fp,Order,"NC"); 
 if ((NC > 10000)||(NC<0)){  /* If NC is too large or negative, try Inverse byte order */
   Order = 'I';
   Inverse_4bytes(&NC); 
  }

 /** Read Headers **/
 NR = Read_Int(fp,Order,"NR"); 
 NS = Read_Int(fp,Order,"NS"); 
 MODE = Read_Int(fp,Order,"MODE"); 
 /* printf("#MODE %d\n",MODE); */

 if ((MODE!=0) && (MODE!=1) && (MODE != 2)){
   printf("#ERROR in density map (%s). MODE=%d is not proper\n",ifname,MODE);
   exit(1);
 }
 NCSTART = Read_Int(fp,Order,"NCSTART"); 
 NRSTART = Read_Int(fp,Order,"NRSTART"); 
 NSSTART = Read_Int(fp,Order,"NSSTART"); 
 NX = Read_Int(fp,Order,"NX"); 
 NY = Read_Int(fp,Order,"NY"); 
 NZ = Read_Int(fp,Order,"NZ"); 
 Xlength = Read_Float(fp,Order,"Xlength"); 
 Ylength = Read_Float(fp,Order,"Ylength"); 
 Zlength = Read_Float(fp,Order,"Zlength"); 
 Alpha = Read_Float(fp,Order,"Alpha"); 
 Beta = Read_Float(fp,Order,"Beta"); 
 Gamma = Read_Float(fp,Order,"Gamma"); 
 if (stdout_type != 'F'){
   printf("ANGLE_Alpha %f\n",Alpha);
   printf("ANGLE_Beta  %f\n",Beta);
   printf("ANGLE_Gamma %f\n",Gamma);
 }
 if ((Alpha !=90.0) || (Beta !=90.0) || (Gamma !=90.0)){
   printf("#ERROR: density map(%s) has non-orthogonal axis angles (Alpha %f Beta %f Gamma %f)\n",ifname,Alpha,Beta,Gamma); 
   exit(1);
 }
 MAPC = Read_Int(fp,Order,"MAPC"); 
 MAPR = Read_Int(fp,Order,"MAPR"); 
 MAPS = Read_Int(fp,Order,"MAPS"); 
 AMIN = Read_Float(fp,Order,"AMIN"); 
 AMAX = Read_Float(fp,Order,"AMAX"); 
 AMEAN = Read_Float(fp,Order,"AMEAN"); 
 ISPG = Read_Int(fp,Order,"ISPG"); 
 NSYMBT = Read_Int(fp,Order,"NSYMBT"); 

 /** CCP4 **/
 if (FileType=='C'){
   LSKFLG = Read_Int(fp,Order,"LSKFLG"); 
   for (i=0;i<3;++i){ 
     for (j=0;j<3;++j){ 
       sprintf(buff,"SKWMAT%d%d",i+1,j+1); 
       SKWMAT[i][j] = Read_Float(fp,Order,buff);
     }
   }

    for (i=0;i<3;++i){ 
      sprintf(buff,"SKWTRN%d",i+1); 
      SKWTRN[i] = Read_Float(fp,Order,buff);
    }
    for (i=0;i<15;++i){Read_Int(fp,Order,"-");} 

    ORIGIN[0] = ORIGIN[1] = ORIGIN[2] = 0.0;
 }

 /** MRC **/
  else if (FileType=='M'){
     for (i=0;i<25;++i){ Read_Int(fp,Order,"EXTRA");}
     ORIGIN[0] =  Read_Float(fp,Order,"ORIGIN_X"); 
     ORIGIN[1] =  Read_Float(fp,Order,"ORIGIN_Y"); 
     ORIGIN[2] =  Read_Float(fp,Order,"ORIGIN_Z"); 
  }
  else{
    printf("#IMPROPER FILETYPE '%c' for Read_Density_File().\n",FileType);
    exit(1);
  }

  fread(MAP,1,4,fp); MAP[4] = '\0';
  fread(MACHST,1,4,fp); MACHST[4] = '\0';
  
  ARMS = Read_Float(fp,Order,"ARMS"); 
  NLABL = Read_Int(fp,Order,"NLABL"); 

  for (i=0;i<10;++i){
    sprintf(buff,"LABEL%d",i+1); 
    Read_String(fp,LABEL[i],80,buff);
  }

  /** Malloc Voxel **/
  vox->N[0] = NC;   vox->N[1] = NR;   vox->N[2] = NS;
  if ((vox->N[0]<=0) || (vox->N[1]<=0) || (vox->N[2]<=0)){
    printf("#ERROR(voxel_malloc): density map ('%s') has negative voxel size(%d %d %d).\n",
       ifname,vox->N[0],vox->N[1],vox->N[2]); 
    exit(1); 
  }

  byte = (double)vox->N[0]*(double)vox->N[1]*(double)vox->N[2]*(double)sizeof(float);
  Mbyte = byte/1024/1024;
  if ((maxMbyte>0.0) && (Mbyte > maxMbyte)){
     printf("#ERROR:memory %lf Mbyte (%d %d %d)is over maxMbyte %lf.\n",Mbyte,vox->N[0],vox->N[1],vox->N[2],maxMbyte);
     fclose(fp);
     exit(1);
  }
 
  malloc_ok = Malloc_Voxel(vox,vox->N[0],vox->N[1],vox->N[2]); 
  if (malloc_ok==0){
    printf("#ERROR(voxel_malloc): density map ('%s') is too large for malloc(%d %d %d).\n",
       ifname,vox->N[0],vox->N[1],vox->N[2]); 
    exit(1); 
  }
  vox->grid_width = (float)Xlength/(float)NX;
  if (stdout_type != 'F'){
   printf("GRID_SIZE_X %d\n",vox->N[0]);
   printf("GRID_SIZE_Y %d\n",vox->N[1]);
   printf("GRID_SIZE_Z %d\n",vox->N[2]);
   printf("GRID_WIDTH %f\n",vox->grid_width);
  }
  /*
  printf("#voxN %d %d %d grid_width %f\n",vox->N[0],vox->N[1],vox->N[2],vox->grid_width);
  */
  vox->OrigPos[0] = NCSTART * vox->grid_width + ORIGIN[0];
  vox->OrigPos[1] = NRSTART * vox->grid_width + ORIGIN[1];
  vox->OrigPos[2] = NSSTART * vox->grid_width + ORIGIN[2];
  if (stdout_type != 'F'){
    printf("NCSTART %d\n",NCSTART);
    printf("NRSTART %d\n",NRSTART);
    printf("NSSTART %d\n",NSSTART);
    printf("ORIGIN0 %f\n",ORIGIN[0]);
    printf("ORIGIN1 %f\n",ORIGIN[1]);
    printf("ORIGIN2 %f\n",ORIGIN[2]);
    printf("OrigPos0 %f\n",vox->OrigPos[0]);
    printf("OrigPos1 %f\n",vox->OrigPos[1]);
    printf("OrigPos2 %f\n",vox->OrigPos[2]);
  }  


  /** Read NSYMBT characters **/
  for (i=0;i<NSYMBT;++i) fgetc(fp);
 
  /** Read Voxel **/
  for (z=0;z<vox->N[2];++z){
    for (y=0;y<vox->N[1];++y){
      for (x=0;x<vox->N[0];++x){ 
              if (MODE==0){ vox->dat[x][y][z] =  (float)Read_Char(fp,Order,"-");}
         else if (MODE==1){ vox->dat[x][y][z] =  (float)Read_Short(fp,Order,"-");}
         else if (MODE==2){ vox->dat[x][y][z] =  Read_Float(fp,Order,"-");}
        /* vox->dat[x][y][z] =  Read_Float(fp,'I',"-"); */ 
        /*
        if (vox->dat[x][y][z] != 0.0){
        printf("#%e isnan %d #isfinite %d isinf %d isnormal %d\n",
          vox->dat[x][y][z],
          isnan(vox->dat[x][y][z]), isfinite(vox->dat[x][y][z]), isinf(vox->dat[x][y][z]), isnormal(vox->dat[x][y][z]));
        */
          if ((isnan(vox->dat[x][y][z])!=0)||(isinf(vox->dat[x][y][z])!=0)){
            printf("#ERROR(density map:'%s'): density[%d][%d][%d]:'%e' is not normal number !!\n",ifname,x,y,z,vox->dat[x][y][z]);
            printf("#byte::"); Show_4bytes(&(vox->dat[x][y][z]));  printf("\n");
            exit(1);
          }
       }
    }
  } 
 
   fclose(fp);

} /* end of Read_Density_File() */




int Write_MapCCP4_File(vox,fname)
  struct VOXEL *vox;
  char   *fname;
{
 FILE *fp;
 int i,x,y,z;
 float v;
 int NCSTART,NRSTART,NSSTART;
 float AMIN,AMAX,AMEAN,ARMS;
 char MAP[5],MACHST[5],LABEL[10][81];
 double SUM;

 

 if (strcmp(fname,"-")==0){
   fp = stdout;
 }
 else{
   fp = fopen(fname,"w");
   if (fp==NULL) {printf("#ERROR:Can't write to CCP4file \"%s\"\n",fname); return(0);}
   printf("#Write_MapCCP4_File(%d %d %d)-->\"%s\"\n",vox->N[0],vox->N[1],vox->N[2],fname);
 }
 /** Calculate AMIN, AMAX, AMEAN,ARMS **/
 AMIN = AMAX = vox->dat[0][0][0];
 SUM = 0.0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){ 
       v = vox->dat[x][y][z];
       if (v<AMIN) v = AMIN;
       if (v>AMAX) v = AMAX;
       SUM += v;
     } 
   }
 }
 
  AMEAN = SUM/vox->N[0]/vox->N[1]/vox->N[2]; 

 SUM = 0.0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){ 
       v = vox->dat[x][y][z];
       SUM += (v-AMEAN)*(v-AMEAN);
     }
   }
 }  
 
 ARMS = SUM/vox->N[0]/vox->N[1]/vox->N[2]; 
 if (ARMS>0.0) ARMS = sqrt(ARMS);

 /* NC, NR, NS */
 fwrite(&(vox->N[0]),1,4,fp); 
 fwrite(&(vox->N[1]),1,4,fp); 
 fwrite(&(vox->N[2]),1,4,fp); 
 /* Mode */
 i = 2; fwrite(&i,1,4,fp); 
 /* NCSTART, NRSTART, NSSTART */
 NCSTART = (int)(vox->OrigPos[0]/vox->grid_width);
 fwrite(&NCSTART,1,4,fp); 
 NRSTART = (int)(vox->OrigPos[1]/vox->grid_width);
 fwrite(&NRSTART,1,4,fp); 
 NSSTART = (int)(vox->OrigPos[2]/vox->grid_width);
 fwrite(&NSSTART,1,4,fp); 

 
 /* NX, NY, NZ */
 fwrite(&(vox->N[0]),1,4,fp); 
 fwrite(&(vox->N[1]),1,4,fp); 
 fwrite(&(vox->N[2]),1,4,fp); 

 /* X length, Y length, Z length */
 v = vox->N[0]*vox->grid_width; fwrite(&v,1,4,fp); 
 v = vox->N[1]*vox->grid_width; fwrite(&v,1,4,fp); 
 v = vox->N[2]*vox->grid_width; fwrite(&v,1,4,fp); 
 
 /* Alpha, Beta, Gamma */
 v = 90.0;
 fwrite(&v,1,4,fp); 
 fwrite(&v,1,4,fp); 
 fwrite(&v,1,4,fp); 

 /* MAPC, MAPR, MAPS */
 i = 1; fwrite(&i,1,4,fp); 
 i = 2; fwrite(&i,1,4,fp); 
 i = 3; fwrite(&i,1,4,fp); 

 /* AMIN, AMAX, AMEAN */
 fwrite(&AMIN,1,4,fp); 
 fwrite(&AMAX,1,4,fp); 
 fwrite(&AMEAN,1,4,fp); 

 /* ISPG, NSYMBT, LSKFLG */
 i = 1; fwrite(&i,1,4,fp); 
 i = 0; fwrite(&i,1,4,fp); 
 i = 0; fwrite(&i,1,4,fp); 
 
 /* SKWMAT11, SKWMAT12, ..., SKWMAT33 */
 v = 0.0;
 for (x=1;x<=3;++x){
  for (y=1;y<=3;++y){ fwrite(&v,1,4,fp); }
 } 
 /* SKWTRN1, SKWTRN2, SKWTRN2 */
 fwrite(&(vox->OrigPos[0]),1,4,fp); 
 fwrite(&(vox->OrigPos[1]),1,4,fp); 
 fwrite(&(vox->OrigPos[2]),1,4,fp); 
 /*
 printf("vox->min %f %f %f\n",vox->min[0],vox->min[1],vox->min[2]);
 */
 
 /* future use (from 38 to 52 words) */
 i = 0;
 for (x=38;x<=52;++x) fwrite(&i,1,4,fp); 

 /* MAP, MACHST */
 sprintf(MAP,"MAP "); fwrite(MAP,1,4,fp);
 sprintf(MACHST,"DA"); fwrite(MACHST,1,4,fp);
 
 /* ARMS */
 fwrite(&ARMS,1,4,fp); 

 /* NLABL */
 i = 10; fwrite(&i,1,4,fp); 

 /* LABEL */
 for (x=0;x<10;++x){
   for (y=0;y<80;++y){LABEL[x][y] = ' ';}
 }

 sprintf(LABEL[0],"THIS IS MADE BY gmconvert");
 fwrite(LABEL[0],1,80,fp);
 for (x=1;x<10;++x) fwrite(LABEL[x],1,80,fp);
 
 /* write voxel values */

 for (z=0;z<vox->N[2];++z){
  for (y=0;y<vox->N[1];++y){
   for (x=0;x<vox->N[0];++x){ fwrite(&(vox->dat[x][y][z]),1,4,fp);}
  }
 }
 
 if (strcmp(fname,"-")!=0){
   fclose(fp);  
 }                 
 return(1);
} /* end of Write_MapCCP4_File() */




int Write_MapCCP4_File_HtoN(vox,fname)
  struct VOXEL *vox;
  char   *fname;
{
 FILE *fp;
 int x,y,z;
 short i,j;
 float v;
 int NCSTART,NRSTART,NSSTART;
 float AMIN,AMAX,AMEAN,ARMS;
 char MAP[5],MACHST[5],LABEL[10][81];
 double SUM;
 float f;
 

 if (strcmp(fname,"-")==0){
   fp = stdout;
 }
 else{
   fp = fopen(fname,"w");
   if (fp==NULL) {printf("#ERROR:Can't write to CCP4file \"%s\"\n",fname); return(0);}
   printf("#Write_MapCCP4_File(%d %d %d)-->\"%s\"\n",vox->N[0],vox->N[1],vox->N[2],fname);
 }
 /** Calculate AMIN, AMAX, AMEAN,ARMS **/
 AMIN = AMAX = vox->dat[0][0][0];
 SUM = 0.0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){ 
       v = vox->dat[x][y][z];
       if (v<AMIN) v = AMIN;
       if (v>AMAX) v = AMAX;
       SUM += v;
     } 
   }
 }
 
  AMEAN = SUM/vox->N[0]/vox->N[1]/vox->N[2]; 

 SUM = 0.0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){ 
       v = vox->dat[x][y][z];
       SUM += (v-AMEAN)*(v-AMEAN);
     }
   }
 }  
 
 ARMS = SUM/vox->N[0]/vox->N[1]/vox->N[2]; 
 if (ARMS>0.0) ARMS = sqrt(ARMS);

 /* NC, NR, NS */
 i = htons(vox->N[0]); fwrite(&i,1,4,fp);
 i = htons(vox->N[1]); fwrite(&i,1,4,fp);
 i = htons(vox->N[2]); fwrite(&i,1,4,fp);
 /* Mode */
 i = 2; fwrite(&i,1,4,fp); 
 /* NCSTART, NRSTART, NSSTART */
 NCSTART = (int)(vox->OrigPos[0]/vox->grid_width);
 i = htons(NCSTART); fwrite(&i,1,4,fp); 
 NRSTART = (int)(vox->OrigPos[1]/vox->grid_width);
 i = htons(NRSTART); fwrite(&i,1,4,fp); 
 NSSTART = (int)(vox->OrigPos[2]/vox->grid_width);
 i = htons(NSSTART); fwrite(&i,1,4,fp); 

 
 /* NX, NY, NZ */
 i = htons(vox->N[0]); fwrite(&i,1,4,fp);
 i = htons(vox->N[1]); fwrite(&i,1,4,fp);
 i = htons(vox->N[2]); fwrite(&i,1,4,fp);

 /* X length, Y length, Z length */
 v = vox->N[0]*vox->grid_width; f = htons(v); fwrite(&f,1,4,fp); 
 v = vox->N[1]*vox->grid_width; f = htons(v); fwrite(&f,1,4,fp); 
 v = vox->N[2]*vox->grid_width; f = htons(v); fwrite(&f,1,4,fp); 
 
 /* Alpha, Beta, Gamma */
 v = 90.0;
 f = htons(v);
 fwrite(&f,1,4,fp); 
 fwrite(&f,1,4,fp); 
 fwrite(&f,1,4,fp); 

 /* MAPC, MAPR, MAPS */
 j = 1; i = htons(j); fwrite(&i,1,4,fp); 
 j = 2; i = htons(j); fwrite(&i,1,4,fp); 
 j = 3; i = htons(j); fwrite(&i,1,4,fp); 


 /* AMIN, AMAX, AMEAN */
 i = htons(AMIN);  fwrite(&i,1,4,fp);
 i = htons(AMAX);  fwrite(&i,1,4,fp);
 i = htons(AMEAN); fwrite(&i,1,4,fp);

 /* ISPG, NSYMBT, LSKFLG */
 j = 1; i = htons(j); fwrite(&i,1,4,fp); 
 j = 0; i = htons(j); fwrite(&i,1,4,fp); 
 j = 0; i = htons(j); fwrite(&i,1,4,fp); 
 
 /* SKWMAT11, SKWMAT12, ..., SKWMAT33 */
 v = 0.0;
 for (x=1;x<=3;++x){
  for (y=1;y<=3;++y){ 
    f = htons(v);
    fwrite(&f,1,4,fp); 
   }
 } 
 /* SKWTRN1, SKWTRN2, SKWTRN2 */
 f = htons(vox->OrigPos[0]); fwrite(&f,1,4,fp); 
 f = htons(vox->OrigPos[1]); fwrite(&f,1,4,fp); 
 f = htons(vox->OrigPos[2]); fwrite(&f,1,4,fp); 
 /*
 printf("vox->min %f %f %f\n",vox->min[0],vox->min[1],vox->min[2]);
 */
 
 /* future use (from 38 to 52 words) */
 j = 0;
 for (x=38;x<=52;++x){
   i = htons(j);  
   fwrite(&i,1,4,fp); 
 }
 /* MAP, MACHST */

 sprintf(MAP,"MAP "); fwrite(MAP,1,4,fp);
 sprintf(MACHST,"DA"); fwrite(MACHST,1,4,fp);
 
 /* ARMS */
 i = htons(ARMS); fwrite(&i,1,4,fp); 

 /* NLABL */
 j = 10; i = htons(j); fwrite(&i,1,4,fp); 

 /* LABEL */
 for (x=0;x<10;++x){
   for (y=0;y<80;++y){LABEL[x][y] = ' ';}
 }

 sprintf(LABEL[0],"THIS IS MADE BY gmconvert");
 fwrite(LABEL[0],1,80,fp);
 for (x=1;x<10;++x){ fwrite(LABEL[x],1,80,fp);}
 
 /* write voxel values */

 for (z=0;z<vox->N[2];++z){
   for (y=0;y<vox->N[1];++y){
     for (x=0;x<vox->N[0];++x){ 
       f = htons(vox->dat[x][y][z]);
       fwrite(&f,1,4,fp);
     }
   }
 }
 
 if (strcmp(fname,"-")!=0){
   fclose(fp);  
 }                 
 return(1);
} /* end of Write_MapCCP4_File_HtoN() */


int Write_MapCCP4_File_NtoH(vox,fname)
  struct VOXEL *vox;
  char   *fname;
{
 FILE *fp;
 int x,y,z;
 short i,j;
 float v;
 int NCSTART,NRSTART,NSSTART;
 float AMIN,AMAX,AMEAN,ARMS;
 char MAP[5],MACHST[5],LABEL[10][81];
 double SUM;
 float f;
 

 if (strcmp(fname,"-")==0){
   fp = stdout;
 }
 else{
   fp = fopen(fname,"w");
   if (fp==NULL) {printf("#ERROR:Can't write to CCP4file \"%s\"\n",fname); return(0);}
   printf("#Write_MapCCP4_File(%d %d %d)-->\"%s\"\n",vox->N[0],vox->N[1],vox->N[2],fname);
 }
 /** Calculate AMIN, AMAX, AMEAN,ARMS **/
 AMIN = AMAX = vox->dat[0][0][0];
 SUM = 0.0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){ 
       v = vox->dat[x][y][z];
       if (v<AMIN) v = AMIN;
       if (v>AMAX) v = AMAX;
       SUM += v;
     } 
   }
 }
 
  AMEAN = SUM/vox->N[0]/vox->N[1]/vox->N[2]; 

 SUM = 0.0;
 for (x=0;x<vox->N[0];++x){
   for (y=0;y<vox->N[1];++y){
     for (z=0;z<vox->N[2];++z){ 
       v = vox->dat[x][y][z];
       SUM += (v-AMEAN)*(v-AMEAN);
     }
   }
 }  
 
 ARMS = SUM/vox->N[0]/vox->N[1]/vox->N[2]; 
 if (ARMS>0.0) ARMS = sqrt(ARMS);

 /* NC, NR, NS */
 i = ntohs(vox->N[0]); fwrite(&i,1,4,fp);
 i = ntohs(vox->N[1]); fwrite(&i,1,4,fp);
 i = ntohs(vox->N[2]); fwrite(&i,1,4,fp);
 /* Mode */
 i = 2; fwrite(&i,1,4,fp); 
 /* NCSTART, NRSTART, NSSTART */
 NCSTART = (int)(vox->OrigPos[0]/vox->grid_width);
 i = ntohs(NCSTART); fwrite(&i,1,4,fp); 
 NRSTART = (int)(vox->OrigPos[1]/vox->grid_width);
 i = ntohs(NRSTART); fwrite(&i,1,4,fp); 
 NSSTART = (int)(vox->OrigPos[2]/vox->grid_width);
 i = ntohs(NSSTART); fwrite(&i,1,4,fp); 

 
 /* NX, NY, NZ */
 i = ntohs(vox->N[0]); fwrite(&i,1,4,fp);
 i = ntohs(vox->N[1]); fwrite(&i,1,4,fp);
 i = ntohs(vox->N[2]); fwrite(&i,1,4,fp);

 /* X length, Y length, Z length */
 v = vox->N[0]*vox->grid_width; f = ntohs(v); fwrite(&f,1,4,fp); 
 v = vox->N[1]*vox->grid_width; f = ntohs(v); fwrite(&f,1,4,fp); 
 v = vox->N[2]*vox->grid_width; f = ntohs(v); fwrite(&f,1,4,fp); 
 
 /* Alpha, Beta, Gamma */
 v = 90.0;
 f = ntohs(v);
 fwrite(&f,1,4,fp); 
 fwrite(&f,1,4,fp); 
 fwrite(&f,1,4,fp); 

 /* MAPC, MAPR, MAPS */
 j = 1; i = ntohs(j); fwrite(&i,1,4,fp); 
 j = 2; i = ntohs(j); fwrite(&i,1,4,fp); 
 j = 3; i = ntohs(j); fwrite(&i,1,4,fp); 


 /* AMIN, AMAX, AMEAN */
 i = ntohs(AMIN);  fwrite(&i,1,4,fp);
 i = ntohs(AMAX);  fwrite(&i,1,4,fp);
 i = ntohs(AMEAN); fwrite(&i,1,4,fp);

 /* ISPG, NSYMBT, LSKFLG */
 j = 1; i = ntohs(j); fwrite(&i,1,4,fp); 
 j = 0; i = ntohs(j); fwrite(&i,1,4,fp); 
 j = 0; i = ntohs(j); fwrite(&i,1,4,fp); 
 
 /* SKWMAT11, SKWMAT12, ..., SKWMAT33 */
 v = 0.0;
 for (x=1;x<=3;++x){
  for (y=1;y<=3;++y){ 
    f = ntohs(v);
    fwrite(&f,1,4,fp); 
   }
 } 
 /* SKWTRN1, SKWTRN2, SKWTRN2 */
 f = ntohs(vox->OrigPos[0]); fwrite(&f,1,4,fp); 
 f = ntohs(vox->OrigPos[1]); fwrite(&f,1,4,fp); 
 f = ntohs(vox->OrigPos[2]); fwrite(&f,1,4,fp); 
 /*
 printf("vox->min %f %f %f\n",vox->min[0],vox->min[1],vox->min[2]);
 */
 
 /* future use (from 38 to 52 words) */
 j = 0;
 for (x=38;x<=52;++x){
   i = ntohs(j);  
   fwrite(&i,1,4,fp); 
 }
 /* MAP, MACHST */

 sprintf(MAP,"MAP "); fwrite(MAP,1,4,fp);
 sprintf(MACHST,"DA"); fwrite(MACHST,1,4,fp);
 
 /* ARMS */
 i = ntohs(ARMS); fwrite(&i,1,4,fp); 

 /* NLABL */
 j = 10; i = ntohs(j); fwrite(&i,1,4,fp); 

 /* LABEL */
 for (x=0;x<10;++x){
   for (y=0;y<80;++y){LABEL[x][y] = ' ';}
 }

 sprintf(LABEL[0],"THIS IS MADE BY gmconvert");
 fwrite(LABEL[0],1,80,fp);
 for (x=1;x<10;++x){ fwrite(LABEL[x],1,80,fp);}
 
 /* write voxel values */

 for (z=0;z<vox->N[2];++z){
   for (y=0;y<vox->N[1];++y){
     for (x=0;x<vox->N[0];++x){ 
       f = ntohs(vox->dat[x][y][z]);
       fwrite(&f,1,4,fp);
     }
   }
 }
 
 if (strcmp(fname,"-")!=0){
   fclose(fp);  
 }                 
 return(1);
} /* end of Write_MapCCP4_File_NtoH() */
























void Read_String(fp,instr,length,str)
 FILE *fp;
 char *instr;
 int  length;
 char *str;
{ int i;
  int ret;

  ret = fread(instr,1,length,fp); 
  if (ret<length) {printf("#ERROR:Can't read %d chars.\n",length); exit(1);}
  for (i=0;i<length;++i)
    if (iscntrl(instr[i])!=0) instr[i] = ' ';
  instr[length] = '\0';
  /* if (str[0]!='-') printf("#%s \"%s\"\n",str,instr); */

} /* end of Read_String() */



int Read_Int(fp,Order,comment)
 FILE *fp;
 char Order; /* 'I'nverse */
 char *comment;
{ int i;
  int ret; 
 ret = fread(&i,1,4,fp); 
 if (ret<4) {printf("#ERROR Read_Int(%s):Can't read 4 chars.\n",comment); exit(1);}
 if (Order=='I'){Inverse_4bytes(&i);}
 /* if (comment[0]!='-') printf("#%s %d\n",comment,i); */
 return(i); 
} /* end of Read_Int() */



short Read_Short(fp,Order,comment)
 FILE *fp;
 char Order; /* 'I'nverse */
 char *comment;
/* short :: 2 bytes(16bits) */
{ short i;
  int ret; 

  ret = fread(&i,1,2,fp); 
  if (ret<2) {printf("#ERROR Read_Short(%s):Can't read 2 chars.\n",comment); exit(1);}
   if (Order=='I'){Inverse_2bytes(&i);} 
 /* if (comment[0]!='-') printf("#%s %d\n",comment,i); */
 return(i); 
} /* end of Read_Short() */


char Read_Char(fp,Order,comment)
 FILE *fp;
 char Order; /* 'I'nverse */
 char *comment;
/* short :: 2 bytes(16bits) */
{ char i;
  int ret; 

  ret = fread(&i,1,1,fp); 
  if (ret<1) {printf("#ERROR:Can't read 1 char.\n"); exit(1);}
 return(i); 
} /* end of Read_Char() */




float Read_Float(fp,Order,comment)
 FILE *fp;
 char Order; /* 'I'nverse */
 char *comment;   
{ float f;
  int ret; 

 ret = fread(&f,1,4,fp); 
 if (ret<4) {printf("#ERROR Read_Float(%s):Can't read 4 chars.\n",comment); exit(1);}
 /* Show_4bytes(&f);  printf(" %e\n",f);  */
 if (Order=='I'){Inverse_4bytes(&f);}
 /* if (comment[0] !='-') printf("#%s %f\n",comment,f); */
 return(f); 

} /* end of Read_Float() */




void Inverse_4bytes(X)
 unsigned char X[4];
{
 /* 
 X0 X1 X2 X3 
    ->
 X3 X2 X1 X0 
 */
 unsigned char buff;
 /*
 printf("%x %x %x %x\n",X[0],X[1],X[2],X[3]);
 */ 
 buff = X[3]; X[3] = X[0]; X[0] = buff;
 buff = X[2]; X[2] = X[1]; X[1] = buff;

}  /* end of Inverse_4bytes() */

void Show_4bytes(X)
 unsigned char X[4];
{
 printf("#%x %x %x %x",X[0],X[1],X[2],X[3]); 
}  /* end of Show_4bytes() */



void Inverse_2bytes(X)
 unsigned char X[2];
{
 /* 
 X0 X1  
    ->
 X1 X0 
 */
 unsigned char buff;
 /*
 printf("%x %x %x %x\n",X[0],X[1],X[2],X[3]);
 */ 
 buff =  X[0];
 X[0] =  X[1];
 X[1] =  buff;
}  /* end of Inverse_2bytes() */

