/*
  <transf_atom.c>

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "io_mmCIF.h"
#include "PDB_from_mmCIF.h"

#define MAX_FILE_LEN 1024

static char LastModDate[] = "Feb 17, 2016";

char MMCIF_DIR[MAX_FILE_LEN];
char PDB_GMM_DIR[MAX_FILE_LEN];
char PDB_CHAIN_GMM_DIR[MAX_FILE_LEN];
char EMDB_GMM_DIR[MAX_FILE_LEN];
char TMPOUT_DIR[MAX_FILE_LEN];
char OMOKAGE_USERDATA_DIR[MAX_FILE_LEN];
char SASBDB_GMM_DIR[MAX_FILE_LEN];
char SASBDB_MODEL_SPLITCIF[MAX_FILE_LEN];

char *ARG_STRING;

static int read_environment_file();
static int read_pdb_and_out_transformed_xyz();
static void get_Rmat_and_Tvec_from_strings();
static void setup_filenames();


int read_environment_file(ienvfile)
  char *ienvfile;
{
/*
 ## >> FORMAT EXAMPLE << ##
 # BASE_DIR  /usr/people/takawaba/work/STRMAT/Matras10
 # SCORE_DIR /usr/people/takawaba/work/STRMAT/Matras10/ROM-00Oct6
 # BSSP_DIR  /lab/DB/BSSP
 # PDB_DIR   /lab/PDB
 # OMOKAGE_USERDATA_DIR /home/pdbj/omokage/userdata
 # TMPOUT_DIR    /var/www/html/gmfit/TMPOUT

*/
  FILE *fp;
  char line[1000],key[1000],value[1000];
  int L;

  /* printf("#read_environment_file('%s')\n",ienvfile); */

  MMCIF_DIR[0]    = '\0'; 
  PDB_GMM_DIR[0]  = '\0'; 
  PDB_CHAIN_GMM_DIR[0]  = '\0'; 
  EMDB_GMM_DIR[0] = '\0'; 
  SASBDB_GMM_DIR[0] = '\0'; 
  SASBDB_MODEL_SPLITCIF[0] = '\0'; 

  fp = fopen(ienvfile,"r");
  if (fp==NULL){
    /* printf("#WARNING: Can't open ienvfile '%s'.\n",ienvfile); */
    return(0);
  }

  while (feof(fp)==0){
    line[0] = '\0';
    fgets(line,1000,fp);
    L = strlen(line);
    /* printf("%s",line); */
    if (line[L-1]=='\n'){ 
      line[L-1] = '\0';
      L = strlen(line);
    }

    if ((L>5) && (line[0] != '#')){
       split_into_two_words(line,0,' ',key,value);  
       if (strcmp(key,"MMCIF_DIR")==0){ sprintf(MMCIF_DIR,"%s",value);}
       if (strcmp(key,"PDB_GMM_DIR")==0){ sprintf(PDB_GMM_DIR,"%s",value);}
       if (strcmp(key,"PDB_CHAIN_GMM_DIR")==0){ sprintf(PDB_CHAIN_GMM_DIR,"%s",value);}
       if (strcmp(key,"EMDB_GMM_DIR")==0){ sprintf(EMDB_GMM_DIR,"%s",value);}
       if (strcmp(key,"TMPOUT_DIR")==0){ sprintf(TMPOUT_DIR,"%s",value);}
       if (strcmp(key,"OMOKAGE_USERDATA_DIR")==0){ sprintf(OMOKAGE_USERDATA_DIR,"%s",value);}
       if (strcmp(key,"SASBDB_GMM_DIR")==0){ sprintf(SASBDB_GMM_DIR,"%s",value);}
       if (strcmp(key,"SASBDB_MODEL_SPLITCIF")==0){ sprintf(SASBDB_MODEL_SPLITCIF,"%s",value);}
    }
  }

  fclose(fp);
  return(1);
} /* end of read_environment_file() */



int read_pdb_and_out_transformed_xyz(ifname,ofname,chainid_specified,residue_select,TransformType,Rmat,Tvec,My_tFactor)
  char *ifname;
  char *ofname;
  char chainid_specified;
  char residue_select;    /* if residue_select = 'T', then write only ' CA ' andn ' P  '. */
  char  TransformType;    /* if TransformType=='T', transform atomic XYZ by Rmat[][] and Tvec[] */
  float Rmat[3][3],Tvec[3];
  float My_tFactor;
{
/*
#          1         2         3         4         5         6         7
#01234567890123456789012345678901234567890123456789012345678901234567890123456789
#ATOM    676  CB  GLU A  85      10.440  29.552  12.788  6.00 16.96           C
#ATOM    680  OE2 GLU A  85      10.230  30.451  16.374  8.00 41.03           O
#ATOM    682  CA  LEU A  86       7.618  29.487   9.238  6.00 12.23           C
#HETATM 1236  O4  SO4   154      33.810  28.815  -4.624  8.00 14.90           O
#HETATM 1237 FE   HEM   155      15.271  27.962   0.622 24.00  7.86          FE
#REMARK   1
#REMARK   1 REFERENCE 1
#REMARK   1  AUTH   T.ALBER,D.W.BANNER,A.C.BLOOMER,G.A.PETSKO,
#REMARK   1  AUTH 2 D.PHILLIPS,P.S.RIVERS,I.A.WILSON
#REMARK   1  TITL   ON THE THREE-DIMENSIONAL STRUCTURE AND CATALYTIC
#REMARK   1  TITL 2 MECHANISM OF TRIOSE PHOSPHATE ISOMERASE
#REMARK   1  REF    PHILOS.TRANS.R.SOC.LONDON,    V. 293   159 1981
*/
  FILE *fpi,*fpo;
  int L;
  char command[MAX_FILE_LEN],line[MAX_FILE_LEN],head[MAX_FILE_LEN],tail[MAX_FILE_LEN];
  char atomname[8],value[16],Occup_str[16],tFactor_str[16],accept;
  float x,y,z,new_x,new_y,new_z,tFactor;

  L = strlen(ifname);
  if ((L>3)&&(ifname[L-3]=='.')&&(ifname[L-2]=='g')&&(ifname[L-1]=='z')){
    sprintf(command,"zcat %s",ifname);
    fpi = popen(command,"r");
  }
  else{
    fpi = fopen(ifname,"r");
  }
  if (fpi==NULL){ return(0);}

  if (ofname[0]=='-'){
    fpo = stdout;
  }
  else{
    fpo = fopen(ofname,"w");
     printf("#read_pdb_and_out_transformed_xyz('%s') --> '%s'",ifname,ofname);
  }

  fprintf(fpo,"REMARK      read_pdb_and_out_transformed_xyz()\n");
  fprintf(fpo,"REMARK      COMMAND transf_atom.cgi %s\n",ARG_STRING);
  fprintf(fpo,"REMARK      R0 %+9.6f %+9.6f %+9.6f T0 %+12.6f\n",Rmat[0][0],Rmat[0][1],Rmat[0][2],Tvec[0]);
  fprintf(fpo,"REMARK      R1 %+9.6f %+9.6f %+9.6f T1 %+12.6f\n",Rmat[1][0],Rmat[1][1],Rmat[1][2],Tvec[1]);
  fprintf(fpo,"REMARK      R2 %+9.6f %+9.6f %+9.6f T2 %+12.6f\n",Rmat[2][0],Rmat[2][1],Rmat[2][2],Tvec[2]);

  while (feof(fpi)==0){
    line[0] = '\0';
    fgets(line,MAX_LENGTH_LINE,fpi);
    L = strlen(line);
    if (line[L-1]=='\n'){
      line[L-1] = '\0';
      L = strlen(line);
    }
    if ((strncmp(line,"HEADER",6)==0)||(strncmp(line,"REMARK",6)==0)){
      fprintf(fpo,"%s\n",line); 
    }
    if ((strncmp(line,"ATOM",4)==0)||(strncmp(line,"HETATM",6)==0)){
      get_part_of_line(atomname,line,12,15);
      accept = 1;
      if (residue_select == 'T'){
        if ((strcmp(atomname," CA ")!=0) && (strcmp(atomname," P  ")!=0)){
        accept = 0;
        }
      }
      if (accept == 1)
/*
#Occup_str  = line[54:60]
#tFactor_str = line[60:66]
#tail        = line[66:]
*/
        get_part_of_line(head,line,0,29);
        get_part_of_line(Occup_str,line,54,59);
        get_part_of_line(tFactor_str,line,60,65);
        get_part_of_line(tail,line,66,L-1);
        tFactor = atof(tFactor_str);
        if (My_tFactor>=0.0){
          tFactor = My_tFactor;
        } 
        get_part_of_line(value,line,30,37); x = atof(value);
        get_part_of_line(value,line,38,45); y = atof(value);
        get_part_of_line(value,line,46,54); z = atof(value);

        if (TransformType=='T'){
          new_x = Rmat[0][0]*x + Rmat[0][1]*y + Rmat[0][2]*z + Tvec[0];
          new_y = Rmat[1][0]*x + Rmat[1][1]*y + Rmat[1][2]*z + Tvec[1];
          new_z = Rmat[2][0]*x + Rmat[2][1]*y + Rmat[2][2]*z + Tvec[2];
          x = new_x;
          y = new_y;
          z = new_z;
        }

        fprintf(fpo,"%s%8.3f%8.3f%8.3f%s%6.2f%s\n",head,x,y,z,Occup_str,tFactor,tail);
      }
    }


  fprintf(fpo,"TER\n");
  fclose(fpi);

  if (ofname[0] != '-'){ 
    fclose(fpo);
  }

  return(1);

} /* end of read_pdb_and_out_transformed_xyz() */
 




void get_Rmat_and_Tvec_from_strings(Rmat,Tvec,Rmat_str,Tvec_str)
  float Rmat[3][3],Tvec[3];
  char *Rmat_str, *Tvec_str;
{
  int Nword,Wsta[100],Wend[100],i,j;
  char value[100];

  split_into_words(Tvec_str,':',&Nword,Wsta,Wend,100);
  for (i=0;i<3;++i){
    get_part_of_line(value,Tvec_str,Wsta[i],Wend[i]);
    Tvec[i] = atof(value);
  }

  split_into_words(Rmat_str,':',&Nword,Wsta,Wend,100);

  for (i=0;i<3;++i){
    for (j=0;j<3;++j){
      get_part_of_line(value,Rmat_str,Wsta[3*i+j],Wend[3*i+j]);
      Rmat[i][j] = atof(value);
    }
  }
} /* end of get_Rmat_and_Tvec_from_strings() */


void setup_filenames(id,type,assembly_id,auth_asym_id,iciffile,ipdbfile,igmmfile)
  char *id;           /* (input) */
  char *type;         /* (input) 'pdb' or 'emdb' or 'gmfit_up_pdb','gmfit_up_map', 'omokage_up_pdb', 'omokage_up_map' */
  char *assembly_id;  /* (input) */
  char *auth_asym_id; /* (input) */
  char *iciffile;     /* (output) */
  char *ipdbfile;     /* (output) */
  char *igmmfile;     /* (output) */
{
  char subdir[100];

  iciffile[0] = '\0';
  ipdbfile[0] = '\0';
  igmmfile[0] = '\0';
 
  if (id[0]!='\0'){
    if (strcmp(type,"pdb")==0){
      get_part_of_line(subdir,id,1,2);
      sprintf(iciffile,"%s/%s/%s.cif.gz",MMCIF_DIR,subdir,id);
      if ((assembly_id[0]!='\0')&&(assembly_id[0]!='-')){
        sprintf(igmmfile,"%s/%s/%s-%s.gmm",PDB_GMM_DIR,subdir,id,assembly_id);
      }
      else if ((auth_asym_id[0]!='\0')&&(auth_asym_id[0]!='-')){
        sprintf(igmmfile,"%s/%s/%s_%s.gmm",PDB_CHAIN_GMM_DIR,subdir,id,auth_asym_id);
      }
      else{
        sprintf(igmmfile,"%s/%s/%s.gmm",PDB_GMM_DIR,subdir,id);
      }
    }
    else if (strcmp(type,"emdb")==0){
      sprintf(igmmfile,"%s/emd_%s.gmm",EMDB_GMM_DIR,id);
    }
    else if (strcmp(type,"sasbdb")==0){
      sprintf(igmmfile,"%s/%s.gmm",SASBDB_GMM_DIR,id);
      sprintf(iciffile,"%s/%s.cif",SASBDB_MODEL_SPLITCIF,id);
    }
    else if (strcmp(type,"gmfit_up_pdb")==0){
      sprintf(igmmfile,"%s/%s.gmm",TMPOUT_DIR,id);
      sprintf(ipdbfile,"%s/%s.pdb",TMPOUT_DIR,id);
    }
    else if (strcmp(type,"gmfit_up_map")==0){
      sprintf(igmmfile,"%s/%s.gmm",TMPOUT_DIR,id);
    }
    else if (strcmp(type,"omokage_up_pdb")==0){
      sprintf(igmmfile,"%s/%s.gmm",OMOKAGE_USERDATA_DIR,id);
      sprintf(ipdbfile,"%s/%s.pdb",OMOKAGE_USERDATA_DIR,id);
    }
    else if (strcmp(type,"omokage_up_map")==0){
      sprintf(igmmfile,"%s/%s.gmm",OMOKAGE_USERDATA_DIR,id);
    }

  }
} /* end of setup_filenames() */


/**************/
/**** MAIN ****/
/**************/

int main(argc,argv)
  int argc;
  char **argv;
{
  struct BLOCK BLKtar,BLKref;
  struct PDB_DATA PDBtar,PDBref;
  struct ASSEMBLY *A; 
  int k;
  char OutType[8],OutGMM,idref[MAX_FILE_LEN],idtar[MAX_FILE_LEN],assembly_id_tar[16],assembly_id_ref[16];
  char auth_asym_id_tar[4],auth_asym_id_ref[4]; 
  char typeref[8],typetar[8],tftar_str[8],tfref_str[8],Rtar_str[MAX_FILE_LEN],Ttar_str[MAX_FILE_LEN],Rgtar_str[MAX_FILE_LEN],Tgtar_str[MAX_FILE_LEN];  
  char iciffile_tar[MAX_FILE_LEN],iciffile_ref[MAX_FILE_LEN];
  char ipdbfile_tar[MAX_FILE_LEN],ipdbfile_ref[MAX_FILE_LEN];
  char igmmfile_tar[MAX_FILE_LEN],igmmfile_ref[MAX_FILE_LEN];
  char download_file[MAX_FILE_LEN];
  char word[128],key[128],value[MAX_FILE_LEN],residue_select; 
  float Rtar[3][3],Ttar[3],Rgtar[3][3],Tgtar[3];   
  float Rref[3][3],Tref[3];
  int Nword,Wsta[100],Wend[100];
  int maxNatom_fullatom, maxNatom_residue;
  int Lquery_string;

  if (read_environment_file("/var/www/html/gmfit/cgi-bin/CONF/webgmfit.conf")==0){
    /* sprintf(MMCIF_DIR,"/DB/mmCIF"); */
    sprintf(MMCIF_DIR,"/kf1/PDBj/ftp/pdbj/pub/pdb/data/structures/divided/mmCIF");
  }  

  maxNatom_fullatom =  100000; 
  maxNatom_residue  =  500000; 

  sprintf(OutType,"plain");
  assembly_id_tar[0] = assembly_id_ref[0] = '\0';
  idref[0] = idtar[0] = '\0';
  sprintf(typeref,"pdb"); 
  sprintf(typetar,"pdb"); 
  sprintf(tfref_str,"0.0");
  sprintf(tftar_str,"1.0");
  sprintf(Rtar_str,"1:0:0:0:1:0:0:0:1");
  sprintf(Rgtar_str,"1:0:0:0:1:0:0:0:1");
  sprintf(Ttar_str,"0:0:0");
  sprintf(Tgtar_str,"0:0:0");
  iciffile_tar[0] = '\0';
  iciffile_ref[0] = '\0';
  auth_asym_id_tar[0] = '\0';
  auth_asym_id_ref[0] = '\0';
  Rref[0][0] = 1.0; Rref[0][1] = 0.0; Rref[0][2] = 0.0;
  Rref[1][0] = 0.0; Rref[1][1] = 1.0; Rref[1][2] = 0.0;
  Rref[2][0] = 0.0; Rref[2][1] = 0.0; Rref[2][2] = 1.0;
  Tref[0] = 0.0; Tref[1] = 0.0; Tref[2] = 0.0;

  OutGMM = 'T';

  if ((argc<2)&&(getenv("QUERY_STRING")==NULL)){
    printf("Content-type: text/html;charset=utf-8\n\n");
    printf("<HTML><BODY><PRE>\n");
    printf("transf_atom.cgi [key1]=[value1] [key2]=[value2] ... \n");
    printf(" for transforming atomic models in PDB or mmCIF format and corresponding GMMs.\n");
    printf(" mainly for the 'pairwise_gmfit' WEB service.\n");   
    printf(" This program is a part of the 'gmconvert' program.\n");   
    printf("  coded by T.Kawabata. LastModDate:%s.\n",LastModDate);
    printf(" ** 'reference' atoms are output with temperature == 0.0.\n");
    printf(" ** 'target'    atoms are output with temperature >  0.0.\n");
    printf("<keys>\n");
    printf(" out              : type of output. 'plain', 'html' or 'download'[%s]\n",OutType);
    printf(" outgmm           : Output corresponding GMM files. ('T' or 'F') [%c]\n",OutGMM);
    printf("<keys for 'reference' (fixed)>\n");
    printf(" idref            : ID (pdb_id, emdb_id, upload_id)  for reference [%s]\n",idref);
    printf(" typeref          : type for reference 'pdb','emdb','sasbdb','gmfit_up_pdb', 'gmfit_up_map', 'omokage_up_map', 'omokage_up_pdb'. [%s]\n",typeref);
    printf(" assembly_id_ref  : assembly_id of biological unit for target. '-':asymmetric unit.[%s]\n",assembly_id_ref);
    printf(" auth_asym_id_ref : auth_asym_id of for target. [%s]\n",auth_asym_id_ref);
    printf(" tfref            : temperature factor for output reference PDB [%s]\n",tfref_str);
    printf("<keys for 'target' (to be transformed)>\n");
    printf(" idtar            : ID (pdb_id, emdb_id, upload_id)  for target [%s]\n",idtar);
    printf(" typetar          : type for target 'pdb','emdb','sasbdb','gmfit_up_pdb', 'gmfit_up_map', 'omokage_up_pdb','omokage_up_map'. [%s]\n",typetar);
    printf(" assembly_id_tar  : assembly_id of biological unit for target. '-':asymmetric unit.[%s]\n",assembly_id_tar);
    printf(" auth_asym_id_tar : auth_asym_id of for target. [%s]\n",auth_asym_id_ref);
    printf(" tftar            : temperature factor for output target PDB [%s]\n",tftar_str);
    printf(" Rtar             : rotation matrix for target[%s]\n",Rtar_str);
    printf(" Ttar             : translational vector for target[%s]\n",Ttar_str);
    printf(" Rgtar            : rotation matrix for target GMM file[%s]\n",Rtar_str);
    printf(" Tgtar            : translational vector for target GMM file[%s]\n",Ttar_str);
    printf("<other keys>\n");
    printf(" natomfull        : maximum atom number for full atom. If over, enforce residue_select = 'T'. [%d]\n",maxNatom_fullatom);
    printf(" natomres         : maximum atom number for residue.   If over, don't show any atoms.         [%d]\n",maxNatom_residue);
    printf("</PRE></BODY></HTML>\n");
    exit(1);
  }
 
 /** [1] Read options **/ 

  if (getenv("QUERY_STRING")!=NULL){
    Lquery_string = strlen(getenv("QUERY_STRING"));
    ARG_STRING = (char *)malloc(sizeof(char)*(Lquery_string+1));
    sprintf(ARG_STRING,"%s",getenv("QUERY_STRING"));
    split_into_words(getenv("QUERY_STRING"),'&',&Nword,Wsta,Wend,100);
  }
  else{
    Nword = argc-1;
    Lquery_string = 0;
    for (k=0;k<argc;++k){ Lquery_string += (strlen(argv[k])+1); }
    ARG_STRING = (char *)malloc(sizeof(char)*(Lquery_string+1));
    ARG_STRING[0] = '\0'; 
    for (k=0;k<argc;++k){ strcat(ARG_STRING,argv[k]); strcat(ARG_STRING," ");}
  }

  /* printf("ARG_STRING:%s\n",ARG_STRING); */

  for (k=0;k<Nword;++k){
    if (getenv("QUERY_STRING")==NULL){
      split_into_two_words(argv[k+1],0,'=',key,value);
      
    }
    else{
      get_part_of_line(word,getenv("QUERY_STRING"),Wsta[k],Wend[k]);
      split_into_two_words(word,0,'=',key,value);
    }
         if (strcmp(key,"out")==0) {sprintf(OutType,"%s",value); }
    else if (strcmp(key,"outgmm")==0) {OutGMM = value[0];}
    else if (strcmp(key,"idref")==0) {sprintf(idref,"%s",value); }
    else if (strcmp(key,"idtar")==0) {sprintf(idtar,"%s",value); }
    else if (strcmp(key,"assembly_id_ref")==0) {sprintf(assembly_id_ref,"%s",value); }
    else if (strcmp(key,"assembly_id_tar")==0) {sprintf(assembly_id_tar,"%s",value); }
    else if (strcmp(key,"auth_asym_id_ref")==0) {sprintf(auth_asym_id_ref,"%s",value); }
    else if (strcmp(key,"auth_asym_id_tar")==0) {sprintf(auth_asym_id_tar,"%s",value); }
    else if (strcmp(key,"typeref")==0) {sprintf(typeref,"%s",value); }
    else if (strcmp(key,"typetar")==0) {sprintf(typetar,"%s",value); }
    else if (strcmp(key,"tfref")==0) {sprintf(tfref_str,"%s",value); }
    else if (strcmp(key,"tftar")==0) {sprintf(tftar_str,"%s",value); }
    else if (strcmp(key,"Rtar")==0) {sprintf(Rtar_str,"%s",value); }
    else if (strcmp(key,"Ttar")==0) {sprintf(Ttar_str,"%s",value); }
    else if (strcmp(key,"Rgtar")==0) {sprintf(Rgtar_str,"%s",value); }
    else if (strcmp(key,"Tgtar")==0) {sprintf(Tgtar_str,"%s",value); }
    else if (strcmp(key,"natomfull")==0) {maxNatom_fullatom = atoi(value); }
    else if (strcmp(key,"natomres")==0)  {maxNatom_residue  = atoi(value); }
    else { printf("#ERROR:Can't understand the option '%s'.\n",key); exit(1);}
  }

  get_Rmat_and_Tvec_from_strings(Rtar, Ttar, Rtar_str, Ttar_str);
  get_Rmat_and_Tvec_from_strings(Rgtar,Tgtar,Rgtar_str,Tgtar_str);

  setup_filenames(idref,typeref,assembly_id_ref,auth_asym_id_ref,iciffile_ref,ipdbfile_ref,igmmfile_ref);
  setup_filenames(idtar,typetar,assembly_id_tar,auth_asym_id_tar,iciffile_tar,ipdbfile_tar,igmmfile_tar);


  /*
  printf("#REFERENCE id %s type '%s' assembly_id '%s' auth_asym_id '%s' iciffile '%s' ipdbfile '%s' igmmfile '%s'\n",
                       idref,typeref,assembly_id_ref,auth_asym_id_ref,iciffile_ref,ipdbfile_ref,igmmfile_ref);
  printf("#TARGET   id %s type '%s' assembly_id '%s' auth_asym_id '%s' iciffile '%s' ipdbfile '%s' igmmfile '%s'\n",
                       idtar,typetar,assembly_id_tar,auth_asym_id_tar,iciffile_tar,ipdbfile_tar,igmmfile_tar);
   */

  /*** write headers ***/
  if (strcmp(OutType,"html")==0){
    printf("Content-type: text/html;charset=utf-8\n\n");
    printf("<HTML><BODY><PRE>\n");
  }
  else if (strcmp(OutType,"download")==0){
    printf("Content-type: text/download\n");
    download_file[0] = '\0';
    if ((idref[0]=='\0')&&(idtar[0]!='\0')&&(strcmp(typetar,"pdb")==0)){
        strcat(download_file,idtar);
        if (assembly_id_tar[0] != '\0'){
          strcat(download_file,"-");
          strcat(download_file,assembly_id_tar);
        }
       strcat(download_file,"_fit.pdb");
    }
    else{
      sprintf(download_file,"download.pdb");
    }
    printf("Content-Disposition: attachment; filename=%s\n\n",download_file);
  }
  else{
    printf("Content-type: text/html;charset=utf-8\n\n");
  }


  /** [ref-GMM] read and output GMM file for the reference **/ 
  if (OutGMM=='T'){
    read_pdb_and_out_transformed_xyz(igmmfile_ref,"-",' ','F','F',Rref,Tref,atof(tfref_str));
    fflush(stdout);
  }

  /** [ref-mmCIF] read and output mmCIF file for the reference **/ 
  if ((iciffile_ref[0] != '\0')&& ((strcmp(typeref,"pdb")==0) || (strcmp(typeref,"sasbdb")==0))){
    residue_select = 'F';
    read_mmCIF_file_into_BLOCK(iciffile_ref,&BLKref);
    make_UNITMOL_from_BLOCK(&PDBref,&BLKref,'T');
    make_ASSEMBLY_from_BLOCK(&PDBref,&BLKref);
    make_OPERATION_from_BLOCK(&PDBref,&BLKref);
 
   if ((assembly_id_ref[0]!='-')&&(assembly_id_ref[0]!='\0')){
      make_ASMBLMOL_from_BLOCK(&PDBref,&BLKref,assembly_id_ref);
      A = find_assembly(&PDBref,assembly_id_ref);
      if (A->Natom < maxNatom_residue){
        if (A->Natom > maxNatom_fullatom) {residue_select = 'T';} else {residue_select = 'F';}
        write_ASSEMBLY_in_pdb("-",&PDBref,&BLKref,assembly_id_ref,residue_select,'F',Rtar,Ttar,atof(tfref_str));
      }
    }
    else if ((assembly_id_ref[0]=='-')||(assembly_id_ref[0]=='\0')){
      if (PDBref.Natom < maxNatom_residue){
        if (PDBref.Natom > maxNatom_fullatom) {residue_select = 'T';} else {residue_select = 'F';}
        write_UNITMOLs_in_pdb("-",&PDBref,&BLKref,residue_select,auth_asym_id_ref,'F',Rtar,Ttar,atof(tfref_str));
      }
    }
  }
  fflush(stdout);
  
  /** [ref-PDB] read and output PDB file for the reference **/ 
  if (ipdbfile_ref[0] != '\0'){
    read_pdb_and_out_transformed_xyz(ipdbfile_ref,"-",' ','F','F',Rref,Tref,atof(tfref_str));
    fflush(stdout);
  }

  /** [tar-GMM] read and output GMM file for the target **/ 
  if (OutGMM=='T'){
    read_pdb_and_out_transformed_xyz(igmmfile_tar,"-",' ','F','T',Rgtar,Tgtar,atof(tftar_str));
    fflush(stdout);
  }

  /** [tar-mmCIF] read and output mmCIF file for the target **/ 
  if ((iciffile_tar[0] != '\0')&& ((strcmp(typetar,"pdb")==0) || (strcmp(typetar,"sasbdb")==0))){
    residue_select = 'F';

    read_mmCIF_file_into_BLOCK(iciffile_tar,&BLKtar);

    make_UNITMOL_from_BLOCK(&PDBtar,&BLKtar,'T');
    make_ASSEMBLY_from_BLOCK(&PDBtar,&BLKtar);
    make_OPERATION_from_BLOCK(&PDBtar,&BLKtar);
 
   if ((assembly_id_tar[0]!='-')&&(assembly_id_tar[0]!='\0')){
      make_ASMBLMOL_from_BLOCK(&PDBtar,&BLKtar,assembly_id_tar);
      A = find_assembly(&PDBtar,assembly_id_tar);
      if (A->Natom < maxNatom_residue){
        if (A->Natom > maxNatom_fullatom) {residue_select = 'T';} else {residue_select = 'F';}
        write_ASSEMBLY_in_pdb("-",&PDBtar,&BLKtar,assembly_id_tar,residue_select,'T',Rtar,Ttar,atof(tftar_str));
      }
    }
    else if ((assembly_id_tar[0]=='-')||(assembly_id_tar[0]=='\0')){
      if (PDBtar.Natom < maxNatom_residue){
        if (PDBtar.Natom > maxNatom_fullatom) {residue_select = 'T';} else {residue_select = 'F';}
        write_UNITMOLs_in_pdb("-",&PDBtar,&BLKtar,residue_select,auth_asym_id_tar,'T',Rtar,Ttar,atof(tftar_str));
      }
    }
  }

  fflush(stdout);
  
  /** [tar-PDB] read and output PDB file for the target **/ 
  if (ipdbfile_tar[0] != '\0'){
    read_pdb_and_out_transformed_xyz(ipdbfile_tar,"-",' ','F','T',Rtar,Ttar,atof(tftar_str));
    fflush(stdout);
  }

  if (strcmp("out","html")==0){
    printf("</PRE></BODY></HTML>\n");
  }

  return(1);
}
