/*

<BasicIO.c>

 Baisc functions for input/output of files

*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <dirent.h>

#include "globalvar.h" 

/*** Functions (GLOBAL) ***/
void Find_Filename_Core();
void Get_Part_Of_Line();
char *Get_Date_String();
double Get_Time_in_Second_by_Double();
void Set_END_DATE();
void Split_to_Words();
void RGBTstr_to_floatRGBT();



void Get_Part_Of_Line(part,line,s,e)
  char *part;
  char *line;
  int  s,e;
{
 int i,E,L;
 L = strlen(line)-1;
 if (line[L] == '\n') L -= 1;
 if (e>L) E = L; else E = e;
 for (i=s;i<=E;++i) part[i-s] = line[i];
 part[E-s+1] = '\0';

} /* end of Get_Part_of_Line() */





char *Get_Date_String()
{
 time_t      now_t;
 struct tm  *loc_t;
 static char Mon[][4] = {"Jan","Feb","Mar","Apr","May","Jun",
                         "Jul","Aug","Sep","Oct","Nov","Dec"};
 static char string[64];
 now_t = time(NULL);
 loc_t = localtime(&now_t);

 sprintf(string,"%s %d,%d %d:%d:%d",
  Mon[loc_t->tm_mon],loc_t->tm_mday,loc_t->tm_year+1900,
  loc_t->tm_hour,loc_t->tm_min,loc_t->tm_sec);
 return(string);

} /* end of Get_Date_String() */


double Get_Time_in_Second_by_Double()
{
 struct timeval tv;
 gettimeofday(&tv,NULL);
 return(tv.tv_sec + (double)tv.tv_usec*1e-6);
} /* end of Get_Time_in_Second_by_Double() */



void Set_END_DATE(){
 double end_sec;
 sprintf(PAR.END_DATE,"%s",Get_Date_String());
 end_sec = Get_Time_in_Second_by_Double();
 PAR.COMP_TIME_SEC = end_sec - PAR.START_TIME_SEC;
}





void Find_Filename_Core(core,fname)
 char *core,*fname;
{
 int Wsta[100],Wend[100],Nword;
 char lastfile[128],hit;
 int i,dot_pos;

 Split_to_Words(fname,'/',&Nword,Wsta,Wend,100);
 Get_Part_Of_Line(lastfile,fname,Wsta[Nword-1],Wend[Nword-1]);
 i = strlen(lastfile)-1; 
 hit = 0;
 while ((i>0)&&(hit==0))
 {
  if ((lastfile[i] == '.') && (lastfile[i-1] != '.')) 
   {hit = 1; dot_pos = i-1;}
  --i; 
 }

 if (hit==1) Get_Part_Of_Line(core,lastfile,0,dot_pos);
 else sprintf(core,"%s",lastfile);

} /* end of Find_Filename_Core() */
 

void Split_to_Words(str,splsym,Nword,Wsta,Wend,Nwordmax)
 char *str;          /* Input String */
 char splsym;        /* Symbol for split */
 int *Nword;         /* Number of words  */
 int *Wsta;         /* Start point of str for a wowd */
 int *Wend;         /* End point of str for a wowd */
 int Nwordmax;       /* Maxium number of Nword  */
{
 /* [Example]
  str = "abc/d/ef//ghi//"
  splsym = '/'
   -->
  Nword 4
  (Wsta,Wend) = {(0,2), (4,4), (6,7), (10,12)}
 */
 int i,L;
/* 
 printf("#Split_to_Words(str '%s' splsym '%c' Nwordmax %d)\n",str,splsym,Nwordmax);
 */
 L = strlen(str);
 *Nword = 0; i = 0;

 while ((i<L)&&(*Nword < Nwordmax)){
   if (str[i]!=splsym){ 
     Wsta[*Nword] = i;
     while ((str[i]!=splsym)&&(i<=(L-1))) { ++i; }
     Wend[*Nword] = i-1;
     ++(*Nword);
   }
   ++i;
 }

} /* end of Split_to_Words() */



void RGBTstr_to_floatRGBT(RGBTstr,RGBT)
  char *RGBTstr;
  float RGBT[4];
{
  int Nword;
  int Wsta[10],Wend[10];
  char buff[100];
 
  Split_to_Words(RGBTstr,':',&Nword,Wsta,Wend,10);
  if (Nword==4){
    Get_Part_Of_Line(buff,RGBTstr,Wsta[0],Wend[0]); RGBT[0] = atof(buff); /* printf("'%s'\n",buff); */
    Get_Part_Of_Line(buff,RGBTstr,Wsta[1],Wend[1]); RGBT[1] = atof(buff); /* printf("'%s'\n",buff); */
    Get_Part_Of_Line(buff,RGBTstr,Wsta[2],Wend[2]); RGBT[2] = atof(buff); /* printf("'%s'\n",buff); */
    Get_Part_Of_Line(buff,RGBTstr,Wsta[3],Wend[3]); RGBT[3] = atof(buff); /* printf("'%s'\n",buff); */
  }
  else{
    printf("#ERROR: format error of RGBTstr ('%s').\n",RGBTstr); exit(1);
  }
  printf("#'%s' -> %lf %lf %lf %lf\n",RGBTstr,RGBT[0],RGBT[1],RGBT[2],RGBT[3]);
} /* end of RGBTstr_to_floatRGBT() */
