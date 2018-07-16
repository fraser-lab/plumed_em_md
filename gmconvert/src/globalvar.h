/*

<globalvar.h>

*/

#define LAST_MOD_DATE   "Feb 26, 2016"
/* #define MAX_MEMORY_MEGA 1024 */
#define MAX_MEMORY_MEGA 2500 
#define MAX_FILENAME    512


/** PARAMETERS: for globally defined important parameters **/
struct PARAMETERS{
  char   COMMAND[MAX_FILENAME];
  char   HETtype;   /* Read HETATM ('T' or 'F') */
  char   MODEL;     /* 'S':read only 1st model in the file, 'M':read multiple model  */

  /* Logfile */
  char ologfile[MAX_FILENAME];
  FILE *fp_olog;

  char OutStdLog;  /* output convergence log as stdout ('T' or 'F') */
  /** Table for exp(-x) **/
  char TabExpMinX;  /* Use table for exp(-x) : 'T' or 'F' */
 
  /* Others */
  char   START_DATE[MAX_FILENAME];  /* Date for starting calculation */
  char   END_DATE[MAX_FILENAME];    /* Date for end      calculation */
  double START_TIME_SEC;            /* Date on starting time ( in second) */
  double COMP_TIME_SEC;             /* Computational time (in second); */
  char   SubType[128];              /* Type for Various Purpuse */
  float  eps;                       /* Very small number for numerical calculation (Angstrom) */
  char   Vtype;                     /* Visible for progress of calculation */ 

};


/** GLOBAL VARIABLES **/
extern struct PARAMETERS PAR;
extern struct TAB_EXP_MIN_X TABEXP;
