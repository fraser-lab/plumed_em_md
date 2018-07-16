/*

 <GridMap.h>

*/

struct CHAR2DMAP{
 int  N[2];           /* Number of grids */
 unsigned char **dat; /* m[N[0]][N[1]] (malloc later) */
};

struct CHAR3DMAP{
 int  N[3];            /* Number of grids */
 unsigned char ***dat; /* m[N[0]][N[1]][N[2]] (malloc later) */
};

struct CHAR4DMAP{
 int  N[4];             /* Number of grids */
 unsigned char ****dat; /* m[N[0]][N[1]][N[2]][N[3]] (malloc later) */
};

struct FLOAT1DMAP{
 int    N;         /* Number of grids */
 float *dat;    /* m[N] (malloc later) */
};

struct FLOAT2DMAP{
 int  N[2];       /* Number of grids */
 float **dat;    /* m[N[0]][N[1]] (malloc later) */
};

struct FLOAT3DMAP{
 int  N[3];       /* Number of grids */
 float ***dat;     /* m[N[0]][N[1]][N[2]] (malloc later) */
};

struct FLOAT4DMAP{
 int  N[4];        /* Number of grids */
 float ****dat;    /* m[N[0]][N[1]][N[2]][N[3]] (malloc later) */
};

struct DOUBLE1DMAP{
 int  N;         /* Number of grids */
 double *dat;    /* m[N] (malloc later) */
};

struct DOUBLE2DMAP{
 int  N[2];       /* Number of grids */
 double **dat;    /* m[N[0]][N[1]] (malloc later) */
};

struct DOUBLE3DMAP{
 int  N[3];       /* Number of grids */
 double ***dat;     /* m[N[0]][N[1]][N[2]] (malloc later) */
};

struct DOUBLE4DMAP{
 int  N[4];        /* Number of grids */
 double ****dat;    /* m[N[0]][N[1]][N[2]][N[3]] (malloc later) */
};

/*** FUNCIONS (GLOBAL) ***/
extern int Malloc_CHAR2DMAP();
extern void Initialize_CHAR2DMAP();
extern void Free_CHAR2DMAP();

extern int Malloc_CHAR3DMAP();
extern void Initialize_CHAR3DMAP();
extern void Free_CHAR3DMAP();

extern int Malloc_CHAR4DMAP();
extern void Initialize_CHAR4DMAP();
extern void Free_CHAR4DMAP();

extern int Malloc_FLOAT1DMAP();
extern void Initialize_FLOAT1DMAP();
extern void Free_FLOAT1DMAP();

extern int Malloc_FLOAT2DMAP();
extern void Initialize_FLOAT2DMAP();
extern void Free_FLOAT2DMAP();

extern int Malloc_FLOAT3DMAP();
extern void Initialize_FLOAT3DMAP();
extern void Free_FLOAT3DMAP();

extern int Malloc_FLOAT4DMAP();
extern void Initialize_FLOAT4DMAP();
extern void Free_FLOAT4DMAP();

extern int Malloc_DOUBLE1DMAP();
extern void Initialize_DOUBLE1DMAP();
extern void Free_DOUBLE1DMAP();

extern int Malloc_DOUBLE2DMAP();
extern void Initialize_DOUBLE2DMAP();
extern void Free_DOUBLE2DMAP();

extern int Malloc_DOUBLE3DMAP();
extern void Initialize_DOUBLE3DMAP();
extern void Free_DOUBLE3DMAP();

extern int Malloc_DOUBLE4DMAP();
extern void Initialize_DOUBLE4DMAP();
extern void Free_DOUBLE4DMAP();

