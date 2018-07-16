/*

<mc_verface.h>

 Definition of Vertex and Face for Marching Cube,
 and MIDPNT and MIDMAP for overlap checking map
 
  When we want to know a midpoint between [x0][y0][z0] and [x1][y1][z1] is
  already constructed or not, 
  we check MIDMAP[x0][y0][z0]->ver[x1-x0+1][y1-y0+1][z1-z0+1] is NULL or not.

*/

#define MAX_MC_VERTEX_NEIGHBOR 5

struct MC_VERTEX{
 int num;                 /* Number of Vertex [0..Nvertex-1]*/
 float pos[3];            /* XYZ coordinates */
  /* 
      pos[3] is described in the "lattice" unit. 
      The "real" xyz coordinate real_pos is calculated as follows:
         real_pos[0] = vox->OrigPos[0] + vox->grid_width * mc_vertex->pos[0]
         real_pos[1] = vox->OrigPos[1] + vox->grid_width * mc_vertex->pos[1]
         real_pos[2] = vox->OrigPos[2] + vox->grid_width * mc_vertex->pos[2]
  */
 struct MC_VERTEX *next;       /* for linked list  */
 struct MC_VERTEX *neighbors[MAX_MC_VERTEX_NEIGHBOR];   /* this is for checking existing MC_EDGE or not */
};


struct MC_FACE{
 int num;                   /* Number of Face [0..Nface-1]*/
 struct MC_VERTEX *a,*b,*c; /* Three MC vertices which composed of triangle */
 struct MC_FACE *next;      /* for linked list  */
};

struct MC_EDGE{
 struct MC_VERTEX *a,*b; /* Two MC vertices which composed of edge */
 struct MC_EDGE *next;   /* for linked list  */
};


/*
 
 MIDPNT and MIDMAP_YZ are for checking overlap of midpoints quickly. 

*/


struct MIDPNT{
 /*
  Pointer to the midpoint MC vertex.
   dx = (-1,0,+1) -> [0,1,2] 
   dy = (-1,0,+1) -> [0,1,2] 
   dz = (-1,0,+1) -> [0,1,2] 

  For the point (x,y,z),
    midpoint bwn [x,y,z] and [x-1,y-1,z-1] : ver[0][0][0] 
    midpoint bwn [x,y,z] and [x+1,y+1,z+1] : ver[2][2][2] 
    midpoint bwn [x,y,z] and [x-1,y  ,z+1] : ver[0][1][2] 
 */
 struct MC_VERTEX *ver[3][3][3];
}; 



struct MIDMAP_YZ{
 int Y,Z; 
 struct MIDPNT **map; /* map[Y][Z] */
};
