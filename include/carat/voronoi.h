#ifndef CARAT_VORONOI_H
#define CARAT_VORONOI_H

#include "typedef.h"

typedef struct {
   matrix_TYP *gram;
   int SV_no;
   int min;
   int pdet;
   int prime;
   int vert_no;
   wall_TYP **vert;
   polyeder_TYP *pol;
   matrix_TYP *red_inv;
   matrix_TYP *T;
   matrix_TYP *Ti;
   bravais_TYP *Gtrred;
   matrix_TYP *SVi;
   bravais_TYP *stab;
   bravais_TYP *linstab;
   matrix_TYP *dir_reps;
} voronoi_TYP;


/*********************************************************************\
| FILE:  all_vor_neighbours.c
\*********************************************************************/
extern matrix_TYP *all_voronoi_neighbours(matrix_TYP *P, bravais_TYP *G,
     matrix_TYP **Ftr, matrix_TYP *tr_bifo);

/*********************************************************************\
| FILE:  calc_vor_data.c
\*********************************************************************/
extern void calc_voronoi_basics(voronoi_TYP *V, bravais_TYP *G,
     bravais_TYP *Gtr, int prime);
extern void calc_voronoi_pol(voronoi_TYP *V, matrix_TYP *bifo);
extern void calc_voronoi_good_inv(voronoi_TYP *V, bravais_TYP *Gtr);
extern void calc_voronoi_stab(voronoi_TYP *V, bravais_TYP *G, bravais_TYP *Gtr,
     matrix_TYP *bifo);
extern matrix_TYP *calc_voronoi_isometry(voronoi_TYP *V1, voronoi_TYP *V2,
     bravais_TYP *G, bravais_TYP *Gtr, matrix_TYP *bifo);
extern void calc_voronoi_dir_reps(voronoi_TYP *V, bravais_TYP *G,
     bravais_TYP *Gtr, matrix_TYP *bifo);
extern void calc_voronoi_complete(voronoi_TYP *V, bravais_TYP *G,
     bravais_TYP *Gtr, matrix_TYP *bifo, int prime);

/*********************************************************************\
@  FILE: first_perfect.c
\*********************************************************************/
extern matrix_TYP *first_perfect(matrix_TYP *A, bravais_TYP *G,
     matrix_TYP **Ftr, matrix_TYP *trbifo, int *min);

/*********************************************************************\
@  FILE: init_voronoi.c
\*********************************************************************/
extern voronoi_TYP *init_voronoi();
extern void clear_voronoi(voronoi_TYP *V);
extern void put_voronoi(voronoi_TYP *V);

/****************************************************************************\
@  FILE: normalizer.c
\****************************************************************************/
extern voronoi_TYP **normalizer(matrix_TYP *P, bravais_TYP *G,
     bravais_TYP *Gtr, int prime, int *V_no);

/************************************************************************\
|  FILE: pair_red_inv.c
\************************************************************************/
extern matrix_TYP *pair_red_inv(matrix_TYP *A, matrix_TYP *T);

/****************************************************************************\
@  FILE: vor_neighbour.c:
\****************************************************************************/
extern matrix_TYP *voronoi_neighbour(matrix_TYP *A, matrix_TYP *X, 
     int Amin, int *lc, int *rc);

/****************************************************************************\
@  FILE: vor_vertices.c
\****************************************************************************/
extern matrix_TYP **voronoi_vertices(matrix_TYP *form, bravais_TYP *grp,
   int *anz, int *form_min, int *SV_no);

/****************************************************************************\
@  FILE: bravais_flok.c
\****************************************************************************/
extern matrix_TYP *is_z_equivalent(bravais_TYP *G,bravais_TYP *Gtr,
                                   bravais_TYP *H,bravais_TYP *Htr);

matrix_TYP *extends_to_isometry(
matrix_TYP **hforms,matrix_TYP *HSV,int anz_hneighbours,
matrix_TYP **gforms,matrix_TYP *GSV,int anz_gneighbours,
int fdim,int offset);

void transform_pair(bravais_TYP *H,bravais_TYP *Htr,matrix_TYP *x);

int max_diagonal_entry(matrix_TYP *A);

int neighbours(matrix_TYP ***perf,bravais_TYP *G,matrix_TYP **Ftr,
               matrix_TYP *tr_bifo,matrix_TYP *SV,int min);

/****************************************************************************\
@  FILE: bravais_flok_datei.c
\****************************************************************************/
extern matrix_TYP *is_z_equivalent_datei(bravais_TYP *G,bravais_TYP *Gtr,
           bravais_TYP *H,bravais_TYP *Htr,voronoi_TYP ***gp,int *anz_gperfect);

/****************************************************************************\
@  FILE: red_normal.c
\****************************************************************************/
extern void red_normal(bravais_TYP *G);

#endif
