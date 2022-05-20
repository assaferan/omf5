#ifndef CARAT_POLYEDER_H
#define CARAT_POLYEDER_H

#include "typedef.h"


/*-------------------------------------------------------------*\
| FILE: first_polyeder.c 
\*-------------------------------------------------------------*/
extern polyeder_TYP *first_polyeder(wall_TYP **mauern, int anz);

/*-------------------------------------------------------------*\
| FILE: polyeder_tools.c 
\*-------------------------------------------------------------*/
extern vertex_TYP *init_vertex(int dim, int wall_no) ;
extern word_TYP *init_word(int dim);
extern wall_TYP *init_wall(int dim);
extern polyeder_TYP *init_polyeder(int vert_no, int wall_no);
extern polyeder_TYP *get_polyeder(const char *file_name);
extern void put_polyeder(polyeder_TYP *F);
extern int wall_times_vertex(wall_TYP *w, vertex_TYP *v);
extern void free_vertex(vertex_TYP **v);
extern void free_word(word_TYP *word);
extern void free_wall(wall_TYP **v);
extern wall_TYP *mat_to_wall(matrix_TYP *M);
extern void normal_wall(wall_TYP *v);
extern void normal_vertex(vertex_TYP *v);
extern int is_vertex_of_wallno(vertex_TYP *v, int w);
extern word_TYP *copy_word(word_TYP *w);
extern wall_TYP *copy_wall(wall_TYP *w);
extern void free_polyeder(polyeder_TYP *P);

/*-------------------------------------------------------------*\
| FILE: polyeder_to_vecs.c 
\*-------------------------------------------------------------*/
extern matrix_TYP **polyeder_to_vecs(polyeder_TYP *P);

/*-------------------------------------------------------------*\
| FILE: refine_polyeder.c 
\*-------------------------------------------------------------*/
extern int refine_polyeder(polyeder_TYP *F, wall_TYP *h);

/*-------------------------------------------------------------*\
| FILE: fuber_tools.c 
\*-------------------------------------------------------------*/
extern vertex_TYP *init_vertex_fuber (int dim, int wall_no);
extern wall_TYP *init_wall_fuber (int dim);
extern fund_domain *init_fund_domain(int vert_no, int wall_no);
extern void free_vertex_fuber(vertex_TYP**);
extern void free_wall_fuber(wall_TYP**);
extern int wall_times_vertex_fuber(wall_TYP*, vertex_TYP*);

#endif
