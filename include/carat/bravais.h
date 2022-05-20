#ifndef CARAT_BRAVAIS_H
#define CARAT_BRAVAIS_H

#include "typedef.h"

/*-------------------------------------------------------------*\
| FILE: bravais_tools.c
\*-------------------------------------------------------------*/
bravais_TYP *bravais_group(bravais_TYP *H,int FLAG);
bravais_TYP *copy_bravais(bravais_TYP *H);

/*-------------------------------------------------------------*\
| FILE: formspace.c 
\*-------------------------------------------------------------*/
extern matrix_TYP **formspace(matrix_TYP **B, int Banz, int sym_opt, int *fdim);

/*-------------------------------------------------------------*\
| FILE: init_bravais.c 
\*-------------------------------------------------------------*/
extern bravais_TYP *init_bravais(int dim);

/*-------------------------------------------------------------*\
| FILE: invar_space.c 
\*-------------------------------------------------------------*/
extern matrix_TYP **invar_space(matrix_TYP **B, int Banz, int fodim,
     int symm_opt, int epsilon, int *anz);

/*-------------------------------------------------------------*\
| FILE: konj_bravais.c 
\*-------------------------------------------------------------*/
extern bravais_TYP *konj_bravais(bravais_TYP *B, matrix_TYP *T);

/*-------------------------------------------------------------*\
| FILE: lincomb.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *vec_to_form(int *v, matrix_TYP **F, int Fanz);
extern void form_to_vec(int *erg, matrix_TYP *A, matrix_TYP **F,
     int Fanz, int *denominator);
extern vertex_TYP *form_to_vertex(matrix_TYP *A, matrix_TYP **F,
     int Fanz, int *denominator);
extern void form_to_vec_modular(int *erg, matrix_TYP *A, matrix_TYP **F,
     int Fanz);

/*-------------------------------------------------------------*\
| FILE: normlin.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *normlin(matrix_TYP **Fo, matrix_TYP *N, int fdim);

/*-------------------------------------------------------------*\
| FILE: p_formspace.c 
\*-------------------------------------------------------------*/
extern matrix_TYP **p_formspace(matrix_TYP **B, int Banz, int prime,
      int sym_opt, int *fdim);

/*-------------------------------------------------------------*\
| FILE: rform.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *rform(matrix_TYP **B, int Banz, matrix_TYP *Fo, int epsilon);

/*-------------------------------------------------------------*\
| FILE: trace_bifo.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *trace_bifo(matrix_TYP **F1, matrix_TYP **F2, int anz);

/*-------------------------------------------------------------*\
| FILE: tr_bravais.c 
\*-------------------------------------------------------------*/
extern bravais_TYP *tr_bravais(bravais_TYP *B, int calcforms, int invert);

/*-------------------------------------------------------------*\
| FILE: normalisator.c
\*-------------------------------------------------------------*/
void normalisator(bravais_TYP *H,
                  bravais_TYP *Gtr,
                  matrix_TYP *A,
                  int prime,
                  boolean b_option,
                  boolean o_option);

#endif
