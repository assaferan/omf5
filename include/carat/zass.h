#ifndef CARAT_ZASSEN_H
#define CARAT_ZASSEN_H

#include "typedef.h"

typedef struct {
       long last;
       long speicher;
       long faktor;
       long *pointer;
} word;


/**********************************************************************\
| FILE: cobundary.c
\**********************************************************************/

void coboundary(bravais_TYP *G,
                matrix_TYP *C,
                matrix_TYP *T);

/**********************************************************************\
| FILE: cong_solve.c
\**********************************************************************/

matrix_TYP **cong_solve(matrix_TYP *A);

/**********************************************************************\
| FILE: convert_cocycle_to_column.c
\**********************************************************************/

void convert_cocycle_to_column(matrix_TYP **Y,
                               int number,
                               int dim,
                               int gen_no);

/**********************************************************************\
| FILE: convert_to_cocycle.c
\**********************************************************************/

extern matrix_TYP *convert_to_cocycle(matrix_TYP *x,
                                      matrix_TYP *cocycle,
                                      matrix_TYP *D);

/**********************************************************************\
| FILE: put_cocycle.c
\**********************************************************************/

void put_cocycle(matrix_TYP *COZ,
                 int dim,
                 int number,
                 const char *file,
                 const char *comment);

/*********************************************************************\
| FILE: normalop.c
\*********************************************************************/

extern matrix_TYP *normalop(matrix_TYP *cocycle,
                            matrix_TYP *D,
                            matrix_TYP *R,
                            bravais_TYP *G,
                            matrix_TYP *N,
                            int opt);

extern matrix_TYP **extensions(matrix_TYP *cocycle,
                               matrix_TYP *D,
                               matrix_TYP *R,
                               bravais_TYP *G,
                               int **lengths,
                               MP_INT **names,
                               int *number_of_orbits,
                               int option);

extern void no_of_extensions(matrix_TYP *cocycle,matrix_TYP *D,
                        matrix_TYP *R,bravais_TYP *G,MP_INT *no);

extern matrix_TYP **identify(matrix_TYP *cocycle,
                             matrix_TYP *D,
                             matrix_TYP *R,
                             bravais_TYP *G,
                             matrix_TYP **extension,
                             MP_INT *a,int number,
                             int transform_flag,
                             int ***WORDS,
                             int *NUMBER_OF_WORDS);

extern void translation(matrix_TYP *TR,
                        matrix_TYP *rep,
                        matrix_TYP *ext,
                        matrix_TYP *cocycle,
                        matrix_TYP *D,
                        bravais_TYP *G);

extern void valuation(matrix_TYP *x,matrix_TYP *D,MP_INT *val);

extern matrix_TYP *orbit_rep(matrix_TYP *x,
                             matrix_TYP **N,
                             int nanz,
                             matrix_TYP *D,
                             int option,
                             char *B,
                             MP_INT *l,
                             int *anz,
                             int **word,
                             int word_flag,
                             int ***WORDS,
                             int *NUMBER_OF_WORDS);

/**********************************************************************\
| FILE: reverse_valuation.c
\**********************************************************************/

extern matrix_TYP *reverse_valuation(MP_INT *val,matrix_TYP *D);

/**********************************************************************\
| FILE: reget_gen.c
\**********************************************************************/

extern matrix_TYP *reget_gen(matrix_TYP **map,int number,bravais_TYP *G,
                             int **words, int word_flag);


/**********************************************************************\
| FILE: zass.c
\**********************************************************************/

extern matrix_TYP *matrizen_in_word(matrix_TYP **mat,
                                    matrix_TYP **matinv,
                                    word g);

extern void matrix_2_word(matrix_TYP *matrix,
                          word *relator,
                          long zeile);

extern matrix_TYP **cohomology(long *dim,
                               matrix_TYP **mat,
                               matrix_TYP **matinv,
                               word *relator,
                               int erzeuger,
                               int relatoren);

extern int wordfree(word *a);

extern matrix_TYP *scalar(long n,
                          long a);

extern matrix_TYP *calc_B(matrix_TYP **mat,
                          long anz_erzeuger);

#endif
