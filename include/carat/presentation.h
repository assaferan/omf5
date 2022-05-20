#ifndef CARAT_PRESENTATION_H
#define CARAT_PRESENTATION_H

/**********************************************************************
|  FILE: presentation.c
***********************************************************************/
extern matrix_TYP *pres(bahn **s,
                        bravais_TYP *G,
                        int *OPT);

extern void normalize_word(int *w);

/**********************************************************************
|  FILE: mapped_word.c
***********************************************************************/
matrix_TYP *mapped_word(int *w,
                        matrix_TYP **A,
                        matrix_TYP **AINV);

/**********************************************************************
|  FILE: check_base.c
***********************************************************************/
void check_base(bahn **s,
                bravais_TYP *G);

/**********************************************************************
|  FILE: put_word.c
***********************************************************************/
void put_word(int *w,
              const char *O);

#endif
