#ifndef CARAT_DATEI_H
#define CARAT_DATEI_H

#include "typedef.h"

typedef struct{
    bravais_TYP *grp;
    char *symbol;
    int almost;
    int zclass;
    int alpha;
    int N_orbits;
    matrix_TYP **TR;
} lattice_element;

/*-------------------------------------------------------------*\
| FILE: brav_from_datei.c
\*-------------------------------------------------------------*/
bravais_TYP *brav_from_datei(const char *symb,int almost,int zclass);

/*-------------------------------------------------------------*\
| FILE: free_bravais.c 
\*-------------------------------------------------------------*/
extern void free_bravais( bravais_TYP *grp);

/*-------------------------------------------------------------*\
| FILE: get_symbol.c 
\*-------------------------------------------------------------*/
extern symbol_out *get_symbol ( const char *file_name);

/*-------------------------------------------------------------*\
| FILE: get_zentr.c 
\*-------------------------------------------------------------*/
extern void get_zentr( symbol_out *B);

/*-------------------------------------------------------------*\
| FILE: gittstab.c 
\*-------------------------------------------------------------*/
extern bravais_TYP *gittstab( bravais_TYP *grp, matrix_TYP *X);

extern bravais_TYP *Z_class( bravais_TYP *B, matrix_TYP *zen);

/*-------------------------------------------------------------*\
| FILE: gittstabneu.c
\*-------------------------------------------------------------*/
extern bravais_TYP *gittstabneu( bravais_TYP *grp, matrix_TYP *X);

/*-------------------------------------------------------------*\
| FILE: read_symbol.c 
\*-------------------------------------------------------------*/
extern symbol_out *read_symbol(const char *file_name);

/*-------------------------------------------------------------*\
| FILE: right_order.c
\*-------------------------------------------------------------*/
extern void right_order(char *string);

/*-------------------------------------------------------------*\
| FILE: lattice_tools.c
\*-------------------------------------------------------------*/

lattice_element *init_lattice_element();
void free_lattice_element(lattice_element *x);
lattice_element *fget_lattice_element(FILE *F,int OPTION);
void fput_lattice_element(lattice_element *E,FILE *F);

/*-------------------------------------------------------------*\
| FILE: lattice.c
\*-------------------------------------------------------------*/

lattice_element **lattice(const char *symb,int dim,int almost,int zclass,int *no,
                          int OPTION);

/*------------------------------------------------------------*\
| FILE: super_lattice.c
\*------------------------------------------------------------*/

lattice_element **super_lattice(const char *symb,int dim,int almost,int zclass,
		int *no, int OPTION);

/*------------------------------------------------------------*\
| FILE: get_data_dir.c
\*------------------------------------------------------------*/

void setup_carat_location(const char * argv0);
const char *get_data_dir(void);

#endif
