#ifndef CARAT_GETPUT_H
#define CARAT_GETPUT_H

#include "typedef.h"

/*-------------------------------------------------------------*\
| FILE: get_bravais.c 
\*-------------------------------------------------------------*/
extern bravais_TYP *get_bravais (const char *file_name);

/*-------------------------------------------------------------*\
| FILE: get_mat.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *fget_mat (FILE *infile);
extern matrix_TYP *get_mat (const char *file_name);
extern matrix_TYP **mget_mat (const char *file_name, int *anz);
extern matrix_TYP **fmget_mat (FILE *infile, int *anz);

/*-------------------------------------------------------------*\
| FILE: put_bravais.c 
\*-------------------------------------------------------------*/
extern void fput_bravais(FILE *outfile, bravais_TYP *G, const char *comment);
extern void put_bravais(bravais_TYP *G, const char *filename, const char *comment);

/*-------------------------------------------------------------*\
| FILE: put_mat.c 
\*-------------------------------------------------------------*/
#define PM_RATIONAL_BIT 0
#define PM_SHORTCUT_BIT 1
#define PM_RATIONAL (1UL << PM_RATIONAL_BIT)
#define PM_SHORTCUT (1UL << PM_SHORTCUT_BIT)

extern void put_mat ( matrix_TYP *mat, const char file_name[], const char comment[],
     unsigned long options);
extern void fput_mat (FILE *outfile, matrix_TYP *mat, const char comment[],
     unsigned long options);

/*-------------------------------------------------------------*\
| FILE: put_order.c 
\*-------------------------------------------------------------*/
extern void fput_order( FILE *outfile,  int *divisors,  int ord);

/*-------------------------------------------------------------*\
| FILE: read_header.c 
\*-------------------------------------------------------------*/
extern char ** FILENAMES;
extern int FILEANZ;

extern void read_header( int argc, char *argv[]);
extern int is_option( char c);
extern int optionnumber(char c);

#endif
