#ifndef CARAT_TYPEDEF_H
#define CARAT_TYPEDEF_H

/* enthaelt die globalen Variablen - und Typendeklarationen */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include <ctype.h>
#include <signal.h>
#include <stdarg.h>

#define TRUE 1
#define FALSE 0

#define EXT_SIZE 100
#define extsize1 30
#define MAXDIM 6

typedef int boolean;
typedef struct { int z; int n; } rational;
typedef struct { int ggt, f1, f2; int alt, g1, g2; } pair;

typedef struct {
  int **SZ;
  int **N;
} array_TYP;

typedef struct {
   int     Integral ,
           Symmetric,
           Diagonal ,
           Scalar    ;
} flag_TYP;

typedef struct {
   flag_TYP flags;
   int cols, rows;
   int kgv;
   int prime;
   array_TYP array;
} matrix_TYP;


typedef struct {
   int dim;
   int gen_no;
   int order;
   int divisors[100];
   int form_no;
   int zentr_no;
   int normal_no;
   int cen_no;
   matrix_TYP **gen;
   matrix_TYP **form;
   matrix_TYP **zentr;
   matrix_TYP **normal;
   matrix_TYP **cen;
} bravais_TYP;

typedef struct {
  bravais_TYP *grp;
  char *fn;
} symbol_out;

typedef struct {
  int *v;
  int dim;
  int *wall;
  int wall_no;
  int wall_SIZE;
} vertex_TYP;

/* inserted tilman 16/07/97 to add more functionality for anne's programs */
typedef struct{
  int dim;
  int *word;
  matrix_TYP *trans;
} word_TYP;

typedef struct{
  int *gl;
  int dim;
  int *product;
  int nproduct;
  int norm;          /* next 4 lines inserted by anne, 8/10/97 */
  int next_no;       /* Anzahl der Nachbarwaende */
  int **next;        /* Gleichungen der Nachbarwaende */
  int ext_no;        /* Anzahl der virtuellen Waende */
  int **extra;       /* Gleichungen der virtuellen Waende */
  int neu;           /* neu[i] = 0 falls die Wand neu hinzugekommen ist */
  int paar;          /* No. der zu dieser Wand gepaarten Wand */
  matrix_TYP *mat;   /* Seitentransformation */
  word_TYP *word;    /* Wort der Seitentrafo in den Gruppenerzeugern */
}wall_TYP;

			/* next 5 lines inserted by anne, 8/10/98 */
typedef struct{
	int w[2];
	int *v;
	int v_anz;
	} corner_TYP;
 
typedef struct {
   vertex_TYP **vert;
   int vert_no;
   int vert_SIZE;
   wall_TYP **wall;
   int wall_no;
   int wall_SIZE;
   corner_TYP *corner;  /* next 3 lines inserted by anne, 8/10/98 */
   int corner_no;
   int corner_SIZE;
   int is_closed;
   int is_degenerate;
} polyeder_TYP;

/* the programs of joerg kock will need this setting */
typedef struct {
   vertex_TYP **vert;
   int vert_no;
   int vert_SIZE;
   wall_TYP **wall;
   int wall_no;
   int wall_SIZE;
   int is_closed;
   int is_finite;
} fund_domain;


struct baum{
   int no;
   struct baum *left;
   struct baum *right;
};


struct tree{
   int no;
   struct tree *left;
   struct tree *right;
};

/* the setting for programs which sit in ...../functions/Base */
#define MIN_SPEICHER 256

typedef struct {
        int length;
        int speicher;
        matrix_TYP **orbit;
        matrix_TYP **representatives;
        matrix_TYP **rep_invs;
        matrix_TYP **generators;
        int **words;
        int gen_no;
        struct tree *hash;
} bahn;

typedef struct {
  int dim;
  int *v;
  int kgv;
} vector_TYP;

/* for QtoZ */
typedef struct {
  int *s;
  matrix_TYP ****Delta;
  int k;
  int r;
} QtoZ_konst_TYP;

typedef struct {
   int anz;
   int *I;
   int *J;
   int *flag;
   matrix_TYP **lattice;
   matrix_TYP **lsf;		/* standard form of the lattices */
} QtoZ_entry_TYP;

typedef struct {
   matrix_TYP **gitter;
   matrix_TYP **tr_gitter;
   matrix_TYP **inv_tr_gitter;
   int anz;
   QtoZ_entry_TYP **entry;
   matrix_TYP ***zoogitter;
} QtoZ_TYP;

// TODO: move the following to a more appropriate location
extern int INFO_LEVEL;

#endif /* CARAT_TYPEDEF_H */
