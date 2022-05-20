#ifndef CARAT_NAME_H
#define CARAT_NAME_H

#include "typedef.h"

#include "gmp.h"


#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

/* I tried to make the program as modular as possible.
   So I hope that it will be easy to expand the entries
   of the database. ...
   */

#define ALL_MATCH 65535

#define NR_OF_ELEMENTS_IN_EACH_ENTRY 10

#define COND_ABBREVIATION 0
#define COND_DEGREE 1
#define COND_SYMBOL 2
#define COND_ORDER 3
#define COND_DISCRIMINANT 4
#define COND_ZCLASSES 5
#define COND_AFFINE 6
#define COND_TORSIONFREE 7
#define COND_NO_CONCLASS 8
#define COND_NO_IDEM 9

#define SET_COND 0
#define DEL_COND 1
#define DISPLAY_POSSIBLE 2

typedef struct
{
  char *abbreviation;  /* filename of the group*/
  int degree;   /* dimension of the group */
  char *symbol;  /* family symbol */
  int order;
  char *discriminant;
  int zclasses;  /* number of Z-classes in the Q-class */
  int affine;  /* number of affine classes in the Q-class */
  int torsionfree; /* number of torsionfree affine classes */
  int no_conclass;  /*  number of conjugacy classes in group. */
  int no_idem;  /*  number of conjugacy classes in group. */
} entry;

typedef struct
{
  long nr;
  entry *entry;
} database;

typedef struct
{
  entry entry;
  int *exists;
} conditions;

/*************************************************************************
|  FILE : HM_symbol.c
**************************************************************************/

void display_HM_symbol(const char *qname,
                       int zname1,
                       int zname2,
                       MP_INT *aff_name);

/*************************************************************************
|  FILE : Q_catalog.c
**************************************************************************/
extern void (*(display_element [NR_OF_ELEMENTS_IN_EACH_ENTRY])) (entry *data);
extern void (*(load_element [NR_OF_ELEMENTS_IN_EACH_ENTRY])) (const char *string, entry *data);
extern void (*(delete_element [NR_OF_ELEMENTS_IN_EACH_ENTRY])) (entry *data);
extern int (*(compare_element [NR_OF_ELEMENTS_IN_EACH_ENTRY])) (entry *data1, entry *data2
);
extern const char *name_element [NR_OF_ELEMENTS_IN_EACH_ENTRY];

void apply_cond_to_display_list (conditions *cond,
                                 database *database,
                                 int display_list[],
                                 int new_condition);

void unapply_cond_to_display_list (database *database,
                                   int display_list[],
                                   int unset_condition);

void display_data_list (database *datas,
                        int display_list[]);

database *load_database (const char *filename,
                         int degree);

void free_database (database *datas);

/*************************************************************************
|  FILE : aff_class_inf.c
**************************************************************************/

matrix_TYP *aff_class_inf(bravais_TYP *R,
                          bravais_TYP *DATAZ,
                          matrix_TYP *PRES,
                          MP_INT *aff_name,
                          bravais_TYP **RC);

void extend(matrix_TYP *T);

bravais_TYP *space_group_from_matrix(bravais_TYP *G,
                                     matrix_TYP *x,
                                     matrix_TYP *cocycle,
                                     matrix_TYP *D);


/*************************************************************************
|  FILE: compute_q_matrix.c
**************************************************************************/

matrix_TYP *compute_q_matrix (bravais_TYP *G);

/*************************************************************************
|  FILE: point_group.c
**************************************************************************/

bravais_TYP *point_group(bravais_TYP *R,
                         int opt);

/*************************************************************************
|  FILE: q_class_inf.c
**************************************************************************/

matrix_TYP *q_class_inf (bravais_TYP *G,
                         database *database,
                         char *qclass_name,
                         char *symbol,
                         bravais_TYP **OUT,
                         matrix_TYP **PRES,
                         int transformation);

/*************************************************************************
|  FILE: z_class_inf.c
**************************************************************************/

matrix_TYP *z_class_inf(bravais_TYP *G,
                        bravais_TYP *DATABASEGROUP,
                        bravais_TYP **RES,
                        int *name);
			
/*************************************************************************
|  FILE: reverse_name_fct.c
**************************************************************************/

bravais_TYP *get_qclass_by_name(const char *name,
                                matrix_TYP **PRES,
                                int dim);

bravais_TYP *get_zclass_by_name(bravais_TYP *G,
                                int *first,
                                int *second,
                                int ignore);

bravais_TYP *split_extension(bravais_TYP *G);

bravais_TYP *get_affine_class_by_name(bravais_TYP *G,
                                     matrix_TYP *PRES,
                                     MP_INT *aff_name,
                                     int check);

bravais_TYP *reverse_name(const char *qname,
                          int zname[2],
			  MP_INT aff_name,
			  int i,
			  boolean iflag,
			  char **affstring);

#endif
