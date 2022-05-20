#ifdef __cplusplus
   extern "C" {
#endif

#ifndef CARAT_GRAPH_H
#define CARAT_GRAPH_H

#include "typedef.h"

#include "zass.h"

#define TWOTO21 2097152

/* global variables */
extern boolean GRAPH_DEBUG;


/* new structures */
typedef struct {
   bravais_TYP *G;                      /* group we started with */
   bravais_TYP **Z;                     /* Z-classes */
   matrix_TYP ***norm_inv;		/* inverses of normalizer for each Z-class */
   int Z_no;                            /* number of Z-classes */
   bravais_TYP ***aff;                  /* affine classes */
   int *aff_no;                         /* number of affine classes in each Z-class */
   int *first_aff;			/* position of the first affine class in each Z-class
   					   in the list (aff[0][0], aff[0][1], ... aff[1][0], ...) */
   matrix_TYP ***X;			/* matrices with informations about affine classes */
   matrix_TYP **X_2_inv;		/* inverse of ...->X[i][2] for each Z-class i */
   matrix_TYP ***coz;                   /* cocycles of the affine classes */
   MP_INT **names;			/* names for affine classes */
   int **names_int;			/* names for affine classes in integers */
   int *coho_size;			/* size of the cohomology group for each Z-class */
   int all;				/* total number of affine classes */
   matrix_TYP *pres;                    /* presentation */
   QtoZ_TYP *INZ;			/* graph for Z-classes */
   matrix_TYP ****stab_coz;             /* stabilizer in N_{Gl_n(Z)}(G) for each cocycle in each Z-class */
   int **stab_gen_no;                   /* number of generators from this stabilizer */
   int ****WORDS;			/* Words for stabilizer of cocycles for each Z-class */
   int **NUMBER_OF_WORDS;		/* number of words */
   matrix_TYP ***N;			/* representation of normalizer on H^1 */
   matrix_TYP ***gen_inv;		/* inverses of generators for each Z-class */
   boolean l_option;			/* TRUE iff is_option(l) */
   int **list_of_names;			/* for each cocycle in each Z-class: number of the
                                           smallest representative in the orbit of this
                                           cocycle */
} Q_data_TYP;


typedef struct {
   matrix_TYP *D;
   matrix_TYP **M;
   matrix_TYP *i;
   int erz_no;
   int flag;
   int D_first;
} H1_mod_ker_TYP;



typedef struct {
   int con_no;			/* number of conjugacy classes */
   bravais_TYP **groups;	/* representatives for each conjugacy class */
   bahn ***strong;		/* strong generating set for each represenative */
   int ***worte;		/* words for generators of the representatives in the
                                   generators of the group in each conjugacy class */
   int *elem_no;		/* number of elements in the conjugacy classes */
   int *orders;			/* orders of the groups in each conjugacy class */
   int total;			/* total number of subgroups */
} t_sub_TYP;				


/* --------------------------------------------------------------------- */
/* some information about a cocycle                                      */
/* --------------------------------------------------------------------- */
typedef struct
{
   matrix_TYP *coz;		/* cocycle in Q^n */
   matrix_TYP *std_coz;		/* standard form for this cocycle */
   matrix_TYP *diff;		/* coz - std_coz */
   MP_INT name;			/* name for this cocycle (number of element with lowest
                                   number in the orbit) */
   MP_INT number;		/* number for this cocycle */
   int **WORDS;			/* Words for stabilizer of cocycle */
   int WORDS_no;		/* number of words for stabilizer of cocycle */
   matrix_TYP **Stab;		/* Stablizer of the cocycle / Pointgroup !!! */
   int Stab_no;			/* number of elements in the stabilizer */
   matrix_TYP *darst;		/* cocycle represented in C_a x C_b x ... */
   matrix_TYP **aff_transl;	/* return value of transl_aff_normal */
   int aff_transl_no;		/* # of elements in transl_aff_normal */
   matrix_TYP **N_darst;	/* representation of G->normal on H^1(G,Q^n/Z^n) */
   matrix_TYP **N_inv;		/* inverses of G->normal */
   int N_darst_no;		/* G->normal_no */
   matrix_TYP **X;		/* return value of cohomology for H^1(G,Q^n/Z^n) */
   matrix_TYP *GLS;		/* inverse of X[2] */
} coz_TYP;




/* functions */
Q_data_TYP *get_Q_data(bravais_TYP *G,
                       matrix_TYP *pres,
                       boolean l_option);

void free_Q_data(Q_data_TYP *data);

void put_Q_data(Q_data_TYP *data,
                char *groupname,
                int printflag);

bravais_TYP *extract_r(bravais_TYP *G,
                       matrix_TYP *X);

matrix_TYP **extensions_o(matrix_TYP *cocycle,
                          matrix_TYP *D,
                          matrix_TYP *R,
                          bravais_TYP *G,
                          int **lengths,
                          MP_INT **names,
                          int *number_of_orbits,
                          int ****WORDS,
                          int **NUMBER_OF_WORDS,
                          matrix_TYP ***N,
                          MP_INT coho_size,
                          int option,
                          int *list_of_names);

matrix_TYP **all_cocycles(matrix_TYP *relator_input,
                          bravais_TYP *G,
                          int *anzahl,
                          matrix_TYP **matinv,
                          matrix_TYP ***X,
                          MP_INT **names,
                          int ****WORDS,
                          int **NUMBER_OF_WORDS,
                          matrix_TYP ***N,
                          int *coho_size,
                          int **list_of_names,
                          boolean l_option);

matrix_TYP *subgroupgraph(Q_data_TYP *data,
                          boolean oflag);

matrix_TYP **stab_coz(int **words,
                      int no_of_words,
                      matrix_TYP **N,
                      matrix_TYP **N_inv,
                      int anz,
                      int flag,
                      int *stab_gen_no);

matrix_TYP *matrix_on_diagonal(matrix_TYP *mat,
                               int anz);

matrix_TYP *standard_rep(matrix_TYP *coz,
                         matrix_TYP *GLS,
                         matrix_TYP *D);

MP_INT cohomology_size(matrix_TYP *D);

int *aff_classes_in_image(MP_INT *names,
                          int aff_no,
                          int coho_size,
                          matrix_TYP **orbit,
                          matrix_TYP **image_gen,
                          int gen_no,
                          matrix_TYP *D,
                          int *aff_in_image);

int equal_zero(matrix_TYP *m);

int yet_there(matrix_TYP *m,
              matrix_TYP **list,
              int no);

int **stab_lattice(matrix_TYP *lattice,
                   matrix_TYP **N,
                   int N_gen_no,
                   int *no,
                   matrix_TYP **LIST,
                   int anz,
                   boolean *lattice_orbit);

void kernel_and_image(matrix_TYP *phi,
                      matrix_TYP *A,
                      int A_first,
                      matrix_TYP *B,
                      matrix_TYP **kernel,
                      matrix_TYP **image);

void calculate_phi(matrix_TYP *diag,
                   matrix_TYP *coz,
                   matrix_TYP **Xi,
                   matrix_TYP **Xj,
                   matrix_TYP *GLS,
                   matrix_TYP **phi,
                   matrix_TYP **kernel,
                   matrix_TYP ***image,
                   int *image_gen_no,
                   H1_mod_ker_TYP *H1_mod_ker);

matrix_TYP **calculate_S1(matrix_TYP *lattice,
                          matrix_TYP **N,
                          int anz,
                          int *no,
                          matrix_TYP **dataN,
                          matrix_TYP **dataNinv,
                          matrix_TYP *dataX);

int *aufspannen(int coho_size,
                matrix_TYP **elements,
                matrix_TYP **M,
                int gen_no,
                matrix_TYP *D,
                int *anz);

matrix_TYP **col_to_list(matrix_TYP *M);

matrix_TYP **orbit_ker(matrix_TYP **ker_elements,
                       int ker_order,
                       matrix_TYP *D,
                       matrix_TYP **N,
                       int no,
                       int coho_size,
                       int *ker_list,
                       int *anz,
                       int **length);

matrix_TYP **orbit_ksi_plus_ker(matrix_TYP *ksi,
                                matrix_TYP **ker_elements,
                                int ker_order,
                                matrix_TYP *D,
                                matrix_TYP **N,
                                int no,
                                int coho_size,
                                int *ker_list,
                                int *anz,
                                int **length);

void free_H1_mod_ker_TYP(H1_mod_ker_TYP H1_mod_ker);

matrix_TYP **new_representation(matrix_TYP **S,
                                int S_no,
                                H1_mod_ker_TYP H1_mod_ker,
                                matrix_TYP *A);

matrix_TYP *graph_mapped_word(int *w,
                              matrix_TYP **A,
                              matrix_TYP **AINV,
                              matrix_TYP *D);

matrix_TYP ***H1_mod_ker_orbit_alg(H1_mod_ker_TYP H1_mod_ker,
                                   matrix_TYP **S,
                                   int S_no,
                                   int *anz,
                                   int **length,
                                   int ****WORDS,
                                   int **WORDS_no,
				   matrix_TYP *start);

void kernel_elements_2_affine(matrix_TYP **elem,
                              int anz);

matrix_TYP *graph_mat_inv(matrix_TYP *A,
                          matrix_TYP *D,
                          int first);

int word_already_there(int **WORDS,
                       int n);

matrix_TYP *H1_of_standard_to_GL(bravais_TYP *GL,
                                 bravais_TYP *standard,
                                 matrix_TYP **X);

int orbit_on_lattices(matrix_TYP **gitter,
                      int gitter_no,
                      matrix_TYP **N,
                      int N_no,
                      int *list,
                      int *length,
                      int *smallest,
                      matrix_TYP **conj);

int number_of_affine_class(Q_data_TYP *data,
                           matrix_TYP *coz,
                           int i,
                           int flag,
			   boolean wortflag,
			   int **wort);

matrix_TYP **transl_aff_normal(matrix_TYP **erzeuger,
                               int erzanz,
                               int *anzahl);

matrix_TYP **kernel_factor_fct(matrix_TYP **translationen,
                               int translanz,
                               int erz_no,
                               matrix_TYP *lattice,
                               int *no);

bravais_TYP *p_group(bravais_TYP *G);

bravais_TYP **min_k_super(bravais_TYP *P,
                          bravais_TYP *R,
			  matrix_TYP *pres,
                          int *anz,
			  boolean bahnflag,
                          boolean debugflag,
			  int **orbitlength);

bravais_TYP **max_k_sub(bravais_TYP *P,
                        bravais_TYP *R,
                        matrix_TYP *pres,
                        int *anz,
			boolean aflag,
                        boolean debugflag,
			int **orbitlength);

void plus_translationen(bravais_TYP *G,
                        matrix_TYP *mat);

matrix_TYP *extract_c(bravais_TYP *R);

boolean is_k_subgroup(bravais_TYP *S,
                      bravais_TYP *R,
                      bravais_TYP *P,
                      int n,
                      matrix_TYP *pres);

matrix_TYP *sg(bravais_TYP *R,
               bravais_TYP *P);

matrix_TYP *add_mod_D(matrix_TYP *A,
                      matrix_TYP *B,
                      matrix_TYP *D,
                      int first,
                      int diff);

matrix_TYP *to_aff_normal_element(matrix_TYP *lin,
                                  matrix_TYP *coz,
                                  int flag,
                                  bravais_TYP *P,
                                  bravais_TYP *R);

void my_translation(matrix_TYP *TR,
                    matrix_TYP *preimage,
                    matrix_TYP *coz,
                    bravais_TYP *G);

bravais_TYP ****t_subgroups(bravais_TYP *G,
                            matrix_TYP **mats,
                            int no,
                            matrix_TYP *pres,
                            int *aff_no,
                            int **aff_cons,
                            int ***aff_cons_no,
                            bravais_TYP ***R,
                            int flag);
			    
coz_TYP identify_coz(bravais_TYP *G,
                     matrix_TYP *coz,
                     matrix_TYP **X);

void free_coz_TYP(coz_TYP coz_info);

matrix_TYP **calculate_H1(bravais_TYP *G,
                          word *relator,
                          int relatoranz);

void cen_to_norm(bravais_TYP *G);

matrix_TYP **max_sublattices(bravais_TYP *P,
                             int *anz,
                             int **trivialflag,
                             boolean debugflag);

matrix_TYP **calculate_representatives(bravais_TYP *G,
                                       bravais_TYP *GL,
				       matrix_TYP **H_G_Z,
				       matrix_TYP **H_GL_Z,
				       matrix_TYP *lattice,
                                       matrix_TYP *phi,
                                       H1_mod_ker_TYP H1_mod_ker,
				       matrix_TYP *preimage,
				       matrix_TYP *kernel_mat,
				       int *orbit_no,
				       int **orbit_length);

int obergruppenzahl(matrix_TYP *L,
                    matrix_TYP **Norm,
                    matrix_TYP **StabStdCoz,
		    int Stab_anz,
		    int *wort);

#endif /* CARAT_GRAPH_H */

#ifdef __cplusplus
   }
#endif
