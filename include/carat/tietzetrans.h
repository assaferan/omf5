#ifndef CARAT_TIETZETRANS_H
#define CARAT_TIETZETRANS_H

typedef struct {
        matrix_TYP* element;
        int* product;
        int nproduct;
        int left, right;
        } derived_TYP;
typedef struct{
        derived_TYP** list;
        int sizemult;
        int firstfree;} derivedsg_TYP;
#define NIL (-1)

typedef struct
        {
        int* lhsproduct;
        int lhsnproduct;
        int* rhsproduct;
        int rhsnproduct;
        } relator_TYP;

typedef struct 
        {
        derivedsg_TYP* generators;
        relator_TYP* relators;
        int norelators;
        int ext_factor;
        } presentation_TYP;


typedef struct
        {
        int* lhsproduct;
        int lhs_no;
        int* rhsproduct;
        int rhs_no;
        } anne_relator_TYP;

typedef struct 
        {
        matrix_TYP** generators;
        anne_relator_TYP* relators;
    	int gen_no; 
        int rel_no;
        } anne_presentation_TYP;


typedef struct {
int dim;
int P_no;
bravais_TYP *R;
matrix_TYP *T;
matrix_TYP **PRV;
} RG_TYP;

typedef struct {
int dim;
polyeder_TYP *P;
matrix_TYP **BK;
} Fube_TYP;

typedef struct{
matrix_TYP *elm;
int *schreier_vec;
matrix_TYP *trans;
} ele_TYP;

typedef struct {
int *sword;
matrix_TYP *trans;
}Stab_word_TYP;

#endif
