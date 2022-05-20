#ifndef CARAT_TOOLS_H
#define CARAT_TOOLS_H

#include "typedef.h"

/*-------------------------------------------------------------*\
| FILE: carat_exit.c 
\*-------------------------------------------------------------*/
extern void carat_exit(const char *fmt, ...);

/*-------------------------------------------------------------*\
| FILE: chin_remainder.c 
\*-------------------------------------------------------------*/
extern int chin_remainder(int x1, int x2, int p1, int p2);

/*-------------------------------------------------------------*\
| FILE: intpow.c 
\*-------------------------------------------------------------*/
extern int intpow(int a,int b);

/*-------------------------------------------------------------*\
| FILE: itoasc.c 
\*-------------------------------------------------------------*/
extern void itoasc(int n, char s[]);

/*-------------------------------------------------------------*\
| FILE: malloc2dim.c 
\*-------------------------------------------------------------*/
extern char **calloc2dim(int r,int c,int size);
extern char **malloc2dim(int r,int c,int size);
extern void memcpy2dim(char **dest, const char **src, int r,int c,int size);
extern void memset2dim(char **dest, int r, int c, int size, const char *value);
extern void free2dim(char **old, int rows);

/*-------------------------------------------------------------*\
| FILE: mindiv.c 
\*-------------------------------------------------------------*/
extern int min_div( int a, int b) ;

/*-------------------------------------------------------------*\
| FILE: ovfl_mul.c 
\*-------------------------------------------------------------*/
extern int ovfl_mul( int a, int b);

/*-------------------------------------------------------------*\
| FILE: prime_tools.c 
\*-------------------------------------------------------------*/
extern int act_prime;

extern int (*S)(int, int);
extern int (*P)(int, int);
extern void cleanup_prime();
extern void init_prime ( int prime);

/*-------------------------------------------------------------*\
| FILE: ramdom.c
\*-------------------------------------------------------------*/
extern int random_own();

/*-------------------------------------------------------------*\
| FILE: tools.c 
\*-------------------------------------------------------------*/
extern rational Zero;
extern rational One;

extern int GGT (int _a, int _b);
extern int KGV(  int a,  int b );
extern void Normal (rational *a);
extern void Normal2 ( int *_z, int *_n );
extern void rat_add( int *az, int *an, int bz, int bn );
extern int *factorize_new( int zahl, int *erg);
extern int *factorize( int zahl);
extern void gcd_darstell(int a1, int a2, int *v1, int *v2, int *gcd);
extern int p_inv(int a, int p);
extern int signum(int a);

#endif
