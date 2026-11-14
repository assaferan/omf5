omf5 - Orthogonal Modular Forms of Rank 5
=========================================

A library for computing Hecke matrices and eigensystems for positive definite rational quinary quadratic forms.

The code in this repository is based on the work of many, among them:

1. Gonzalo Tornaria - His fast implementation for discriminant 61 is the basis for this code.
2. Jeffery P. Hein - his repository of ternary_birch contains the case of rank 3 and we have widely used ideas from there.
3. Matthew Greenberg and John Voight - their magma code for p-neighbors is the algorithm on which this code is based. (https://github.com/assaferan/ModFrmAlg)

We also use the CARAT package for working with lattices

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Data](#data)
- [Contributing](#contribution)

## Requirements

- [autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html) (autoconf, automake and libtool)
- a C compiler, such as [gcc](https://gcc.gnu.org/) or [clang](https://clang.llvm.org/)
- [GMP](https://gmplib.org/)
- [flint](https://github.com/flintlib/flint), version 3.0 or later
  (flint 2.9 also works, together with [antic](https://github.com/flintlib/antic);
  since flint 3.0 the number field code of antic is part of flint itself.
  Versions of flint older than 2.9 are not supported)
- [CARAT](https://github.com/lbfm-rwth/carat)

On x86-64 the code is compiled with `-msse3 -mavx`, so it needs a processor with AVX.

The continuous integration builds and tests the code on:

- Ubuntu 24.04 with flint 3.0.1, compiled with both gcc and clang
- Debian 12 with flint 2.9.0 and antic 0.2.5, compiled with gcc
- macOS with flint from Homebrew, compiled with Apple clang

It has also been tested on macOS 26 on Apple M4 Pro, with flint 3.6.0 and Apple clang 21.

On macOS, the requirements other than CARAT can be installed with [Homebrew](https://brew.sh/):

    brew update && brew install autoconf automake libtool gmp flint

On Ubuntu 24.04 or later, with:

    sudo apt-get install autoconf automake libtool make gcc libgmp-dev libflint-dev

On Debian 12, whose flint is 2.9, add `libantic-dev`. Older releases package versions of flint older than 2.9, so flint has to be built from source there.

CARAT is built from source, for example from the folder containing it:

    cd carat && ./autogen.sh && ./configure && make && make install && cp libfunctions.a /usr/local/lib && cd ..

CARAT does not install libfunctions.a, hence the copy above. Alternatively, skip the copy and pass the CARAT source directory to omf5 as CARAT_PATH (see below). omf5 only needs libfunctions.a, which `make libfunctions.a` builds on its own.

## Installation

To build as a C library on a standard Linux/Mac OS system with autotools:

    ./autogen.sh
    ./configure --prefix=DIR
    make
    make install

The locations of the dependencies default to /usr/local, and can be set with the variables FLINT_PATH, GMP_PATH, CARAT_PATH (and ANTIC_PATH for flint 2.9). For example, on macOS with Homebrew on Apple silicon:

    FLINT_PATH=/opt/homebrew GMP_PATH=/opt/homebrew CARAT_PATH=/path/to/carat ./configure

Or all at once:

    ./autogen.sh && ./configure && make && make install

To check the build, run

    make check

This runs the standard tests (`src/omf5 -tests`) and checks the examples below (`scripts/check_examples.sh`).

## Usage

The executable is src/omf5.
It runs with command-line arguments as follows

src/omf5 [-tests] [-quad=Q] [-format=f] [-prec=L] [-hecke] [-row] [-p=p] [-genus=g] [-isom] [-disc=d] [-cond=c] [-nonlifts] [-idxs=form_idxs]

where the  arguments are:

The genus can be specified in one of two ways. Either via Q and f - 

[Q] is the quinary quadratic form (lattice) given as 15 comma-separated integers in a format specified by f,

[f] is either 'GG' or 'A', the former for the format in Rama-Tornaria webpage, the latter for the Magma format,

in which case, the genus will be computed using p-neighbors, or via g and d - 

[g] is the name of a file containing the list of genera,

[d] is the discriminant of the lattice, so that g[d] is the relevant genus,

[p] is a prime at which to compute the Hecke matrix/eigenvalue, 

[L] is the precision up to which to compute the hecke matrices/eigenvalues (a_p for p <= L and a_{p^2} for p^2 <= L),

[c] is the conductor of the spinor norm character. If not specified, the program will compute all of them. At the moment, only relevant for computing a column of the Hecke matrix.

[form_idxs] is a comma-separated list of indices of forms in the space for which to compute the hecke eigenvalues. If not specified, computes for all of them.

If either L or p is not supplied, only decomposes the space, and finds Hecke eigenvectors.

If the flag -hecke is supplied, computes a column of the Hecke matrix, otherwise computes the Hecke eigenvalues of forms that are non-lifts. If p is not supplied, computes the Hecke matrix of 
the first prime not dividing the discriminant.

If the flag -row is supplied in addition to -hecke, computes a single row.

If the flag -isom is supplied in addition to -genus=g, read from the genus file, in addition to the list of genera, a list of fixed isometries (over QQ) between the different lattices in the genus.

The genus file is a list of genera, and the genus read is the one at index d, counting from 0. A genus is a list of lattices, and may be an empty list if d is not 0. Each lattice is a list of 15 integers in the format 'GG', as in this file with a genus of two lattices at index 1

[[],[[1,0,-1,-1,1,1,-1,1,0,2,1,-1,2,-2,2],[1,-1,0,1,1,1,0,-1,0,1,1,0,2,1,3]]]

With -isom, each lattice is instead a list of two 5x5 matrices, given by their rows - the Gram matrix of the lattice, and an isometry from the first lattice of the genus to it, whose entries are rational numbers such as -2/3.

The lattices have to be a full set of representatives for the classes of a single genus. The program exits with an error if one of them is not positive definite and even, if they do not all have the same discriminant, if two of them are isometric, if their mass is not the mass of the genus, or if one of them is not in the genus of the first lattice. With -isom the last of these is not checked, and each isometry is checked to be an isometry from the first lattice, with a denominator which is prime to the discriminant.

If the flag -tests is supplied, additionally runs standard tests.

If the flag -nonlifts is supplied, computes eigenvalues only for the forms which are not lifts.

Each form is printed with its automorphic type, `aut_rep_type`, which is one of:

- `G` - the form is not a lift,
- `P` - a Saito-Kurokawa lift,
- `Y` - a Yoshida lift,
- `F` - the Eisenstein series, which is the constant function,
- `O` - an old form, for which no eigenvector was found.

A form is taken to be a lift if its L-polynomial factors, in the same way, at each of the primes p not dividing the discriminant with p^2 <= L. If there are fewer than two such primes, or L is not supplied, the first two primes not dividing the discriminant are used. So the types are more reliable when L is larger. For example, for discriminant 249 there is a form which is not a lift, whose L-polynomial factors at 2 and at 5, but not at 7. It is reported as a Yoshida lift, and so is left out when -nonlifts is supplied, unless L is at least 49.

If form_idxs is supplied the lifts are not detected. Instead, the forms at these indices are given the type `G`, and the eigenvalues are computed only for them, whether or not -nonlifts is supplied. The indices are counted over all the conductors, in the order in which the forms are printed, starting from 0. They can be given in any order. An index at which there is no form is ignored, with a warning.
    
The progress messages, the eigenvectors and the traces of the Hecke eigenvalues are printed to stderr. The Hecke eigenvalues themselves are printed to stdout, as a Python dictionary mapping each conductor to the list of its forms.

Example runs:

1. Decomposition:

```
> src/omf5 -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG > /dev/null
computing genus took 0.026693 sec
recomputing genus took 0.033674 sec
The possible conductors are:
1 61
The corresponding dimensions are: 8 0
constructing eigenvector for poly
x+7
constructing eigenvector for poly
x-15
constructing eigenvector for poly
x^6-29*x^5+322*x^4-1714*x^3+4471*x^2-5205*x+2026
computing eigenvectors took 0.006443 sec
For conductor 1, eigenvectors are:
0 -6 6 4 12 0 -12 -3 over Number field with defining polynomial x+7
1 1 1 1 1 1 1 1 over Number field with defining polynomial x-15
4*a^5 - 88*a^4 + 700*a^3 - 2448*a^2 + 3648*a - 1880 3*a^4 - 50*a^3 + 272*a^2 - 558*a + 397 4*a^4 - 73*a^3 + 453*a^2 - 1091*a + 867 -2*a^4 + 32*a^3 - 150*a^2 + 172*a + 12 14*a^3 - 194*a^2 + 786*a - 862 -a^5 + 20*a^4 - 147*a^3 + 475*a^2 - 620*a + 241 4*a^3 - 40*a^2 + 96*a - 92 -a^4 + 19*a^3 - 123*a^2 + 329*a - 288 over Number field with defining polynomial x^6-29*x^5+322*x^4-1714*x^3+4471*x^2-5205*x+2026
For conductor 61, eigenvectors are:
Found nonlift at total_idx = 0, c = 0, form_idx = 0
Found lift of type F at total_idx = 1, c = 0, form_idx = 1
Found lift of type P at total_idx = 2, c = 0, form_idx = 2
For conductor 1:
traces of hecke eigenvalues of T_{p^1} are:

computing eigenvalues for k = 1, took 0.000003 sec
traces of hecke eigenvalues of T_{p^2} are:

computing eigenvalues for k = 2, took 0.000002 sec
traces of hecke eigenvalues of T_{p^1} are:

computing eigenvalues for k = 1, took 0.000003 sec
traces of hecke eigenvalues of T_{p^2} are:

computing eigenvalues for k = 2, took 0.000002 sec
traces of hecke eigenvalues of T_{p^1} are:

computing eigenvalues for k = 1, took 0.000002 sec
traces of hecke eigenvalues of T_{p^2} are:

computing eigenvalues for k = 2, took 0.000003 sec
For conductor 61:
```

Its output on stdout is

```
{1 : [{ 'aut_rep_type' : 'G', 'field_poly' : '2  7 1', 'lambda_p' : [], 'lambda_p_square' : [], }, { 'aut_rep_type' : 'F', 'field_poly' : '2  -15 1', 'lambda_p' : [], 'lambda_p_square' : [], }, { 'aut_rep_type' : 'P', 'field_poly' : '7  2026 -5205 4471 -1714 322 -29 1', 'lambda_p' : [], 'lambda_p_square' : [], }, ], 61 : [], }
```

so that the first form is not a lift, the second is the Eisenstein series, and the third is a Galois orbit of 6 Saito-Kurokawa lifts.

2. Hecke eigenvalues, for primes up to 100. The output to stderr starts as in the previous example, and continues with:

```
> src/omf5 -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG -prec=100 > /dev/null
...
For conductor 61, eigenvectors are:
Found nonlift at total_idx = 0, c = 0, form_idx = 0
Found lift of type F at total_idx = 1, c = 0, form_idx = 1
Found lift of type P at total_idx = 2, c = 0, form_idx = 2
For conductor 1:
traces of hecke eigenvalues of T_{p^1} are:
-7 -3 3 -9 -4 -3 37 -75 10 212 -6 -88 -3 547 -147 -108 -45 -632 -650 859 -978 931 -571 453
computing eigenvalues for k = 1, took 2.498523 sec
traces of hecke eigenvalues of T_{p^2} are:
7 -9 -9 -42
computing eigenvalues for k = 2, took 0.045981 sec
traces of hecke eigenvalues of T_{p^1} are:
15 40 156 400 1464 2380 5220 7240 12720 25260 30784 52060 70644 81400 106080 151740 208920 305320 363024 394420 499360 578760 712980 922180
computing eigenvalues for k = 1, took 2.310008 sec
traces of hecke eigenvalues of T_{p^2} are:
30 120 780 2800
computing eigenvalues for k = 2, took 0.046518 sec
traces of hecke eigenvalues of T_{p^1} are:
29 59 164 309 612 1038 1786 2287 3235 4757 5717 8337 9491 11542 13464 16843 19388 27780 29755 34353 38362 40541 47818 61861
computing eigenvalues for k = 1, took 2.335693 sec
traces of hecke eigenvalues of T_{p^2} are:
-3 -4 48 72
computing eigenvalues for k = 2, took 0.047114 sec
For conductor 61:
```

The primes dividing the discriminant are skipped, so that above p runs over the primes up to 100 other than 61, and p^2 over the squares of 2, 3, 5 and 7.

With -nonlifts added, only the first of the three forms is computed.

3. Hecke eigenvalues, for discriminant 167, where there are forms of conductor 167 as well:

```
> src/omf5 -quad=1,1,0,1,0,1,0,0,0,1,1,0,1,1,34 -format=GG -prec=31 > /dev/null
computing genus took 0.158128 sec
recomputing genus took 0.200698 sec
The possible conductors are:
1 167
The corresponding dimensions are: 19 1
constructing eigenvector for poly
x+2
constructing eigenvector for poly
x-15
constructing eigenvector for poly
x^2+3*x-1
constructing eigenvector for poly
x^15-83*x^14+3136*x^13-71421*x^12+1094296*x^11-11921809*x^10+95169836*x^9-565289347*x^8+2511011707*x^7-8311069876*x^6+20243201600*x^5-35448149896*x^4+42917301824*x^3-33621134688*x^2+15050560000*x-2845000192
constructing eigenvector for poly
x+8
computing eigenvectors took 0.028038 sec
For conductor 1, eigenvectors are:
0 0 -2 6 -4 1 0 4 0 6 4 -4 -24 0 6 0 -4 0 0 over Number field with defining polynomial x+2
-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 over Number field with defining polynomial x-15
0 0 3*a-6 -9*a+18 -12*a+6 3*a+3 18 -6*a-6 0 -9*a-9 12*a-6 -3*a-3 36*a+36 -54 18*a+18 0 -12*a-12 0 0 over Number field with defining polynomial x^2+3*x-1
a vector of length 19 over a number field of degree 15
For conductor 167, eigenvectors are:
1 over Number field with defining polynomial x+8
Found nonlift at total_idx = 0, c = 0, form_idx = 0
Found lift of type F at total_idx = 1, c = 0, form_idx = 1
Found nonlift at total_idx = 2, c = 0, form_idx = 2
Found lift of type P at total_idx = 3, c = 0, form_idx = 3
Found nonlift at total_idx = 4, c = 1, form_idx = 0
For conductor 1:
traces of hecke eigenvalues of T_{p^1} are:
-2 0 -2 2 -14 -34 -15 16 155 40 -152
computing eigenvalues for k = 1, took 0.135079 sec
traces of hecke eigenvalues of T_{p^2} are:
2 -17 16
computing eigenvalues for k = 2, took 0.015046 sec
traces of hecke eigenvalues of T_{p^1} are:
15 40 156 400 1464 2380 5220 7240 12720 25260 30784
computing eigenvalues for k = 1, took 0.557155 sec
traces of hecke eigenvalues of T_{p^2} are:
30 120 780
computing eigenvalues for k = 2, took 0.076117 sec
traces of hecke eigenvalues of T_{p^1} are:
-3 -9 2 3 92 -41 95 -80 189 -220 -2
computing eigenvalues for k = 1, took 0.136283 sec
traces of hecke eigenvalues of T_{p^2} are:
-3 12 -28
computing eigenvalues for k = 2, took 0.015261 sec
traces of hecke eigenvalues of T_{p^1} are:
83 167 424 755 1941 2588 4436 5485 8130 12589 14431
computing eigenvalues for k = 1, took 0.560784 sec
traces of hecke eigenvalues of T_{p^2} are:
24 68 204
computing eigenvalues for k = 2, took 0.077537 sec
For conductor 167:
traces of hecke eigenvalues of T_{p^1} are:
-8 -10 -4 -14 -22 -4 -47 -12 41 50 -504
computing eigenvalues for k = 1, took 2.972527 sec
traces of hecke eigenvalues of T_{p^2} are:
10 11 -44
computing eigenvalues for k = 2, took 0.035555 sec
```

## Data

After building the executable, it can be used to generate data. The scripts folder contains several instances of data that can be produced.

Most of them require updating the "export LD_LIBRARY_PATH" statement and the existence of folders named **data** and **logs** for keeping the data and log files produced.

The shell script **ev_run.sh** can be used to compute all hecke eigenvalues for discriminant up to 1000 and primes up to 100.

The python script **ev_run_nonlifts.py** can be used to compute all hecke eigenvalues for specific forms inside each space, specified in the file **nonlifts_idx.dat**.
This was used after discarding the oldforms to compute additional eigenvalues only for the newforms which are non lifts. 

The shell script **test_runs.sh** runs sample tests, whose timings were recorded in the paper. 

The shell script **check_examples.sh** checks the Hecke eigenvalues computed for the examples above against their known values. It is run by `make check`, after the standard tests.

The python script **hecke_run.py** computes the full hecke matrices for primes up to 10 and its first row for primes up to 100 for all discriminants up to 1000.

The shell script **parse_omf_runs.sh** parses the output of the runs computing eigenvalues in order to produce data to be transferred to the LMFDB.

## Contribution

If you want to help develop this project, please create your own fork on Github and submit a pull request. I will do my best to integrate any additional useful features as necessary. Alternatively, submit a patch to me via email at assaferan@gmail.com.
