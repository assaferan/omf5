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

- [autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html)
- [gcc](https://gcc.gnu.org/)

This code has only been tested so far on:

macOS Monterey 12.3.1, compiled with clang 12.0.0 on Intel i9 8-Core Processor.

macOS Catalina 10.15.7, compiled with clang 12.0.0 on Intel Core i5 Quad-Core Processor.

linux Ubuntu 22.04.2, compiled with gcc 11.3.0 on AMD Ryzen ThreadRipper 2970WX 24-Core Processor.

The main (and develop) branch require the following libraries:

- [CARAT](https://github.com/lbfm-rwth/carat)
- [flint](https://github.com/wbhart/flint2)
- [antic](https://github.com/wbhart/antic)

Installation of all the requirements on either linux or macOS can be done as follows, from a folder containing the above 3 libraries.

    brew update && brew install autoconf automake libtool
    cd flint2 && ./configure && make && make check && make install && cd ..
    cd antic && ./configure && make && make install && make check && cd ..
    cd carat && ./autogen.sh && ./configure && make && make install && cp libfunctions.a /usr/local/lib && cd ..

## Installation

To build as a C library on a standard Linux/Mac OS system with autotools:

    ./autogen.sh
    ./configure --prefix=DIR
    make
    make install

Or all at once:

   ./autogen.sh && ./configure && make && make install

## Usage

The executable is src/omf5.
It runs with command-line arguments as follows

src/omf5 [-tests] [-quad=Q] [-format=f] [-prec=L] [-hecke] [-row] [-p=p] [-genus=g] [-isom] [-disc=d] [-cond=c] [-nonlifts] [-idxs=form_idxs]

where the  arguments are:

The genus can be specified in one of two ways. Either via Q and f - 

[Q] is the quinary quadratic form (lattice) given as 15 comma-separated integers in a format specified by f,

[f] is either 'GG' or 'A', the former for the format in Rama-Toranria webpage, the latter for the Magma format,

in which case, the genus will be computed using p-neighbors, or via g and d - 

[g] is the name of a file containing the list of genera,

[d] is the discriminant of the lattice, so that g[d] is the relevant genus,

[p] is a prime at which to compute the Hecke matrix/eigenvalue, 

[L] is the preicision up to which to compute the hecke matrices/eigenvalues (a_p for p <= L and a_{p^2} for p^2 <= L),

[c] is the conductor of the spinor norm character. If not specified, the program will compute all of them. At the moment, only relevant for computing a column of the Hecke matrix.

[form_idxs] is a list of indices of forms in the space for which to compute the hecke eigenvalues. If not specified, computes for all of them.

If either L or p is not supplied, only decomposes the space, and finds Hecke eigenvectors.

If the flag -hecke is supplied, computes a column of the Hecke matrix, otherwise computes the Hecke eigenvalues of forms that are non-lifts. If p is not supplied, computes the Hecke matrix of 
the first prime not dividing the discriminant.

If the flag -row is supplied in addition to -hecke, computes a single row.

If the flag -isom is supplied in additoin to -genus=g, read from the genus file, in addition to the list of genera, a list of fixed isometries (over QQ) between the different lattices in the genus.

If the flag -tests is supplied, additionally runs standard tests.

If the flag -nonlifts is supplied, computes eigenvalues only for the forms which are not lifts.
    
Example runs:

1. Decomposition:
> src/omf5 -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG
computing genus took 0.088564 sec
recomputing genus took 0.120246 sec
The possible conductors are: 
1 61 
The corresponding dimensions are: 8 0 
computing eigenvectors took 0.041764 sec
For conductor 1, eigenvectors are:
0 -6 6 4 12 0 -12 -3 over Number field with defining polynomial x+7
1 1 1 1 1 1 1 1 over Number field with defining polynomial x-15
4*a^5-88*a^4+700*a^3-2448*a^2+3648*a-1880 3*a^4-50*a^3+272*a^2-558*a+397 4*a^4-73*a^3+453*a^2-1091*a+867 -2*a^4+32*a^3-150*a^2+172*a+12 14*a^3-194*a^2+786*a-862 -a^5+20*a^4-147*a^3+475*a^2-620*a+241 4*a^3-40*a^2+96*a-92 -a^4+19*a^3-123*a^2+329*a-288 over Number field with defining polynomial x^6-29*x^5+322*x^4-1714*x^3+4471*x^2-5205*x+2026
For conductor 61, eigenvectors are:
For conductor 1:
traces of hecke eigenvalues of T_{p^1} are:

computing eigenvalues for k = 1, took 0.000013 sec
traces of hecke eigenvalues of T_{p^2} are:

computing eigenvalues for k = 2, took 0.000002 sec
For conductor 61:

2. Hecke eigenvalues:
> src/omf5 -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG -prec=100
computing genus took 0.090954 sec
recomputing genus took 0.120304 sec
The possible conductors are: 
1 61 
The corresponding dimensions are: 8 0 
computing eigenvectors took 0.041667 sec
For conductor 1, eigenvectors are:
0 -6 6 4 12 0 -12 -3 over Number field with defining polynomial x+7
1 1 1 1 1 1 1 1 over Number field with defining polynomial x-15
4*a^5-88*a^4+700*a^3-2448*a^2+3648*a-1880 3*a^4-50*a^3+272*a^2-558*a+397 4*a^4-73*a^3+453*a^2-1091*a+867 -2*a^4+32*a^3-150*a^2+172*a+12 14*a^3-194*a^2+786*a-862 -a^5+20*a^4-147*a^3+475*a^2-620*a+241 4*a^3-40*a^2+96*a-92 -a^4+19*a^3-123*a^2+329*a-288 over Number field with defining polynomial x^6-29*x^5+322*x^4-1714*x^3+4471*x^2-5205*x+2026
For conductor 61, eigenvectors are:
For conductor 1:
traces of hecke eigenvalues of T_{p^1} are:
-7 -3 3 -9 -4 -3 37 -75 10 212 -6 -88 -3 547 -147 -108 -45 145 -632 -650 859 -978 931 -571 453 
computing eigenvalues for k = 1, took 111.827940 sec
traces of hecke eigenvalues of T_{p^2} are:
7 -9 -9 -42 
computing eigenvalues for k = 2, took 0.178717 sec
For conductor 61:

> src/omf5 -quad=1,1,0,1,0,1,0,0,0,1,1,0,1,1,34 -format=GG -prec=31  
computing genus took 0.470110 sec
recomputing genus took 0.640326 sec
The possible conductors are: 
1 167 
The corresponding dimensions are: 19 1 
computing eigenvectors took 0.146885 sec
For conductor 1, eigenvectors are:
0 0 2 -6 -6 0 0 4 -1 -4 -4 4 24 0 -6 0 0 4 0 over Number field with defining polynomial x+2
-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 over Number field with defining polynomial x-15
0 0 -3*a+6 9*a-18 9*a+9 -18 0 12*a-6 -3*a-3 6*a+6 -12*a+6 3*a+3 -36*a-36 54 -18*a-18 0 0 12*a+12 0 over Number field with defining polynomial x^2+3*x-1
a vector of length 19 over a number field of degree 15
For conductor 167, eigenvectors are:
1 over Number field with defining polynomial x+8
For conductor 1:
traces of hecke eigenvalues of T_{p^1} are:
-2 0 -2 2 -14 -34 -15 16 155 40 -152 
computing eigenvalues for k = 1, took 9.674285 sec
traces of hecke eigenvalues of T_{p^2} are:
2 -17 16 
computing eigenvalues for k = 2, took 0.119902 sec
traces of hecke eigenvalues of T_{p^1} are:
-3 -9 2 3 92 -41 95 -80 189 -220 -2 
computing eigenvalues for k = 1, took 9.684862 sec
traces of hecke eigenvalues of T_{p^2} are:
-3 12 -28 
computing eigenvalues for k = 2, took 0.119482 sec
For conductor 167:
traces of hecke eigenvalues of T_{p^1} are:
-8 -10 -4 -14 -22 -4 -47 -12 41 50 -504 
computing eigenvalues for k = 1, took 27.231407 sec
traces of hecke eigenvalues of T_{p^2} are:
10 11 -44 
computing eigenvalues for k = 2, took 0.292487 sec

## Data

After building the executable, it can be used to generate data. The scripts folder contains several instances of data that can be produced.

Most of them require updating the "export LD_LIBRARY_PATH" statement and the exisence of folders named **data** and **logs** for keeping the data and log files produced.

The shell script **ev_run.sh** can be used to compute all hecke eigenvalues for discriminant up to 1000 and primes up to 100.

The python script **ev_run_nonlifts.py** can be used to compute all hecke eigenvalues for specific forms inside each space, specified in the file **nonlifts_idx.dat**.
This was used after discarding the oldforms to compute additional eigenvalues only for the newforms which are non lifts. 

The shell script **test_runs.sh** runs sample tests, whose timings were recorded in the paper. 

The python script **hecke_run.py** computes the full hecke matrices for primes up to 10 and its first row for primes up to 100 for all discriminants up to 1000.

The shell script **parse_omf_runs.sh** parses the output of the runs computing eigenvalues in order to produce data to be transferred to the LMFDB.

## Contribution

If you want to help develop this project, please create your own fork on Github and submit a pull request. I will do my best to integrate any additional useful features as necessary. Alternatively, submit a patch to me via email at assaferan@gmail.com.
