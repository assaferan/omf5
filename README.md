omf5 - Orthogonal Modular Forms of Rank 5
=========================================

A library for computing Hecke matrices for positive definite rational quinary quadratic forms.

The code in this repository is based on the work of many, among them:

1. Gonzalo Tornaria - His fast implementation for discriminant 61 is the basis for this code.
2. Jeffery P. Hein - his repository of ternary_birch contains the case of rank 3 and we have widely used ideas from there.
3. Matthew Greenberg and John Voight - their magma code for p-neighbors is the algorithm on which this code is based. (https://github.com/assaferan/ModFrmAlg)

We also use the CARAT package for working with lattices

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contribution)

## Requirements

- [autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html)
- [gcc](https://gcc.gnu.org/)

* This code has only been tested so far on:
- macOS Monterey 12.3.1, compiled with gcc 4.2.1 on Intel i9 8-Core Processor.

The main (and develop) branch require the following libraries:

- [CARAT](https://github.com/lbfm-rwth/carat)
- [flint](https://github.com/wbhart/flint2)
- [antic](https://github.com/wbhart/antic)

Installation of all the requirements on macOS can be done as follows, from a folder containing the above 3 libraries.

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

The executable is bin/omf5.
It runs with command-line arguments as follows

bin/omf5 [-tests] [-quad=Q] [-format=f] [-prec=L] [-form_idx=idx]

where the  arguments are:

    [Q] - the quinary quadratic form (lattice) given as 15 comma-separated integers in a format specified by f, currently one of two.
    [f] - either 'GG' or 'A', the former for the format in Rama-Toranria webpage, the latter for the Magma format, as in the ModFrmAlg package.
    [L] - the preicision up to which to compute the hecke eigenvalues (a_p for p <= L and a_{p^2} for p^2 <= L).
    [idx]- the index of the form in the decomposition to eigenvectors.
    If either L or i is not supplied, only decomposes the space, and finds Hecke eigenvectors.
    If the flag -tests is supplied, additionally runs standard tests.

Example runs:

1. Decomposition:
> src/omf5 -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG
Expected cost of isometries is 72.425000
Expected cost of reduced isometries is 441.550000
Expecting average number of 1.354839 calls to is_isometric.
Expecting average number of 0.129032 calls to is_isometric.
Expecting average number of 0.000000 calls to is_isometric.
Expecting average number of -0.000000 calls to is_isometric.
Recalibrated with theta_prec = 3 and red_on_isom = 0 
computing genus took 0.000000
computing eigenvectors took 0.000000
eigenvectors are:
0 -4 0 3 12 6 -6 -12 over Number field with defining polynomial x+7
1 1 1 1 1 1 1 1 over Number field with defining polynomial x-15
-4*a^5+88*a^4-704*a^3+2548*a^2-4224*a+2584 a^5-24*a^4+212*a^3-842*a^2+1459*a-870 3*a^4-45*a^3+228*a^2-441*a+255 2*a^4-39*a^3+256*a^2-679*a+620 2*a^4-46*a^3+338*a^2-886*a+592 -4*a^4+75*a^3-479*a^2+1205*a-925 -4*a^4+80*a^3-565*a^2+1628*a-1571 -a^5+28*a^4-302*a^3+1522*a^2-3447*a+2728 over Number field with defining polynomial x^6-29*x^5+322*x^4-1714*x^3+4471*x^2-5205*x+2026
traces of hecke eigenvalues are:

traces of hecke eigenvalues T_p^2 are:

computing eigenvalues took 0.000000

2. Hecke eigenvalues:
> src/omf5 -quad=1,0,0,1,1,1,0,1,0,1,0,0,1,0,8 -format=GG -form_idx=0 -prec=100
Expected cost of isometries is 70.025000
Expected cost of reduced isometries is 442.337500
Expecting average number of 1.354839 calls to is_isometric.
Expecting average number of 0.129032 calls to is_isometric.
Expecting average number of 0.000000 calls to is_isometric.
Expecting average number of -0.000000 calls to is_isometric.
Recalibrated with theta_prec = 3 and red_on_isom = 0 
computing genus took 0.000000
computing eigenvectors took 0.000000
traces of hecke eigenvalues are:
-7 -3 3 -9 -4 -3 37 -75 10 212 -6 -88 -3 547 -147 -108 -45 145 -632 -650 859 -978 931 -571 453 
traces of hecke eigenvalues T_p^2 are:
7 -9 -9 -42 
computing eigenvalues took 124.000000

## Contribution

If you want to help develop this project, please create your own fork on Github and submit a pull request. I will do my best to integrate any additional useful features as necessary. Alternatively, submit a patch to me via email at assaferan@gmail.com.
