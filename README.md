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

The build defines `NDEBUG`, which removes the assertions of the code. To build with the assertions:

    ./configure CFLAGS=-UNDEBUG

The assertions inside `#ifdef DEBUG` and `#ifdef DEBUG_LEVEL_FULL` stay off. To build with the code and the assertions inside `#ifdef DEBUG` as well:

    ./configure CFLAGS="-UNDEBUG -DDEBUG"

This build prints debugging output as well, and on the standard output some of it is in the middle of the eigenvalues, so its output is not meant to be parsed. To check it, run

    OMF5_DEBUG_BUILD=yes make check

which runs the same examples and checks their traces, types, errors and Hecke matrices. It does not check the eigenvalues printed to the standard output, or that two bases of a lattice give the same output.

## Usage

The executable is src/omf5.
It runs with command-line arguments as follows

src/omf5 [-tests] [-quad=Q] [-format=f] [-prec=L] [-hecke] [-row] [-p=p] [-square] [-genus=g] [-isom] [-disc=d] [-cond=c] [-nonlifts] [-idxs=form_idxs] [-noreduce] [-transport[=M]] [-weight=a[,b]]

where the  arguments are:

The genus can be specified in one of two ways. Either via Q and f - 

[Q] is the quinary quadratic form (lattice) given as 15 comma-separated integers in a format specified by f,

[f] is either 'GG' or 'A', the former for the format in Rama-Tornaria webpage, the latter for the Magma format,

in which case, the genus will be computed using p-neighbors, or via g and d - 

[g] is the name of a file containing the list of genera,

[d] is the discriminant of the lattice, so that g[d] is the relevant genus,

The lattice Q may be given in any basis. When a reduction of the basis lowers the largest diagonal entry of the Gram matrix, the genus is computed from the reduced basis, and a line on the standard error says so: the time grows quickly with this entry. The eigenvalues do not depend on the basis. The order of the classes of the genus does, and with it the Hecke matrices which -hecke prints and the eigenvectors. The lattices of a genus file are reduced in the same way, and the order of the classes is that of the file. With -isom the isometries which are given are taken to the reduced bases, and the Hecke matrices and the eigenvectors do not depend on the bases. Without -isom the program finds the isometries between the lattices. These are isometries over the rationals, which are determined only up to the isometries of the rational quadratic space, and the ones the program finds depend on the bases: on a space of nontrivial conductor this can change the signs of entries of the Hecke matrices and of the eigenvectors. If the flag -noreduce is supplied, the lattice Q, or the lattices of the genus file, are used in the bases in which they are given. In a basis far from reduced, the searches for isometries and automorphisms are still made in a reduced basis, and what they find is taken back to the basis given, since in such a basis they need much time and memory.

[p] is a prime at which to compute the Hecke matrix/eigenvalue, 

[L] is the precision up to which to compute the hecke matrices/eigenvalues (a_p for p <= L and a_{p^2} for p^2 <= L),

[c] is the conductor of the spinor norm character. If not specified, the program will compute all of them. At the moment, only relevant for computing a column of the Hecke matrix.

[form_idxs] is a comma-separated list of indices of forms in the space for which to compute the hecke eigenvalues. If not specified, computes for all of them.

If either L or p is not supplied, only decomposes the space, and finds Hecke eigenvectors.

If the flag -hecke is supplied, computes a column of the Hecke matrix, otherwise computes the Hecke eigenvalues of forms that are non-lifts. If p is not supplied, computes the Hecke matrix of 
the first prime not dividing the discriminant.

If the flag -row is supplied in addition to -hecke, computes a single row.

If the flag -transport is supplied, the Hecke matrices of T_p with -hecke, and the eigenvalues of T_p with -p=p and -prec=L, are computed by carrying the classes of the p-neighbors along isometries with entries in Z[1/q] (issue #63), for the smallest prime q which does not divide the discriminant and is not p. The program first finds the q-neighbors of every class, with an isometry onto their class (one isometry test for each class and isotropic line modulo q). The q-neighbors in the class itself give isometries g of the lattice with entries in Z[1/q], and the class and an isometry of the p-neighbor at g v are found from those of the p-neighbor at v and from the q-neighbors. Usually only one p-neighbor in each orbit of the group generated by the automorphisms and these isometries is tested for isometry (more when the stack of the lines to visit, of at most 4096, is full), and the spinor characters are carried along with them. The output is the same as without the flag: this was checked for the Hecke matrices of every conductor on the 85 lattices of discriminant at most 100, at the primes up to 13, including those which divide the discriminant, and for -prec=50 on the same lattices.

It trades memory for time. The lines modulo p which have been reached are kept in a bitmap, of 2 bits for each of the about p^3 isotropic lines, so the memory grows as p^3 / 4 bytes (4 MB at p = 251, where the run took 14 MB), while without the flag, when the genus has only the trivial conductor, nothing grows with p. The neighbors are counted by class and spinor character as they are found, without an entry for each. The bitmap solves Q(v) = 0 for a coordinate k of v with Q(e_k) prime to p; for a class whose Gram matrix has no such diagonal entry, the lines are kept in a hash set instead, of up to 32 bytes a line, and its column is computed as without the flag if that does not fit in the bound. -transport=M bounds this memory by M MB (256 by default, which allows p up to 1021; -transport=64 allows p up to 643): T_p at a larger prime, and at q, is computed as without the flag, and a line on the standard error says so. The first line on the standard error gives q and the largest prime allowed. T_{p^2} and the weights are always computed as without the flag. Some measurements, with the same output:

| run | without -transport | with -transport |
|---|---|---|
| -hecke -p=83, D = 83 (8 classes) | 13.6 s, 5.5 MB | 6.4 s, 11 MB |
| -hecke -p=23, D = 389 (80 classes) | 20.3 s, 8.6 MB | 4.0 s, 13 MB |
| -prec=100, D = 389 | 272 s, 28 MB | 29 s, 17 MB |
| -hecke -p=97, D = 61 (8 classes) | 10.7 s, 5.0 MB | 9.7 s, 10 MB |
| -hecke -p=251, D = 61 | 182 s, 5.0 MB | 173 s, 14 MB |
| -prec=50, D = 997 (389 classes) | 388 s | 320 s |

The gain is largest when the isometry tests are costly, as in genera with many classes. With -prec=L for a large L, most of the time can go to T_{p^2}, which the flag does not change.

If the flag -square is supplied in addition to -p=p, computes also the eigenvalues of T_{p^2}, the Hecke operator of the p^2-neighbors, and with -hecke its Hecke matrix instead of that of T_p. It is not used with -prec or with -row.

If the flag -weight=a,b is supplied with a >= b >= 0 and a > 0, the forms have the weight (a, b), the highest weight of an irreducible representation of SO(5); -weight=a is the weight (a, 0). Their values are the polynomials f(x, y) in two vectors of 5 variables, of degree a in x and b in y, which are harmonic for the inverse of the Gram matrix of the first lattice in x, in y and in both, and are killed by the sum of the x_i d/dy_i; an isometry g acts by f -> det(g)^(a+b) f(g^(-1) x, g^(-1) y). These form the representation of highest weight (a, b), of dimension (2a+3)(2b+1)(a-b+1)(a+b+2)/6, and for b = 0 they are the harmonic polynomials of degree a. The factor det(g)^(a+b) makes -1 act trivially, so these are the forms of the special orthogonal group, as in Magma's `OrthogonalModularForms(L : weight := [a,b], Special := true)`. They correspond to paramodular forms of weight (k, j) = (a - b + 3, 2b), det^k Sym^j: the weight (a, 0) gives the scalar weight a + 3. (In the coordinates of the symplectic group, as in Assaf-Ladd-Rama-Tornaria-Voight, the weight is (a + b, a - b).) The representation is computed by linear algebra on the polynomials of bidegree (a, b), in coordinates in which the Gram matrix is diagonal, where the conditions split into 16 blocks; the number of these polynomials grows quickly with a and b. A weight (a, b) with b > 0 for which the largest matrices would have more than 2.5 10^7 entries (about 5 GB), such as (7, 7), is an error (the weights (a, 0) are not limited). With -hecke -p=2 on the lattice of discriminant 83, the weight (4, 4) takes about 15 seconds, (5, 5) about a minute and 1 GB, and (6, 6) about 6 minutes and 4 GB. It is supported with -hecke -p=p, which prints the Hecke matrices for each conductor as without the weight, with rational entries; with -p=p, which prints the eigenforms and their eigenvalues at p; and with -prec=L, which prints them at the primes up to L, in both cases as -prec does without the weight. With -square, -hecke -p=p prints the matrices of T_{p^2} instead of T_p, and -p=p prints the eigenvalues of T_{p^2} as well; -prec=L prints them at the primes p with p^2 <= L, as without the weight. The prime of -p=p may divide the discriminant, while -prec=L, as without the weight, gives only the primes which do not (see below for the primes dividing the discriminant). With -prec=L the types of the forms are determined as without the weight, from the factors at p of their spin L-functions (below) at the primes p with p^2 <= L which do not divide the discriminant, when there are three of them at least: a form is `P` (Saito-Kurokawa) or `Y` (Yoshida) when the factor is the product of two factors 1 - alpha x + p^w x^2 at each of these primes, of the same kind at all of them, and `G` otherwise. A Saito-Kurokawa lift, with alpha = p^(a+1) (p + 1), comes from a newform f of weight 2a + 4, with T_p = a_p(f) + p^(a+1) + p^(a+2), and exists only for b = 0. A Yoshida lift has the factors of a newform f of weight 2a + 4 and of a newform g of weight 2b + 2, with T_p = a_p(f) + p^(a-b+1) a_p(g). The type `F` is that of the constant functions, so it is only of the weight 0. For example, on the lattice of discriminant 19 with the weight (1, 0), the form of conductor 19 with T_2 = 10, T_3 = 35 and T_5 = 126 is the lift of the newform of weight 6 and level 19, and on the lattice of discriminant 39 with the weight (1, 1), the form of conductor 3 with T_2 = -16, T_5 = -29 and T_7 = -131 is the Yoshida lift of the newforms of weight 6 and level 3 and of weight 4 and level 13. These were checked for the 38 forms of the weights (1, 0), (2, 0), (3, 0), (1, 1) and (2, 1) on the lattices of discriminant 19, 39, 61 and 83 with -prec=50. Without -prec, and with -p=p, the types are not determined: they are `U`, or `O` for an old form, whose eigenvalue is `NULL` if T_p is not an element of its field on its space (then the number would depend on the vector which is kept for the form). T_p is normalized by p^a, which makes its eigenvalues algebraic integers: on the lattice of discriminant 83, the first eigenvalues with the weight (1, 0) are T_2 = -17 and T_3 = -23, which are p times those of Magma. T_{p^2}, the operator of the p^2-neighbors, is normalized by p^(2a), p^(2a) times `HeckeOperator(M, p, 2)` in Magma. Its eigenvalues are then algebraic integers, and the factor at p of the spin L-function of a form is, with w = 2a + 3,

    1 - T_p x + p (T_{p^2} + p^(2a) (p^2 + 1)) x^2 - p^w T_p x^3 + p^(2w) x^4,

as for the weight 0. With this formula every root has absolute value p^(w/2), for p = 2, 3 and 5, for all the forms of the weights (1, 1), (2, 1), (2, 2), (3, 1) and (3, 2) on the lattices of discriminant 19, 61 and 83. For example, the paramodular form 2.K.19.3.2.a.a, of weight (1, 1), has T_2 = -9 and T_4 = 0, and the factor 1 + 9x + 40x^2 + 288x^3 + 1024x^4. T_{p^2} tests every neighbor for isometry, not one in each orbit of the automorphisms as T_p does. A line on the standard error gives the dimension of the space of each conductor.

With a weight, -p=p may be a prime dividing the discriminant. The Hecke matrices are then those of the p-neighbors (the p^2-neighbors for T_{p^2}), with the characters of the conductors as without the weight (see below). T_p takes the neighbors by orbits of the automorphisms, as at the other primes, and T_{p^2} tests every neighbor for isometry. At the conductor 1 they are p^a (p^(2a) for T_{p^2}) times those of `HeckeOperator(M, p, k)` in Magma, for the weights (1, 1) and (2, 1) on the lattices of discriminant 26 and 39 at each of their primes, and of discriminant 45, 50, 52, 63, 68, 75 and 76 at the prime whose square divides it, and for (1, 0) and (1, 1) on the lattice of discriminant 83 at p = 83. The eigenvalues are printed as without the weight, with p^a in place of 1, and w = 2a + 3:
  - When p divides the discriminant exactly once, the eigenvalue of T_p is printed as lambda = T_p + p^a, negated when p divides the conductor. The factor at p of the L-function of a newform is then (1 - (lambda + eps p^(a+1)) x + p^w x^2) (1 + eps p^(a+1) x) for a sign eps, the shape of Roberts and Schmidt for the paramodular level p. This was checked against the functional equation for the newforms with rational eigenvalues of the weight (1, 1) on the lattices of discriminant 15, 19, 26, 39, 45, 50, 55 and 75, and of the weight (2, 1) on those of discriminant 13 and 39. The sign eps is the Atkin-Lehner sign at p, as in `parse_omf5` of scripts/parse_omf.py (Theorem 10.1 of Dummigan-Pacetti-Rama-Tornaria): 1 when p does not divide the conductor and -1 when it does, negated when p ramifies in the definite quaternion algebra of the genus. This held at each prime of the discriminants 15, 26, 39, 45, 50, 55 and 75 where the two signs could be told apart (3 and 5 for 15, 2 for 26, 3 for 39, 5 for 45 and 55, 2 for 50, 3 for 75), and among these 5 for 15 and 45, 2 for 50 and 3 for 75 are ramified. It also holds with b = 0 and for the other weights at prime discriminants N, where N is always ramified, through the root number (-1)^k eps_N of the functional equation, k = a - b + 3, for the newforms with rational eigenvalues (or, for 59, both embeddings of a quadratic field) which are not lifts: with -prec=200 the root number is +1 for the weight (1, 0) at 31, 37, 41, 43, 47 and 59 (conductor N, eps_N = 1), +1 for the weight (1, 1) at 19, 23 and 29 (conductor 1, eps_N = -1) and -1 at 43 (conductor 43, eps_N = 1), and +1 for the weight (2, 1) at 11 and 13 (conductor N, eps_N = 1). For N at least 17, N^2 is larger than the 199 coefficients used, so these determine the root number and not eps in the factor at N. For 2.K.19.3.2.a.a, T_19 + 19 = -969.
  - A form which is old at p has the Euler factors of the form of lower level it comes from. For example the forms of conductor 5 and 15 of the weight (1, 1) on the lattice of discriminant 45 with T_2 = -14 come from the form of conductor 15 of the discriminant 15: they have its eigenvalues at the primes other than 3, and T_3 = 15 - 27 and 15 + 27, lambda_3 - theta 3^(a+2), with lambda_3 = 15 that of the form of discriminant 15 and theta = 1 for the conductor 5 and -1 for 15.
  - When p^2 divides the discriminant, the eigenvalue of T_p is printed as it is, and for the newforms with rational eigenvalues of the weight (1, 1) on the lattices of discriminant 45, 50 and 75 the factor at p was 1 - T_p x + (p^w + p T_{p^2} / 2) x^2, as for the weight 0 (with T_{p^2} not 0 for four of them).
  - The eigenvalue of T_{p^2} is printed as it is.

The prime p may divide the discriminant. The Hecke matrix is then that of the p-neighbors, and the character of a conductor divisible by p is found at the isometries to the p-neighbors, which are not integral at p, as follows. When p divides the discriminant exactly once, it is given by the parity of the valuation at p of the spinor norm. When p^2 divides the discriminant and the radical of the lattice modulo p is a line, the lattice and its p-neighbors have the same radical modulo p, the isometry acts on it by a sign, and the character is this sign. The Hecke matrices are then those of HeckeOperator in the package ModFrmAlg of Magma.

The radical of the lattice modulo a prime q dividing the discriminant need not be a line. When q is odd, the character at q of an isometry which is integral at q is the determinant by which it acts on the radical, after the isometry is multiplied by -1 if its determinant is -1. The radical of a q-neighbor is identified with that of the lattice, and the character at an isometry to a q-neighbor is the determinant of the map it induces from the radical to itself. It is 1 at an isometry which is the identity on the orthogonal complement of a hyperbolic plane, as it is when the radical is a line.

At the prime 2 the radical is that of the bilinear form modulo 2, and its dimension is 1, 3 or 5.

- Dimension 3. The character depends on whether the quadratic form is even on the radical.
  - The quadratic form is even on the radical. An isometry which is integral at 2 acts on the quotient of the lattice modulo 2 by the radical, a plane over the field of 2 elements, and its character at 2 is the sign of the permutation of the three nonzero vectors of this plane. This is (-1)^rank(g - 1) for its action g on the plane. When the radical is a line, the same formula for the action on the quotient of the lattice modulo 2 by the radical gives the sign on the radical. When the lattice has 2-neighbors, the plane has two isotropic vectors other than 0. An isometry takes the lattice to its 2-neighbor at x, and its inverse takes it to the 2-neighbor at x'. The character is -1 when x and x' are the same in the plane, and 1 otherwise. It is 1 at an isometry which is the identity on the orthogonal complement of a hyperbolic plane.
  - The quadratic form is odd at a vector of the radical, and 16 does not divide the discriminant. The character at an isometry which is integral at 2 is given by the parity of the valuation at 2 of the spinor norm, as when 2 divides the discriminant exactly once: it is 1 when the valuation is even. At an isometry to a 2-neighbor it is the opposite, so that it is 1 at an isometry which is the identity on the orthogonal complement of a hyperbolic plane. When 4 divides the discriminant exactly, the vectors of the radical at which the quadratic form is even are a plane over the field of 2 elements, and in the lattices in which it was compared the character at an integral isometry is the sign of the permutation of the three nonzero vectors of this plane. omf5 computes this sign there, as (-1)^rank(g - 1) for the action g of the isometry on the radical modulo 2, because it is much cheaper than the spinor norm.
  - The quadratic form is odd at a vector of the radical, and 16 divides the discriminant. The vectors whose products with the lattice are divisible by 2^k are the radical modulo 2^k, and the deepest part of the radical is its image in the lattice modulo 2 for the largest k at which this image is not 0. It is the last component of a Jordan decomposition of the lattice at 2, modulo 2, and here it is a line or a plane. The character at 2 is the sign by which an isometry acts on it. When it is a line spanned by e, an isometry which is integral at 2 takes e to e or to -e modulo 4 times the lattice, and the character is this sign, after the isometry is multiplied by -1 if its determinant is -1. When it is a plane, the character is the sign of the permutation of its three nonzero vectors. A 2-neighbor of the lattice has the same deepest part, modulo 4 times either lattice when it is a line, and the same formula gives the character at an isometry to a 2-neighbor. It is 1 at an isometry which is the identity on the orthogonal complement of a hyperbolic plane. The character of the spinor norm is not used here: with it, the space of a conductor 2d of some of these lattices is the same as that of the conductor d, with the same Hecke matrices.

  A product of an even number of reflections which are integral at 2 is an isometry of determinant 1. In the lattices in which it was compared, each of these two characters at such a product is given by the parity of the number of its reflections in vectors v at which the valuation at 2 of Q(v) is the largest that a reflection integral at 2 has in this lattice.

  In all these cases the lattice has no 4-neighbors, and T_4 is 0.
- Dimension 5. The bilinear form is 0 modulo 2, and there is no character at 2: the conductors are then odd, and a warning says so.

A lattice has no p-neighbors when no vector outside its radical modulo p is isotropic modulo p. The Hecke matrix at p is then 0.

The eigenvalue at a prime p which divides the discriminant exactly once is printed as T_p + 1, and as -(T_p + 1) when p divides the conductor, where T_p is the eigenvalue of the Hecke matrix. When p^2 divides the discriminant it is printed as T_p.

The eigenvalue of T_{p^2} at a prime p dividing the discriminant is printed as it is. The lattice has no p^2-neighbors when the quotient of the lattice modulo p by its radical has no totally isotropic plane, and T_{p^2} is then 0. This is so when the radical is not a line, and it can be so when p divides the discriminant exactly once. When p^2 divides the discriminant and the radical modulo p is a line, the factor at p of the L-function of a form of type G was 1 - T_p x + (p^3 + p T_{p^2} / 2) x^2 in every example checked against the functional equation. In the notation of Roberts and Schmidt for the paramodular group, T_p is lambda and T_{p^2} / 2 is mu. A guess for the factor of a form which is new at p when the radical is a hyperbolic plane, and what supports it, is in doc/euler-factors-bad-primes.md.

If the flag -isom is supplied in addition to -genus=g, read from the genus file, in addition to the list of genera, a list of fixed isometries (over QQ) between the different lattices in the genus.

The genus file is a list of genera, and the genus read is the one at index d, counting from 0. A genus is a list of lattices, and may be an empty list if d is not 0. Each lattice is a list of 15 integers in the format 'GG', as in this file with a genus of two lattices at index 1

[[],[[1,0,-1,-1,1,1,-1,1,0,2,1,-1,2,-2,2],[1,-1,0,1,1,1,0,-1,0,1,1,0,2,1,3]]]

With -isom, each lattice is instead a list of two 5x5 matrices, given by their rows - the Gram matrix of the lattice, and an isometry from the first lattice of the genus to it, whose entries are rational numbers such as -2/3.

The lattices have to be a full set of representatives for the classes of a single genus. The program exits with an error if one of them is not positive definite and even, if they do not all have the same discriminant, if two of them are isometric, if their mass is not the mass of the genus, or if one of them is not in the genus of the first lattice. With -isom the last of these is not checked, and each isometry is checked to be an isometry from the first lattice, with a denominator which is prime to the discriminant. The isometries are multiplied by one another in the computation, as fractions of 64 bit integers, and the program exits with an error if a product in lowest terms, or the inverse of an isometry, does not fit in them.

If the flag -tests is supplied, additionally runs standard tests.

If the flag -nonlifts is supplied, computes eigenvalues only for the forms which are not lifts. It needs -prec=L, which determines the types of the forms, unless form_idxs is supplied.

Each form is printed with its automorphic type, `aut_rep_type`, which is one of:

- `G` - the form is not a lift,
- `P` - a Saito-Kurokawa lift,
- `Y` - a Yoshida lift,
- `F` - the Eisenstein series, which is the constant function,
- `O` - an old form, for which no eigenvector was found. The vector used for it is an eigenvector of one Hecke operator, in a space which the Hecke operators at the primes up to 20 not dividing the discriminant do not split. The eigenvalues printed for it are those of an eigenform only if the vector is an eigenvector of the other operators as well, which is not checked at the primes not dividing the discriminant. It need not be an eigenvector of the operator at a prime p dividing the discriminant, and the number found from it then depends on the vector, which depends on the basis of the lattice. So with -p=p the program computes the matrix of the operator, and prints the eigenvalue of this form only when the operator acts on that space as an element of the field of the form: every vector which could have been used is then an eigenvector, with this eigenvalue. Otherwise no eigenvalue at p is printed for this form, also when form_idxs gives the form the type `G`, and a warning says so. With -square the same is done for T_{p^2}. The matrix takes longer to compute than an eigenvalue does.
- `U` - not determined, since -prec was not supplied (see below).

With -prec=L, a form is taken to be a lift if its L-polynomial factors, in the same way, at each of the primes p not dividing the discriminant with p^2 <= L. If there are fewer than three such primes, the first three primes not dividing the discriminant are used. (With two, forms which are not lifts were reported as Yoshida lifts, for example on the lattices of discriminant 50 and 63 with the weight (1, 1); among the forms of discriminant at most 83 computed with -prec=200, no form which is not a lift has an L-polynomial which factors as that of a lift at three of the first primes not dividing the discriminant.) So the types are more reliable when L is larger. Even then the types `P`, `Y` and `F` say only that the L-polynomial factors in this way at the primes checked: a form which is not a lift can have an L-polynomial which factors at all of them. The function `aut_rep` in scripts/parse_omf.py confirms a type: it looks in the LMFDB for the classical modular forms of which the form would be the lift (for `Y`, newforms f of weight 2 and g of weight 4, of trivial character and levels dividing the discriminant, with T_p = p a_p(f) + a_p(g) at the primes p below a bound, through `check_yoshida`), and returns them. Without -prec, as with -p=p alone, the program computes only what was asked for, and not the eigenvalues at other primes which the types need: the forms which are not old get the type `U`, a warning says so, and -nonlifts, which needs the types, is an error unless -idxs is supplied. For example, for discriminant 249 there is a form which is not a lift, whose L-polynomial factors at 2 and at 5, but not at 7. With two primes it was reported as a Yoshida lift for L below 49, and so left out when -nonlifts was supplied; with three, 7 is always used and it is `G`.

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
Warning: without -prec the types of the forms are not determined: they are printed as U, and as O for the old forms. Use -prec=L to determine them.
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
Found form of type U at total_idx = 0, c = 0, form_idx = 0
Found form of type U at total_idx = 1, c = 0, form_idx = 1
Found form of type U at total_idx = 2, c = 0, form_idx = 2
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
{1 : [{ 'aut_rep_type' : 'U', 'field_poly' : '2  7 1', 'lambda_p' : [], 'lambda_p_square' : [], }, { 'aut_rep_type' : 'U', 'field_poly' : '2  -15 1', 'lambda_p' : [], 'lambda_p_square' : [], }, { 'aut_rep_type' : 'U', 'field_poly' : '7  2026 -5205 4471 -1714 322 -29 1', 'lambda_p' : [], 'lambda_p_square' : [], }, ], 61 : [], }
```

Without -prec the types of the forms are not determined, so all three have the type `U`. With -prec=100, as in the next example, they are `G`, `F` and `P`: the first form is not a lift, the second is the Eisenstein series, and the third is a Galois orbit of 6 Saito-Kurokawa lifts.

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
