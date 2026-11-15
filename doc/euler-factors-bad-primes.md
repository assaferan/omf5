# Euler factors at the primes dividing the discriminant

This note records which factors of the L-function were found at a prime p
dividing the discriminant D, and how well each is supported. They are the
factors of a form which is new at p. In the last case, a radical modulo p
which is a hyperbolic plane, the factor is a guess.

T_p and T_{p^2} are the eigenvalues of the Hecke matrices of the p-neighbors
and of the p^2-neighbors (`-hecke -p=p`, `-hecke -p=p -square`), d is the
conductor of the form, and the L-function is the one whose factor at a prime
not dividing D is

    1 - T_p x + p (T_{p^2} + p^2 + 1) x^2 - p^3 T_p x^3 + p^6 x^4,

with the functional equation s -> 4 - s. The factors below are those of
Section 2 of Johnson-Leung and Roberts (arXiv:1006.5105) for a paramodular
newform of weight 3, which follow Roberts and Schmidt. There lambda and mu
are the eigenvalues of the two Hecke operators at p, and eps is that of the
Atkin-Lehner involution at p. The sign of the functional equation is minus
the product of the eps at the primes dividing the conductor.

A form is called new at p here when the exponent of p in the conductor of
its L-function is the exponent of p in D. The program does not tell which
forms are new at p. The L-function of a lift of type Y is
L(s - 1, f) L(s, g) for newforms f of weight 2 and g of weight 4, and its
conductor is the product of their levels. Below, f and g were found with
PARI from the eigenvalues at the primes not dividing D, where T_p is
p a_p(f) + a_p(g).

## What is known

1. p divides D exactly once, and the form is new at p. Let c be the value
   printed by `-p=p`, which is T_p + 1, and -(T_p + 1) when p divides d. The
   factor is

       (1 - (c + eps p) x + p^3 x^2) (1 + eps p x),

   that of a newform of level p, with c = lambda + eps and
   mu = -eps c - p. In the examples below T_{p^2} was 0 for every form of
   the lattice or for none. At a conductor prime to p, eps was -1 in the
   first case and 1 in the second, and at a conductor divisible by p it was
   the opposite. For the lifts this is read from f and g (D = 150 at 2 and
   3, with d = 1 and d = 30, D = 175 at 7, D = 245 at 5, D = 275 at 11,
   D = 363 at 3, with d = 1 and d = 3), and for the form of type G of
   D = 275 from the functional equation.

   A form which is old at p has another factor. At D = 363 and p = 3 the
   form of type Y with the field x^2 + 6 x + 6 is the lift of two newforms
   of level 11. The conductor of its L-function is 121, and its factor at 3
   has degree 4. With lambda the coefficient of -x in this factor, the
   values printed for this form are lambda + 9 at d = 1 and lambda - 9 at
   d = 3.
2. p^2 divides D and the radical modulo p is a line. For a form of type G
   which is new at p the factor is 1 - T_p x + (p^3 + p T_{p^2} / 2) x^2,
   that of a newform of level p^n, n >= 2, with lambda = T_p and
   mu = T_{p^2} / 2. It was checked against the functional equation. For
   the six lifts new at p of the cyclic comparison below (D = 150, 175,
   245, 275 and 363) it is the factor of L(s - 1, f) L(s, g).

   A form which is old at p has another factor here also. On the lattice of
   D = 275 with a cyclic discriminant group, the lift of the newforms of
   levels 11 and 5 has T_5 = 25 and T_25 = 100 at d = 1. Its factor at 5 is
   (1 - 5 x + 125 x^2) (1 + 5 x), and not 1 - 25 x + 375 x^2.

## The radical modulo p is a plane

Let p be odd, and let the discriminant group at p be (Z/p)^2. The radical
modulo p is a plane, p^2 is the exact power of p dividing D, and
T_{p^2} = 0. The quadratic form divided by p is a quadratic form on the
radical with values modulo p, and the radical is a hyperbolic plane or an
anisotropic one for it.

### A hyperbolic radical

For a form which is new at p the guess is that the factor at p is

    1 - (T_p + e (p - 1)) x + p^2 x^2,

where e = 1 when p does not divide d and e = -1 when it does, and that
eps = 1. This is the factor of case 2 with lambda = T_p + e (p - 1) and
mu = p - p^2. In the table of Hecke eigenvalues in Section 4 of
Johnson-Leung and Roberts, the rows of level p^2 with this mu are those of
the representations of type IIIa and VIa, both with eps = 1. In the
examples below eps = 1 is a_p(f) a_p(g) / p for the lifts, and for the form
which is not a lift it follows from the sign of the functional equation.

The sign e depends on the sign which was chosen for the character at p (see
the README). With the other choice T_p is negated at the conductors divisible
by p, and the guess reads T_p + p - 1, negated when p divides d, as T_p + 1 is
in case 1.

A form which is old at p does not have this factor: see the form of type Y of
D = 275 with T_5 = 21 below. Nothing is guessed for the lifts of type P and
F, or for a form of type O.

### An anisotropic radical

No factor is guessed. On the three such lattices of the comparison below,
D = 50, 147 and 363, T_p is -(p + 1) for every form but the constant
function and one form of type P at D = 147 and at D = 363, so that T_p does
not tell the forms apart. The five forms of type Y at D = 363, which are
four eigensystems, are lifts with p a_p(f) = p = -a_p(g), so that their
factor at 11 is 1 - 121 x^2. This is the factor of case 2 with
lambda = 0 and mu = -p^2 - p, the row of type Va of the same table, which
has eps = -1. The formula of the hyperbolic case would give
1 + 2 x + 121 x^2.

### Evidence

The lattices, as given to `-quad` with `-format=GG`. The radical modulo p is
hyperbolic for the first four and anisotropic for the fifth, and the last has
a cyclic discriminant group:

    D = 150, p = 5:   1,0,-1,0,-1,1,-1,1,0,2,0,1,3,-2,3
    D = 175, p = 5:   1,0,0,0,1,1,-1,1,-1,2,-2,0,2,-1,5
    D = 245, p = 7:   1,1,1,-1,-1,1,1,0,0,3,-2,1,3,2,4
    D = 275, p = 5:   1,-1,0,-1,-1,1,0,0,0,2,-1,-1,4,3,4
    D = 363, p = 11:  2,1,1,-1,-2,2,-1,1,-2,2,-1,0,2,0,3
    D = 275, cyclic:  1,0,-1,-1,0,1,1,0,-1,3,2,0,3,0,3

*Lifts.* The forms of type Y on the first four lattices which are new at p.
The column "primes" is the number of primes below 250 not dividing D at
which T_p = p a_p(f) + a_p(g) was checked.

| D | d | T_p | level of f | level of g | primes | factor of the lift at p |
|---|---|-----|------------|------------|--------|-------------------------|
| 150 | 1 | -14 | 30 | 5 | 50 | (1 + 5 x)^2 |
| 150 | 30 | 14 | 15 | 10 | 50 | (1 - 5 x)^2 |
| 175 | 1 | -14 | 35 | 5 | 51 | (1 + 5 x)^2 |
| 245 | 1 | -20 | 35 | 7 | 51 | (1 + 7 x)^2 |
| 275 | 1 | -14 | 55 | 5 | 51 | (1 + 5 x)^2 |

In every row p divides each of the two levels exactly once and
p a_p(f) = a_p(g), and the factor (1 - a_p(g) x)^2 of the lift is the
guess. The forms of D = 245 and D = 275 have degree 2, as their f has, and
both embeddings were checked. For D = 175 the first 250 coefficients found
from the eigenvalues, with the guess at 5 and the factor of case 1 with
eps = -1 at 7, are those of L(s - 1, f) L(s, g).

*Functional equation, for the form which is not a lift.* The form of type G
of D = 275 has T_5 = -4 at d = 1 and T_5 = 4 at d = 5. These are one
L-function: the eigenvalues at the other primes are the same, and the
guessed factor is 1 + 25 x^2 for both. `lfuncheckfeq` of PARI with 250
coefficients, the gamma factors [-1, 0, 0, 1], weight 4, conductor 275,
sign 1, and the factor of case 1 at 11 (c = 33). The value is the logarithm
in base 2 of the error in the functional equation, so that -30 is an
agreement to 30 bits.

| bits | guess, eps = -1 at 11 | guess, eps = 1 at 11 | factor 1 | 1 + 125 x^2 |
|------|-----------------------|----------------------|----------|-------------|
| 40 | -29 | -29 | -16 | -14 |
| 52 | -46 | -39 | -16 | -14 |
| 64 | -59 | -39 | -16 | -14 |

With the sign -1 the guess gives 0 or 1 at the three precisions, for both
eps at 11. So the sign is 1, and eps at 11 is -1, as in case 1 when T_{p^2}
is 0. The sign being minus the product of the two eps, eps is 1 at 5. The
two eps at 11 differ in the coefficients of 121 and 242 only, and are told
apart from 52 bits on.

At 40 bits the other candidates were 1 - 25 x^2, 1 -+ 4 x and
1 -+ 4 x + 125 x^2 (-6 to -15), and each candidate and the guess with the
conductors 11, 55, 1375 and 6875 (0 to -3).

*Functional equation, for the lifts.* The same check at 40 bits, where the
value does not go below -30, with the conductor D and the factors of case 1
at the other primes dividing D:

| D | d | guessed factor | guess | best other candidate |
|---|---|----------------|-------|----------------------|
| 150 | 1 | 1 + 10 x + 25 x^2 | -30 | -16 |
| 150 | 30 | 1 - 10 x + 25 x^2 | -30 | -17 |
| 175 | 1 | 1 + 10 x + 25 x^2 | -30 | -16 |
| 245 | 1 | 1 + 14 x + 49 x^2 | -28, -29 | -22, -23 |

The other candidates were the guess with -e in place of e, with T_p in place
of T_p + e (p - 1), with p^3 in place of p^2, the factors 1 and 1 - T_p x,
each with the conductors D, D / p and D p, and with both signs eps at the
other primes. The best of them is the guess with p^3 in place of p^2 in
every row. At D = 175 the check at 40 bits does not tell the two signs eps
at 7 apart (-30 and -25); the lift gives eps = -1.

*A lattice outside the comparison.* The lattice A_3 + 7 A_2 of D = 294 has a
hyperbolic radical modulo 7:

    D = 294, p = 7:   1,1,0,0,0,1,1,0,0,1,0,0,7,7,7

Its forms of type Y, at the conductors 1, 2, 3 and 6, are lifts, found from
the 22 primes below 100 not dividing D:

| d | T_7 | level of f | level of g | factor of the lift at 7 |
|---|-----|------------|------------|-------------------------|
| 1, 2 | -20 | 21 | 7 | (1 + 7 x)^2 |
| 1 | -20 | 21 | 14 | (1 + 7 x)^2 |
| 6 | -20 | 42 | 7 | (1 + 7 x)^2 |
| 1 | 8 | 14 | 21 | (1 - 7 x)^2 |
| 3 | 8 | 14 | 21 | (1 - 7 x)^2 |

All are new at 7, and the factor is the guess in every row. The form of the
first row is old at 2 and is found at two conductors, and that of the fourth
has degree 2, as its g has. In the last two rows T_p + e (p - 1) is 2p with
e = 1; in the tables above it is 2p only at d = 30, where e is -1.

*A form which is old at p.* The form of type Y of D = 275 at d = 1 with the
field x + 8 has T_5 = 21. It is the lift of the newforms of levels 11 and 5
(51 primes). The conductor of its L-function is 55 and its factor at 5 is
(1 - 5 x + 125 x^2) (1 + 5 x). The functional equation gives -29 for this
factor with the conductor 55, and -10 for 1 - 25 x + 25 x^2 with the
conductor 275.

*The anisotropic radical.* The five forms of type Y of D = 363, three at
d = 1 and two at d = 3, have T_11 = -12. The one with the field
x^2 + 6 x + 6 is at both conductors, so that they are four eigensystems.
They are lifts of newforms whose levels are 11 or 33, with a_11(f) = 1 and
a_11(g) = -11 (23 primes below 100). For the one with the field x + 9 at
d = 3, with 100 coefficients, the conductor 363 and the sign -1, the
functional equation gives -24 for 1 - 121 x^2 and -7 for 1 + 2 x + 121 x^2.
With 100 coefficients this tests the coefficient of x only.

*The lattice with a cyclic discriminant group.* For 15 lattices with D <= 400
and discriminant group Z/m x Z/p the eigenvalues were compared with those of a
lattice of the same D whose radical modulo p is a line, the forms being
matched by their eigenvalues at three primes not dividing D. They are D = 18,
50, 126, 147, 150, 175, 245, 252, 275, 300 and 363, where p^2 is the exact
power of p dividing D, and D = 135, 162, 189 and 375. Two of them, D = 252 and
300, could not be compared. The lattices with a cyclic discriminant group
whose forms are used in this note, besides that of D = 275 above:

    D = 150:  1,0,1,0,0,1,-1,-1,-1,2,0,1,2,1,4
    D = 175:  1,1,1,0,-1,2,0,1,-2,2,0,-2,2,-1,3
    D = 245:  1,0,0,0,0,1,0,-1,1,1,0,1,4,2,5
    D = 363:  1,0,0,0,1,2,0,1,0,2,-2,1,3,1,3

On the 9 lattices where p^2 is the exact power of p dividing D:

- 9 eigensystems match a form with T_{p^2} = 2 (p - p^2), at D = 150, 175,
  245 and 275. For all of them T_p of the cyclic lattice is T_p + e (p - 1),
  so that the guess is the factor of case 2 for the form of the cyclic
  lattice. These are the lifts of the first table and the form of type G
  of D = 275.
- 7 eigensystems, six of type P and one of type Y (D = 275, T_p = 21), match
  a form with another T_{p^2}. T_p differs by p - 1 for them also, and the
  formula of case 2 on the cyclic lattice does not give the guess. The one
  of type Y is the form which is old at p, and it also matches a second
  form of the cyclic lattice.
- The constant function has T_p smaller by p (p + 1) on all 9 lattices.
- 4 eigensystems at D = 147 and 363, where the radical is anisotropic, match
  a form of the cyclic lattice whose conductor is divisible by p, and T_p
  differs by p + 1. Two of them are a form of type Y of D = 363, which
  matches a form with T_p = 0 and T_{p^2} = -2 (p^2 + p): the factor of
  case 2 for it is 1 - 121 x^2.
- 4 eigensystems at D = 50, 147 and 363 have T_p = -(p + 1), against 2 p^2
  on the cyclic lattice.
- 27 eigensystems are not found on the cyclic lattice at all.

### What is not known

- In every case where the guess was tested, T_p + e (p - 1) is 0, 2p or -2p.
  The shift p - 1 was not tested on any other value, and it was not derived.
  A lift does not test it on another value when p divides each of the two
  levels exactly once, as it does for every lift new at p which was found:
  a_p(f) is then 1 or -1 and a_p(g) is p or -p, so that p a_p(f) + a_p(g)
  is 0, 2p or -2p. A form of type G for which it is not 0 would.
- The form of type G of D = 275 is the only one which is not a lift, and the
  functional equation is all that supports the guess for it.
- Nothing is guessed for the lifts of type P and F, for the forms which are
  old at p, or when the radical is anisotropic.
- When p^3 divides D (D = 135, 162, 189, 375) T_p differs from that of the
  cyclic lattice by p, or -p at a conductor divisible by p, in 7 matches,
  and by p (p + 1) for the constant function. Two more forms, at D = 375,
  are not told apart by three primes. No factor was tested.
- There is no proof. ModFrmAlg of Magma does not find these spaces, and the
  only comparison with another program is that of the lifts with the
  classical newforms of PARI.
