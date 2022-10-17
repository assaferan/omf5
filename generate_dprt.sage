load("dubey_holenstein.sage")
def generate_form_dprt(dminus, dplus):
    d = dminus*dplus
    F = factor(d)
    primes = [ p[0] for p in F]
    assert moebius(dminus) == -1, "μ(d-) must be -1"
    assert gcd(dminus, dplus) == 1, "d-, d+ must be relatively prime"
    genus = list()
    if 2 not in primes:
        primes.append(2)
    for p in primes:
        n = d.valuation(p)
        R = IntegerModRing(p^(n+4))
        M = Matrix(R, 5, 5, 0)
        if dminus % p != 0:
            M[0,0] = 2*d
            M[1,2] = d/(p^n)
            M[2,1] = M[1,2]
            M[3,4] = d/(p^n)
            M[4,3] = M[3,4]
        else:
            if p != 2:
                r = distinguished_nonsquare(p)
                M[0,0] = 2*d*r
                M[1,1]=2*d/p
                M[2,2]=-2*d*r/p
                M[3,4]=d/p
                M[4,3]=d/p
            else:
                M[0,0] = -2*d*R(1/3)
                M[1,2] = -d/2
                M[2, 1] = -d/2
                M[1, 1] = -d
                M[2, 2] = -d
                M[3, 4] = d/2
                M[4, 3] = d/2
        genus.append((p, normal_form(M, p^(n+4))[1]))
    form = QuadraticForm(form_from_genus(genus))
    return form
