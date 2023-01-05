import pickle
from sage.all import (PolynomialRing, QQ, NumberField, prime_range, prime_divisors)
from sage.misc.persist import SagePickler

def parse_omf5(k,j,N,folder,suffix="_nl_200_",hecke_ring=False,max_deg=20,B=200):
    fname = folder + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    Qx = PolynomialRing(QQ, name="x")
    x = Qx.gens()[0]
    al_signs = eval(open(fname).read())
    ret = []
    for al_sign in al_signs:
        forms = al_signs[al_sign]
        for f in forms:
            pol = Qx([int(c) for c in f['field_poly'].split()[1:]])
            F = NumberField(pol, name = "a")
            a = F.gens()[0]
            f['lambda_p'] = [F(lamda) for lamda in f['lambda_p']]
            f['lambda_p_square'] = [F(lamda) for lamda in f['lambda_p_square']]
            is_lft, lift_type = is_lift(f, N, B)
            if (is_lft):
                f['aut_rep_type'] = lift_type
            if (F.degree() > max_deg):
                f['trace_lambda_p'] = [lamda.trace() for lamda in f['lambda_p']]
                f['lambda_p_square'] = [lamda.trace() for lamda in f['lambda_p_square']]
                f['lambda_p'] = []
                f['lambda_p_square'] = []
            # !! TODO : represent the eigenvalues in the polredabs field, currently some things break
            coeffs = [int(c) for c in pol.coefficients(sparse=False)]
            f['field_poly'] = coeffs
            f['atkin_lehner_eigenvals'] = [[p, -1 if al_sign % p == 0 else 1] for p in prime_divisors(N)]
            f['atkin_lehner_string'] = ''.join(['-' if al_sign % p == 0 else '+' for p in prime_divisors(N)])
            if (hecke_ring):
                # !!! This can be very slow
                # hecke_ring = F.order(orbit['lambda_p'])
                index = 0
                p_idx = 0
                nbound = 0
                while ((index != 1) and (p_idx < len(f['lambda_p']))):
                    H = F.order(f['lambda_p'][:p_idx+1])
                    new_index = H.index_in(F.ring_of_integers())
                    if (new_index != index):
                        index = new_index
                        nbound = p_idx
                        p_idx += 1
                f['hecke_ring'] = H
                f['hecke_ring_index'] = index
                f['hecke_ring_generator_nbound'] = nbound
                        
            ret.append(f)
    return ret

def pickle_omf5(k,j,N,folder,suffix="_nl_200_"):
    fname = folder + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    parsed = parse_omf5(k,j,N,folder)
    f = open(fname, "wb")
    gherkin = SagePickler.dumps(parsed)
    _ = f.write(gherkin)
    f.close()
    # check
    f = open(fname, "rb")
    pickled = f.read()
    assert parsed == pickle.loads(pickled)
    return

def unpickle_omf5(k,j,N,folder,suffix="_"):
    fname = folder + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    f = open(fname, "rb")
    pickled = f.read()
    return pickle.loads(pickled)

def is_newform(f, N, G_type, prime_bound):
    divs = divisors(N)[:-1]
    primes_N = [p for p in prime_range(prime_bound) if N % p != 0]
    if (len(f['lambda_p']) > 0):
        ap = [x.trace() for x in f['lambda_p']]
    else:
        ap = f['trace_lambda_p']
    for d in divs:
        primes_d = [p for p in prime_range(prime_bound) if d % p != 0]
        inds = [i for i in range(len(primes_d)) if primes_d[i] in primes_N]
        for g in G_type[d]:
            if (len(g['lambda_p']) > 0):
                bp = [g['lambda_p'][i] for i in inds]
                bp = [x.trace() for x in bp]
            else:
                bp = [g['trace_lambda_p'][i] for i in inds]
            if (ap == bp):
                return False
    return True

def is_lift(f, N, prime_bound):
    primes_N = [p for p in prime_range(prime_bound) if N % p != 0]
    if (len(f['lambda_p_square']) > 1):
        K = f['lambda_p'][0].parent()
        Kx = PolynomialRing(K, name = "x")
        x = Kx.gens()[0]
        lps = [p^6 * x^4 - f['lambda_p'][i]*p^3*x^3 + (f['lambda_p_square'][i] + p^2 + 1)*p*x^2 - f['lambda_p'][i] * x + 1 for i,p in enumerate(primes_N[:len(f['lambda_p_square'])])]
        if any([lp.is_irreducible() for lp in lps]):
            return False, 'G'
        if (len(lps[0].factor()) == 3):
            lift_type = 'P'
        else:
            lift_type = 'Y'
        return True, lift_type
    # Not final, should check these as well
    # There are those where we saved traces, but also multiples of 210
    return False, 'G'

def load_ev_data(k,j,folder,suffix="_"):
    ev_data = [unpickle_omf5(k,j,N,folder,suffix) for N in range(1,1000)]
    ev_data = [[]] + ev_data
    G_type = [[f for f in ev_data[N] if f['aut_rep_type'] == 'G'] for N in range(1000)]
    return ev_data, G_type

def update_ev_data(ev_data, B):
    G_type = [[f for f in ev_data[N] if f['aut_rep_type'] == 'G'] for N in range(1000)]
    for N in range(1000):
        print(N)
        for f in G_type[N]:
            is_new = is_newform(f, N, G_type, B)
            if (not is_new):
                f['aut_rep_type'] = 'O'
            is_lft, lift_type = is_lift(f, N, B)
            if (is_lft):
                f['aut_rep_type'] = lift_type
    return
    
