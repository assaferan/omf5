import pickle
from sage.all import (PolynomialRing, QQ, NumberField, prime_range, prime_divisors, divisors, is_square, ZZ)
from sage.misc.persist import SagePickler

def parse_omf5(k,j,N,folder,suffix="_nl_200_",hecke_ring=True,max_deg=20,B=200):
    fname = folder + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    Qx = PolynomialRing(QQ, name="x")
    x = Qx.gens()[0]
    al_signs = eval(open(fname).read())
    ret = []
    for al_sign in al_signs:
        forms = al_signs[al_sign]
        if (len(forms) > 0):
            f = forms[0]
            bad_ps = [p for p in primes_first_n(len(f['lambda_p'])) if N % p == 0]
            ps = primes_first_n(len(f['lambda_p']) + len(bad_ps))
            good_ps = [p for p in ps if p not in bad_ps]
            assert len(good_ps) == len(f['lambda_p'])
            
        for f in forms:
            pol = Qx([int(c) for c in f['field_poly'].split()[1:]])
            F = NumberField(pol, name = "a")
            a = F.gens()[0]
            f['lambda_p'] = ['NULL' if ps[i] in bad_ps else F(f['lambda_p'][good_ps.index(ps[i])]) for i in range(len(ps))]
            f['lambda_p_square'] = ['NULL' if ps[i] in bad_ps else F(f['lambda_p_square'][good_ps.index(ps[i])])
                                    for i in range(len(f['lambda_p_square']))]
            # At the moment can't do it on toby - needs to set up LMFDB there
            # aut_tp, friends = aut_rep(f,N,B,db)
            # f['aut_rep_type'] = aut_tp
            # if aut_tp in ['Y', 'P']:
            #   f['related_objects'] = friends
            if (F.degree() > max_deg):
                f['trace_lambda_p'] = ['NULL' if ps[i] in bad_ps else f['lambda_p'][good_ps.index(ps[i])].trace() for i in range(len(ps))]
                f['trace_lambda_p_square'] = ['NULL' if ps[i] in bad_ps else f['lambda_p_square'][good_ps.index(ps[i])].trace()
                                    for i in range(len(f['trace_lambda_p_square']))]
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
                while (type(f['lambda_p'][p_idx]) == str):
                    p_idx += 1
                while ((index != 1) and (p_idx < len(f['lambda_p']))):
                    print("p_idx = ", p_idx)
                    H = F.order([x for x in f['lambda_p'][:p_idx+1] if type(x) != str])
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
    f.close()
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
            if (ap[:len(inds)] == bp):
                return False
    return True

def is_lift(f, N, prime_bound, default=False):
    primes_N = [p for p in prime_range(prime_bound) if N % p != 0]
    if (len(f['lambda_p_square']) > 1):
        K = f['lambda_p'][0].parent()
        Kx = PolynomialRing(K, name = "x")
        x = Kx.gens()[0]
        lps = [p**6 * x**4 - f['lambda_p'][i]*p**3*x**3 + (f['lambda_p_square'][i] + p**2 + 1)*p*x**2 - f['lambda_p'][i] * x + 1 for i,p in enumerate(primes_N[:len(f['lambda_p_square'])])]
        if any([lp.is_irreducible() for lp in lps]):
            return False, 'G'
        if (len(lps[0].factor()) == 3):
            lift_type = 'P'
        else:
            lift_type = 'Y'
        return True, lift_type
    # Not final, should check these as well
    # There are those where we saved traces, but also multiples of 210
    return default, 'G'

def load_ev_data(k,j,folder,suffix="_", upTo=1000):
    ev_data = [unpickle_omf5(k,j,N,folder,suffix) for N in range(1,upTo)]
    ev_data = [[]] + ev_data
    G_type = [[f for f in ev_data[N] if f['aut_rep_type'] == 'G'] for N in range(upTo)]
    return ev_data, G_type

def update_ev_data(ev_data, B):
    G_type = [[f for f in ev_data[N] if f['aut_rep_type'] == 'G'] for N in range(len(ev_data))]
    for N in range(len(ev_data)):
        print(N)
        for f in G_type[N]:
            is_new = is_newform(f, N, G_type, B)
            if (not is_new):
                f['aut_rep_type'] = 'O'
            is_lft, lift_type = is_lift(f, N, B)
            if (is_lft):
                f['aut_rep_type'] = lift_type
    return

def check_eis(ap_array, N, B):
    eis_tr = [p^3 + p^2 + p + 1 for p in prime_range(B) if N % p != 0]
    return ap_array == eis_tr, []

def check_sk(ap_array, primes_N, divs_N, dim, db, alpha_q={}):
    orbits = db.mf_newforms.search({"level" : {"$in" : divs_N}, "weight" : 4, "char_order" : 1}, ["hecke_orbit_code", "label"])
    orbits = {orb["hecke_orbit_code"] : orb["label"] for orb in orbits}
    for q in alpha_q.keys():
        tr_orbits = db.mf_hecke_traces.search({"hecke_orbit_code" : {"$in" : list(orbits.keys())}, "n" : q, "trace_an" : alpha_q[q]}, "hecke_orbit_code")
        orbits = { tr : orbits[tr] for tr in tr_orbits}
    aps = {hoc : {a["n"] : int(a["trace_an"]) for a in db.mf_hecke_traces.search({"hecke_orbit_code" : hoc})} for hoc in orbits.keys()}
    lamda_ps = { orbits[hoc] : [dim*p*(p+1) + aps[hoc][p] for p in primes_N] for hoc in aps.keys()}
    res = [label for label in lamda_ps.keys() if lamda_ps[label] == ap_array]
    if len(res) == 0:
        return False, []

    assert(len(res) == 1);
    
    return True, [res[0]]

def check_yoshida(ap_array, primes_N, divs_N, dim, db, alpha_q={}):

    wts = [2,4]

    divs_d = divisors(dim)
    
    orbits = { w : db.mf_newforms.search({"level" : {"$in" : divs_N}, "weight" : w, "char_order" : 1, "dim" : {"$in" : divs_d} }, ["hecke_orbit_code", "label", "dim"]) for w in wts}

    orb_dict = {w : { d : {}  for d in divs_d } for w in wts}
    for w in wts:
        for orb in orbits[w]:
            orb_dict[w][orb["dim"]][orb["hecke_orbit_code"]] = orb["label"]

    orbits = orb_dict
    for q in alpha_q.keys():
        tr_orbits = { w : {d : db.mf_hecke_traces.search({"hecke_orbit_code" : {"$in" : list(orbits[w][d].keys())}, "n" : q, "trace_an" : alpha_q[q][w] * d // dim }, "hecke_orbit_code") for d in divs_d} for w in wts}
    
        orbits = { w : {d : { tr : orbits[w][d][tr] for tr in tr_orbits[w][d]} for d in divs_d} for w in wts}

    aps = { w : {d : {hoc : {a["n"] : int(a["trace_an"]) for a in db.mf_hecke_traces.search({"hecke_orbit_code" : hoc})} for hoc in orbits[w][d].keys()} for d in divs_d} for w in wts}

    lamda_ps = {}
    for d in divs_d:
        lamda_ps_d = {(orbits[2][d][hoc2], orbits[4][dim//d][hoc4]) : [p*aps[2][d][hoc2][p]*dim//d + aps[4][dim//d][hoc4][p]*d for p in primes_N] for hoc2 in aps[2][d].keys() for hoc4 in aps[4][dim//d].keys()}
        lamda_ps.update(lamda_ps_d)

    res = [labels for labels in lamda_ps.keys() if lamda_ps[labels] == ap_array]
    
    if len(res) == 0:
        return False, []

    assert(len(res) == 1)
    
    return True, list(res[0])

def aut_rep(f, N, B, db):
    primes_N = [p for p in prime_range(B) if N % p != 0]
    divs_N = divisors(N)
    dim = len(f['field_poly'])-1
    if (len(f['lambda_p']) > 0):
        q = primes_N[0]
        disc = (f['lambda_p'][0])**2-4*q*(f['lambda_p_square'][0]+1-q**2)

        if not is_square(disc):
            return 'G', []

        x1 = (f['lambda_p'][0] + disc.sqrt())/2
        x2 = (f['lambda_p'][0] - disc.sqrt())/2

        y1 = x1/q
        y2 = x2/q
        
        tr_array = [x.trace() for x in f['lambda_p']]
    
        # check eisenstein
        if (((y2 == q + 1) and (x1 == q^3 + 1)) or ((y1 == q + 1) and (x2 == q^3 + 1))):
            is_eis, labels = check_eis(tr_array, N, B)
            if (is_eis):
                return 'F', labels
        
        # check sk
        if (y1 == q + 1):
            alpha_q = {q : ZZ(x2.trace())}
            is_sk, labels = check_sk(tr_array, primes_N, divs_N, dim, db, alpha_q)
            if (is_sk):
                return 'P', labels
        
        if (y2 == q + 1):
            alpha_q = {q : ZZ(x1.trace())}
            is_sk, labels = check_sk(tr_array, primes_N, divs_N, dim, db, alpha_q)
            if (is_sk):
                return 'P', labels

        # check yoshida
        if (y2.is_integral()):
            alpha_q = {q : {2 : ZZ(y2.trace()), 4 : ZZ(x1.trace())} }
            is_yosh, labels = check_yoshida(tr_array, primes_N, divs_N, dim, db, alpha_q)
            if (is_yosh):
                return 'Y', labels
        if (y1.is_integral()):
            alpha_q = {q : {2 : ZZ(y1.trace()), 4 : ZZ(x2.trace())} }
            is_yosh, labels = check_yoshida(tr_array, primes_N, divs_N, dim, db, alpha_q)
            if (is_yosh):
                return 'Y', labels
    else:
        tr_array = f['trace_lambda_p']
        is_eis, labels = check_eis(tr_array, N, B)
        if (is_eis):
                return 'F', labels
        is_sk, labels = check_sk(tr_array, primes_N, divs_N, dim, db)
        if (is_sk):
            return 'P', labels
        is_yosh, labels = check_yoshida(tr_array, primes_N, divs_N, dim, db)
        if (is_yosh):
            return 'Y', labels
                
    return 'G', []

def update_ev_data_class_old(ev_data, B):
    # G_type = [[f for f in ev_data[N] if f['aut_rep_type'] == 'G'] for N in range(len(ev_data))]
    for N in range(len(ev_data)):
        print(N)
        #for f in G_type[N]:
        for f in ev_data[N]:
            is_new = is_newform(f, N, ev_data, B)
            if (not is_new):
                f['aut_rep_type'] = 'O'
            else:
                if all([a.trace() > 0 for a in f['lambda_p']]):
                    f['aut_rep_type'] = 'P'
                elif (N % 210 != 0) and is_yoshida(f,N,B):
                    f['aut_rep_type'] = 'Y'
    return

def update_ev_data_class(ev_data, B, db, start = 0):
    for N in range(start, len(ev_data)):
        print(N)
        for f in ev_data[N]:
            if (f['aut_rep_type'] == 'O'):
                continue
            is_new = is_newform(f, N, ev_data, B)
            if (not is_new):
                f['aut_rep_type'] = 'O'
            elif (N % 210 != 0):
                aut_tp, friends = aut_rep(f,N,B,db)
                f['aut_rep_type'] = aut_tp
                if aut_tp in ['Y', 'P']:
                    f['related_objects'] = friends
    return
