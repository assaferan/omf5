import pickle
from sage.all import (Matrix, PolynomialRing, QQ, NumberField, prime_range, prime_divisors, divisors, is_square, ZZ, sqrt, sign)
from sage.misc.persist import SagePickler

def integer_squarefree_part(n):
    """ returns the squarefree part of the integer n (uses factor rather than calling pari like sage 9.3+ does) """
    return sign(n)*prod([p**(e%2) for p, e in ZZ(n).factor()])

def get_nf_basis(ev, quad_special=False):
    R = PolynomialRing(QQ, 'x')
    K = NumberField(R(ev['field_poly']), 'nu')
    nu = K.gen(0)
    dim = K.degree()

    if ev['field_poly'] == [1,0,1]:
        basis = inv_basis = [1, nu]
        return basis, inv_basis

    # Checking for the two other cyclotomic polynomials of degree 2
    if ev['hecke_ring_power_basis'] and ev['field_poly'] in [[1,-1,1],[1,1,1]]:
        basis = inv_basis = [1, nu]
        return basis, inv_basis
    
    if (dim == 2) and quad_special:
        c, b, a = map(ZZ, ev['field_poly'])
        D = b**2 - 4*a*c
        d = integer_squarefree_part(D)
        s = (D//d).isqrt()
        if ev['hecke_ring_power_basis']:
            k, l = ZZ(0), ZZ(1)
        else:
            k, l = map(ZZ, ev['hecke_ring_numerators'][1])
            k = k / ev['hecke_ring_denominators'][1]
            l = l / ev['hecke_ring_denominators'][1]
        beta = vector((k - (b*l)/(2*a), ((s*l)/(2*a)).abs()))
        den = lcm(beta[0].denom(), beta[1].denom())
        beta *= den
        sqrt = K(d).sqrt()
        num = beta[0] + beta[1]*sqrt
        frac = num / den
        basis = inv_basis = [1, frac]
        return basis, inv_basis
    
    if ev['hecke_ring_power_basis']:
        basis = inv_basis = [nu**i for i in range(dim)]
        return basis, inv_basis
    
    coeff_data = zip(ev['hecke_ring_numerators'], ev['hecke_ring_denominators'])
    basis = [sum([nums[i] * nu**i for i in range(len(nums))])/den for (nums, den) in coeff_data]

    inv_coeff_data = zip(ev['hecke_ring_inverse_numerators'], ev['hecke_ring_inverse_denominators'])
    inv_basis = [sum([nums[i] * nu**i for i in range(len(nums))])/den for (nums, den) in inv_coeff_data]
    
    return basis, inv_basis

def nf_lists_to_elements(coeffs, basis):
    d = len(basis)
    return [sum([coeff[i]*basis[i] for i in range(d)]) if type(coeff) != str else coeff for coeff in coeffs]

def nf_elts_to_lists(elts, inv_basis):
    d = len(inv_basis)
    def to_list(elt):
        if type(elt) == int:
            return [elt] + [0 for i in range(d-1)]
        elif type(elt) == str:
            return elt
        return list(elt)
    return [ list(sum([to_list(elt)[i]*inv_basis[i] for i in range(d)])) if type(elt) != str else elt for elt in elts]

def parse_omf5(k,j,N,folder,suffix="_nl_200_",hecke_ring=True,B=200,max_deg=13,skip_lines=0):
    fname = folder + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    Qx = PolynomialRing(QQ, name="x")
    x = Qx.gens()[0]

    fl = open(fname)
    r = fl.readlines()
    fl.close()
    if (len(r) <= skip_lines):
        return []
    al_signs = eval(r[skip_lines])
    ret = []
    bad_ps = [p for p in prime_range(B) if N % p == 0]
    ps = prime_range(B)
    good_ps = [p for p in ps if p not in bad_ps]
    square_ps = prime_range(sqrt(B))
    good_square_ps = [p for p in square_ps if p not in bad_ps]
    for al_sign in al_signs:
        forms = al_signs[al_sign]
        for f in forms:
            pol = Qx([int(c) for c in f['field_poly'].split()[1:]])
            F = NumberField(pol, name = "a")
            V, from_V, to_V = F.vector_space()
            a = F.gens()[0]
            assert len(good_ps) == len(f['lambda_p'])
            assert len(good_square_ps) == len(f['lambda_p_square'])
            f['lambda_p'] = ['NULL' if ps[i] in bad_ps else F(f['lambda_p'][good_ps.index(ps[i])]) for i in range(len(ps))]
            f['lambda_p_square'] = ['NULL' if ps[i] in bad_ps else F(f['lambda_p_square'][good_ps.index(ps[i])])
                                    for i in range(len(square_ps))]
            # At the moment can't do it on toby - needs to set up LMFDB there
            # aut_tp, friends = aut_rep(f,N,B,db)
            # f['aut_rep_type'] = aut_tp
            # if aut_tp in ['Y', 'P']:
            #   f['related_objects'] = friends
            # !! TODO : represent the eigenvalues in the polredabs field, currently some things break
            coeffs = [int(c) for c in pol.coefficients(sparse=False)]
            f['field_poly'] = coeffs
            f['atkin_lehner_eigenvals'] = [[p, -1 if al_sign % p == 0 else 1] for p in prime_divisors(N)]
            f['atkin_lehner_string'] = ''.join(['-' if al_sign % p == 0 else '+' for p in prime_divisors(N)])
            if (hecke_ring) and (f['aut_rep_type'] != 'O') and (F.degree() <= max_deg): # for old forms currently we may get non-integral stuff
                # !!! This can be very slow
                # hecke_ring = F.order(orbit['lambda_p'])
                index = 0
                p_idx = 0
                nbound = 0
                num_gens = 0
                is_init_H = False
                while (type(f['lambda_p'][p_idx]) == str):
                    p_idx += 1
                while ((index != 1) and (p_idx < len(f['lambda_p'])) and (p_idx < 5)):
                    print("p_idx = ", p_idx)
                    gens = [x for x in f['lambda_p'][:p_idx+1] if type(x) != str]              
                    mod_gens = [to_V(x) for x in gens]
                    ambient = ZZ**V.dimension()
                    W = ambient.span(mod_gens)
                    if (W.rank() == F.degree()):
                        is_init_H = True
                        H = F.order(gens)
                        new_index = H.index_in(F.ring_of_integers())
                        if (new_index != index):
                            index = new_index
                            nbound = p_idx
                    p_idx += 1
                if (is_init_H):
                    f['hecke_ring'] = H
                    f['hecke_ring_index'] = index
                    f['hecke_ring_generator_nbound'] = nbound
                else:
                    f['hecke_ring'] = F.ring_of_integers()
                    f['hecke_ring_index'] = 1
                    f['hecke_ring_generator_nbound'] = 0

                basis = f['hecke_ring'].basis()
                _ = f.pop('hecke_ring')
                f['hecke_ring_power_basis'] = False
                mat = Matrix([list(b) for b in basis])
                f['hecke_ring_denominators'] = [row.denominator() for row in mat]
                f['hecke_ring_numerators'] = [list(row.denominator()*row) for row in mat]  
                f['hecke_ring_inverse_denominators'] = [row.denominator() for row in mat**(-1)]
                f['hecke_ring_inverse_numerators'] = [list(row.denominator()*row) for row in mat**(-1)]   
                inv_coeff_data = zip(f['hecke_ring_inverse_numerators'], f['hecke_ring_inverse_denominators'])
                inv_basis = [sum([nums[i] * a**i for i in range(len(nums))])/den for (nums, den) in inv_coeff_data]
                f['lambda_p'] = nf_elts_to_lists(f['lambda_p'], inv_basis)
                f['lambda_p_square'] = nf_elts_to_lists(f['lambda_p_square'], inv_basis)
            else:
                f['trace_lambda_p'] = ['NULL' if f['lambda_p'][i] == 'NULL' else f['lambda_p'][i].trace() for i in range(len(f['lambda_p']))]
                f['trace_lambda_p_square'] = ['NULL' if f['lambda_p_square'][i] == 'NULL' else f['lambda_p_square'][i].trace() for i in range(len(f['lambda_p_square']))]
                f['lambda_p'] = []
                f['lambda_p_square'] = []
                        
            ret.append(f)
    fl = open(fname, "w")
    fl.write(str(ret))
    fl.close()
    return ret

def pickle_omf5(k,j,N,folder,suffix="_nl_200_",B=200):
    fname = folder + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    parsed = parse_omf5(k,j,N,folder,suffix,B=B)
    f = open(fname, "wb")
    gherkin = SagePickler.dumps(parsed)
    _ = f.write(gherkin)
    f.close()
    # check
    f = open(fname, "rb")
    pickled = f.read()
    assert parsed == pickle.loads(pickled)
    f.close()
    return

def unpickle_omf5(k,j,N,folder,suffix="_"):
    fname = folder + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    f = open(fname, "rb")
    pickled = f.read()
    f.close()
    return pickle.loads(pickled)

def unplain_omf5(k,j,N,folder,suffix="_"):
    fname = folder + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    f = open(fname, "r")
    plain = f.read()
    f.close()
    forms = eval(plain)
    #Qx = PolynomialRing(QQ, name="x")
    #x = Qx.gens()[0]
    for form in forms:
        form['field_poly'] = eval(form['field_poly'])
        #pol = Qx(form['field_poly'])
        #F = NumberField(pol, name = "a")
        #form['lambda_p'] = [F([QQ(y[0])/QQ(y[1]) for y in x]) if x != 'NULL' else x for x in form['lambda_p']]
        #form['lambda_p_square'] = [F([QQ(y[0])/QQ(y[1]) for y in x]) if x != 'NULL' else x for x in form['lambda_p_square']]
        # if 'hecke_ring' in form:
        #    form['hecke_ring_numerators'] = F.order([F([QQ(y[0]) for y in x]) for x in form['hecke_ring']])
        #    form['hecke_ring_denominators'] = F.order([F([QQ(y[1]) for y in x]) for x in form['hecke_ring']])
    return forms

def unplain_omf5_2(k,j,N,folder,suffix="_"):
    fname = folder + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    f = open(fname, "r")
    plain = f.read()
    f.close()
    forms = eval(plain)
    Qx = PolynomialRing(QQ, name="x")
    x = Qx.gens()[0]
    for form in forms:
        pol = Qx(form['field_poly'])
        F = NumberField(pol, name = "a")
        if 'hecke_ring' in form:
            basis = form['hecke_ring'].basis()
            _ = form.pop('hecke_ring')
            form['hecke_ring_power_basis'] = False
            mat = Matrix([list(b) for b in basis])
            form['hecke_ring_denominators'] = [row.denominator() for row in mat]
            form['hecke_ring_numerators'] = [list(row.denominator()*row) for row in mat]  
            form['hecke_ring_inverse_denominators'] = [row.denominator() for row in mat**(-1)]
            form['hecke_ring_inverse_numerators'] = [list(row.denominator()*row) for row in mat**(-1)]   
            # form['hecke_ring'] = [list(x) for x in form['hecke_ring'].basis()]
            # form['hecke_ring'] = [[[y.numerator(), y.denominator()] for y in x] for x in form['hecke_ring']]
            inv_coeff_data = zip(form['hecke_ring_inverse_numerators'], form['hecke_ring_inverse_denominators'])
            inv_basis = [sum([nums[i] * nu**i for i in range(len(nums))])/den for (nums, den) in inv_coeff_data]
            form['lambda_p'] = nf_elts_to_lists(form['lambda_p'], inv_basis)
            form['lambda_p_square'] = nf_elts_to_lists(form['lambda_p_square'], inv_basis)
    return forms

def nf_elts_to_lists(elts, inv_basis):
    d = len(inv_basis)
    def to_list(elt):
        if type(elt) == int:
            return [elt] + [0 for i in range(d-1)]
        elif type(elt) == str:
            return elt
        return list(elt)
    return [ list(sum([to_list(elt)[i]*inv_basis[i] for i in range(d)])) if type(elt) != str else elt for elt in elts]

def py3_pickle_to_plain(k,j,N,folder,suffix="_"):
    fname = folder + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    f = open(fname, "rb")
    pickled = f.read()
    f.close()
    # py2pickle = pickle.dumps(pickle.loads(pickled),2)
    forms = pickle.loads(pickled)
    Qx = PolynomialRing(QQ, name="x")
    x = Qx.gens()[0]
    for form in forms:
        F = NumberField(Qx(form['field_poly']), name = "nu")
        nu = F.gen(0)
        form['hecke_ring_cyclotomic_generator'] = 0
        form['hecke_ring_rank'] = F.degree()
        form['field_poly'] = '[' + ','.join([str(x) for x in form['field_poly']]) + ']'
        if len(form['lambda_p']) > 0:
            form['maxp'] = nth_prime(len(form['lambda_p']))
            form['maxp_square'] = nth_prime(len(form['lambda_p_square']))
            form['trace_lambda_p'] = ['NULL' if form['lambda_p'][i] == 'NULL' else form['lambda_p'][i].trace() for i in range(len(form['lambda_p']))]
            form['trace_lambda_p_square'] = ['NULL' if form['lambda_p_square'][i] == 'NULL' else form['lambda_p_square'][i].trace() for i in range(len(form['lambda_p_square']))]
        else:
            form['maxp'] = nth_prime(len(form['trace_lambda_p']))
            form['maxp_square'] = nth_prime(len(form['trace_lambda_p_square']))
        # form['lambda_p'] = [list(x) if x != 'NULL' else x for x in form['lambda_p']]
        # form['lambda_p'] = [[[y.numerator(), y.denominator()] for y in x] if x != 'NULL' else x for x in form['lambda_p']]
        # form['lambda_p_square'] = [list(x) if x != 'NULL' else x for x in form['lambda_p_square']]
        # form['lambda_p_square'] = [[[y.numerator(), y.denominator()] for y in x] if x != 'NULL' else x for x in form['lambda_p_square']]
        if 'hecke_ring' in form:
            basis = form['hecke_ring'].basis()
            _ = form.pop('hecke_ring')
            form['hecke_ring_power_basis'] = False
            mat = Matrix([list(b) for b in basis])
            form['hecke_ring_denominators'] = [row.denominator() for row in mat]
            form['hecke_ring_numerators'] = [list(row.denominator()*row) for row in mat]  
            form['hecke_ring_inverse_denominators'] = [row.denominator() for row in mat**(-1)]
            form['hecke_ring_inverse_numerators'] = [list(row.denominator()*row) for row in mat**(-1)]   
            # form['hecke_ring'] = [list(x) for x in form['hecke_ring'].basis()]
            # form['hecke_ring'] = [[[y.numerator(), y.denominator()] for y in x] for x in form['hecke_ring']]
            inv_coeff_data = zip(form['hecke_ring_inverse_numerators'], form['hecke_ring_inverse_denominators'])
            inv_basis = [sum([nums[i] * nu**i for i in range(len(nums))])/den for (nums, den) in inv_coeff_data]
            form['lambda_p'] = nf_elts_to_lists(form['lambda_p'], inv_basis)
            form['lambda_p_square'] = nf_elts_to_lists(form['lambda_p_square'], inv_basis)
        else:
            # For now, if we don't have a basis for the hecke ring, we don't present the lambdas
            form['lambda_p'] = []
            form['lambda_p_square'] = []
    f = open(fname, "w")
    _ = f.write(str(forms))
    f.close()
    return

def unify_loc_and_repo(k,j,N,folder_loc,folder_repo,suffix="_",B=100, max_deg=20):
    forms_loc = unpickle_omf5(k,j,N,folder_loc)
    forms_repo = unplain_omf5(k,j,N,folder_repo)
    bad_ps = [p for p in prime_range(B) if N % p == 0]
    ps = prime_range(B)
    good_ps = [p for p in ps if p not in bad_ps]
    square_ps = prime_range(B.sqrt())
    good_square_ps = [p for p in square_ps if p not in bad_ps]
    for idx in range(len(forms_loc)):
        f = forms_loc[idx]
        g = forms_repo[idx]
        # verifying this is the same form
        assert f['atkin_lehner_string'] == g['atkin_lehner_string']
        assert f['field_poly'] == g['field_poly']
        if len(f['field_poly']) <= max_deg:
            assert [f['lambda_p'][ps.index(good_ps[i])] for i in range(len(good_ps))] == g['lambda_p']
            assert [f['lambda_p_square'][ps.index(good_square_ps[i])] for i in range(len(good_square_ps))] == g['lambda_p_square']
        else:
            assert [f['lambda_p'][ps.index(good_ps[i])].trace() for i in range(len(good_ps))] == g['trace_lambda_p']
            assert [f['lambda_p_square'][ps.index(good_square_ps[i])].trace() for i in range(len(good_square_ps))] == g['trace_lambda_p_square']
        f['aut_rep_type'] = g['aut_rep_type']
        if 'related_objects' in g:
            f['related_objects'] = g['related_objects']
    gherkin = SagePickler.dumps(forms_loc)
    f = open(fname, "wb")
    _ = f.write(gherkin)
    f.close()
    return

def is_newform(f, N, G_type, prime_bound):
    divs = divisors(N)[:-1]
    primes_N = [p for p in prime_range(prime_bound) if N % p != 0]
#    if (len(f['lambda_p']) > 0):
#        ap = [x.trace() for x in f['lambda_p']]
#    else:
#        ap = f['trace_lambda_p']
    ap = [x for x in f['trace_lambda_p'] if x != 'NULL']
    for d in divs:
        primes_d = [p for p in prime_range(prime_bound) if d % p != 0]
        inds = [i for i in range(len(primes_d)) if primes_d[i] in primes_N]
        for g in G_type[d]:
#            if (len(g['lambda_p']) > 0):
#                bp = [g['lambda_p'][i] for i in inds]
#                bp = [x.trace() for x in bp]
#            else:
#                bp = [g['trace_lambda_p'][i] for i in inds]
            g_trace = [x for x in g['trace_lambda_p'] if x != 'NULL']
            #bp = [g['trace_lambda_p'][i] for i in inds]
            bp = [g_trace[i] for i in inds]
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
    # ev_data = [unpickle_omf5(k,j,N,folder,suffix) for N in range(1,upTo)]
    ev_data = [unplain_omf5(k,j,N,folder,suffix) for N in range(1,upTo)]
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
#    if (len(f['lambda_p']) > 0):
#        q = primes_N[0]
#        disc = (f['lambda_p'][0])**2-4*q*(f['lambda_p_square'][0]+1-q**2)
#        if not is_square(disc):
#            return 'G', []
#        x1 = (f['lambda_p'][0] + disc.sqrt())/2
#        x2 = (f['lambda_p'][0] - disc.sqrt())/2
#        y1 = x1/q
#        y2 = x2/q
#        tr_array = [x.trace() for x in f['lambda_p']]
# check eisenstein
#        if (((y2 == q + 1) and (x1 == q^3 + 1)) or ((y1 == q + 1) and (x2 == q^3 + 1))):
#            is_eis, labels = check_eis(tr_array, N, B)
#            if (is_eis):
#                return 'F', labels
# check sk
#        if (y1 == q + 1):
#            alpha_q = {q : ZZ(x2.trace())}
#            is_sk, labels = check_sk(tr_array, primes_N, divs_N, dim, db, alpha_q)
#            if (is_sk):
#                return 'P', labels  
#        if (y2 == q + 1):
#            alpha_q = {q : ZZ(x1.trace())}
#            is_sk, labels = check_sk(tr_array, primes_N, divs_N, dim, db, alpha_q)
#            if (is_sk):
#                return 'P', labels
# check yoshida
#        if (y2.is_integral()):
#            alpha_q = {q : {2 : ZZ(y2.trace()), 4 : ZZ(x1.trace())} }
#            is_yosh, labels = check_yoshida(tr_array, primes_N, divs_N, dim, db, alpha_q)
#            if (is_yosh):
#                return 'Y', labels
#        if (y1.is_integral()):
#            alpha_q = {q : {2 : ZZ(y1.trace()), 4 : ZZ(x2.trace())} }
#            is_yosh, labels = check_yoshida(tr_array, primes_N, divs_N, dim, db, alpha_q)
#            if (is_yosh):
#                return 'Y', labels
#    else:
    tr_array = [ x for x in f['trace_lambda_p'] if x != 'NULL']
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
            # elif (N % 210 != 0):
            else:
                aut_tp, friends = aut_rep(f,N,B,db)
                f['aut_rep_type'] = aut_tp
                if aut_tp in ['Y', 'P']:
                    f['related_objects'] = friends
    return

def unite(k,j,N, folder_orig, folder_new, suffix = "_"):
    orig_fname = folder_orig + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    new_fname = folder_new + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    orig_file = open(orig_fname)
    new_file = open(new_fname)
    orig_forms = eval(orig_file.read())
    new_forms = eval(new_file.read())
    orig_file.close()
    new_file.close()
    assert len(orig_forms) == len(new_forms)
    for i in range(len(orig_forms)):
        orig_f = orig_forms[i]
        new_f = new_forms[i]
        assert orig_f['field_poly'] == new_f['field_poly']
        new_f['aut_rep_type'] = orig_f['aut_rep_type']
        if ('lambda_p' in new_f) and (len(new_f['lambda_p']) > 0):
            basis, inv_basis = get_nf_basis(new_f)
            lambdas = nf_lists_to_elements(new_f['lambda_p'], basis)
            lambdas2 = nf_lists_to_elements(new_f['lambda_p_square'], basis)
            if 'trace_lambda_p' not in new_f:
                new_f['trace_lambda_p'] =  ['NULL' if lam == 'NULL' else lam.trace() for lam in lambdas]
                new_f['trace_lambda_p_square'] =  ['NULL' if lam == 'NULL' else lam.trace() for lam in lambdas2]
    new_file = open(new_fname, "w")
    new_file.write(str(new_forms))
    new_file.close()
    return

def fix_quadratic_bug_and_make_friends(N, db):
    folder = "../omf5_data/hecke_evs_3_0/data/"
    space = unplain_omf5_2(3,0,N,folder)
    k = 3
    j = 0
    suffix = "_"
    fname = folder + "hecke_ev_%d_%d%s%d.dat" %(k,j,suffix,N)
    space_file = open(fname, "r")
    space_data = eval(space_file.read())
    space_file.close()
    for i,f_data in enumerate(space_data):
        f = space[i]
        if (len(f['field_poly']) == 3):
            basis, inv_basis = get_nf_basis(f)
            lambdas = nf_lists_to_elements(f['lambda_p'], basis)
            lambdas2 = nf_lists_to_elements(f['lambda_p_square'], basis)
            f_data['trace_lambda_p'] = ['NULL' if lambdas[i] == 'NULL' else lambdas[i].trace() for i in range(len(f['lambda_p']))]
            f_data['trace_lambda_p_square'] = ['NULL' if lambdas2[i] == 'NULL' else lambdas2[i].trace() for i in range(len(f['lambda_p_square']))]
        if f['aut_rep_type'] in ['Y','P']:
            aut_tp, friends = aut_rep(f_data, N, 100, db)
            assert aut_tp == f['aut_rep_type']
            f_data['related_objects'] = friends
        else:
            f_data['related_objects'] = []
    space_file = open(fname, "w")
    _ = space_file.write(str(space_data))
    space_file.close()
