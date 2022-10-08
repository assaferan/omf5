# -*- mode: sage-shell:sage -*-
# Implementation of Dubey Holenstein method
# Note that they use the modern integral definition, not classically integral
# This makes this largely incompatible with Sage's Quadratic form definition?
# However Sage's Genus class is compatible.

# Did I mention how much I hate this subject?

#The normal form functions are copied from Brandon's code. Oh wait, those use Sage's definitions.. Adaptation time
import numpy
import sage.arith.misc

# This should not be necessary
def sqrt_mod_q(a, q):
    p = prime_factors(q)[0]
    if p==2:
        return sqrt_mod_even(a,q)
    else:
        return sqrt_mod_primepow(a,p, q)

def sqrt_mod_primepow(a, p, q):
    olda = a
    a = a.lift()
    if a==0:
        return 0
    k = valuation(q, p)
    val = valuation(a,p)
    if val %2 == 1:
        assert False
    retval = val/2
    t = a/(p**val)
    k = k-val
    s = Mod(t, p).sqrt().lift()
    for i in range(2, k+1):
        corr = (s*s-t) %(p**i)
        s -= (corr/(2*s)) % (p**i)
    ret = Mod(s*p**(retval), q)
    assert ret*ret == olda
    return ret

def sqrt_mod_even(a, q):
    olda = a
    a = a.lift()
    if a==0:
        return 0
    k = valuation(q, 2)
    val = valuation(a, 2)
    if val %2 == 1:
        assert False
    retval = val/2
    t = a/(2**val)
    if q<=8:
        for i in range(0, 7):
            if i*i % q ==a:
                return Mod(i, q)
    s = Mod(t, 8).sqrt().lift()
    for i in range(4, k+1):
        corr = (s*s-t) %(2**i)
        s -= (corr/(2*s)) % (2**i)
    ret = Mod(s*2**(retval), q)
    assert ret*ret == olda
    return ret

def twonf_with_change_vars_correct(Q):
    #basically copied from the built-in function local_normal_form(), specialized to p=2, but it also computes a change-of-variables to the local normal form. it skips the "Cassels' proof" step. this is accounted for later.
    I = list(range(Q.dim()))
    P = identity_matrix(QQ,Q.dim())
    Q_Jordan = DiagonalQuadraticForm(ZZ,[])
    while Q.dim() > 0:
        n = Q.dim()
        (min_i,min_j) = Q.find_entry_with_minimal_scale_at_prime(2)
        if min_i == min_j:
            min_val = valuation(2 * Q[min_i,min_j],2)
        else:
            min_val = valuation(Q[min_i, min_j],2)
        if (min_i == min_j):
            block_size = 1
            Q.swap_variables(0,min_i,in_place = True)
            P.swap_rows(I[0],I[min_i])
        else:
            Q.swap_variables(0, min_i, in_place = True)
            Q.swap_variables(1, min_j, in_place = True)
            P.swap_rows(I[0],I[min_i])
            P.swap_rows(I[1],I[min_j])
            block_size = 2
        min_scale = 2^(min_val)
        if (block_size == 1):
            a = 2 * Q[0,0]
            for j in range(block_size, n):
                b = Q[0, j]
                g = GCD(a, b)
                Q.multiply_variable(ZZ(a/g), j, in_place = True)
                Q.add_symmetric(ZZ(-b/g), j, 0, in_place = True)
                P.rescale_row(I[j],a/g)
                P.add_multiple_of_row(I[j],I[0],-b/g)
            NewQ = QuadraticForm(matrix([a]))
        else:
            a1 = 2*Q[0,0]
            a2 = Q[0,1]
            b2 = 2*Q[1,1]
            big_det = (a1*b2 - a2^2)
            small_det = big_det / (min_scale * min_scale)
            for j in range(block_size,n):
                a = Q[0,j]
                b = Q[1,j]
                Q.multiply_variable(big_det,j,in_place = True)
                Q.add_symmetric(ZZ(-(a*b2 - b*a2)),j,0,in_place = True)
                Q.add_symmetric(ZZ(-(-a*a2 + b*a1)),j,1,in_place = True)
                Q.divide_variable(ZZ(min_scale^2),j,in_place = True)
                P.rescale_row(I[j],big_det)
                P.add_multiple_of_row(I[j],I[0],-a*b2+b*a2)
                P.add_multiple_of_row(I[j],I[1],a*a2-b*a1)
                P.rescale_row(I[j],min_scale^(-2))
            Qgcd = gcd([Q[0,0],Q[0,1],Q[1,1]])
            Qgcdsqf = Qgcd.squarefree_part()
            fsq = ZZ(Qgcd / Qgcdsqf)
            fsq = fsq / 2^(fsq.valuation(2))
            f = ZZ(sqrt(Qgcd / Qgcdsqf))
            f = f / 2^(f.valuation(2))
            NewQ = QuadraticForm(matrix([[ZZ(2*Q[0,0]/fsq),ZZ(Q[0,1]/fsq)],[ZZ(Q[0,1]/fsq),ZZ(2*Q[1,1]/fsq)]]))
            P.rescale_row(I[0],f^(-1))
            P.rescale_row(I[0],f^(-1))
        Q_Jordan = Q_Jordan + NewQ
        for j in range(block_size):
            I.remove(I[0])
        Q = Q.extract_variables(range(block_size, n))
    return [Q_Jordan,P]

def local_normal_form_with_change_vars(Q,p):
    #basically copied from the built-in function local_normal_form(), but it also computes a change-of-variables to the local normal form
    I = list(range(Q.dim()))
    assert p != 2
    M = identity_matrix(QQ,Q.dim())
    Q_Jordan = DiagonalQuadraticForm(ZZ,[])
    while Q.dim() > 0:
        n = Q.dim()
        (min_i, min_j) = Q.find_entry_with_minimal_scale_at_prime(p)
        if min_i == min_j:
            min_val = valuation(2 * Q[min_i, min_j], p)
        else:
            min_val = valuation(Q[min_i, min_j], p)
        if (min_i == min_j):
            Q.swap_variables(0, min_i, in_place = True)
            M.swap_rows(I[0],I[min_i])
        else:
            Q.swap_variables(0, min_i, in_place = True)
            Q.swap_variables(1, min_j, in_place = True)
            M.swap_rows(I[0],I[min_i])
            M.swap_rows(I[1],I[min_j])
            Q.add_symmetric(1, 0, 1, in_place = True)
            M.add_multiple_of_row(I[0],I[1],1)
        min_scale = p ** min_val
        a = 2 * Q[0,0]
        for j in range(1, n):
            b = Q[0, j]
            g = GCD(a, b)
            Q.multiply_variable(ZZ(a/g), j, in_place = True)
            Q.add_symmetric(ZZ(-b/g), j, 0, in_place = True)
            M.rescale_row(I[j],a/g)
            M.add_multiple_of_row(I[j],I[0],-b/g)
        Q_Jordan = Q_Jordan + Q.extract_variables(range(1))
        I.remove(I[0])
        Q = Q.extract_variables(range(1, n))
    return [Q_Jordan,M]

def make_similar_quadratic_form(A):
    M=Matrix(ZZ, A.dimensions()[0], A.dimensions()[1], 0)
    for i in range(0, A.dimensions()[0]):
        for j in range(0, i):
            M[i,j]=A[i,j]+A[j,i]
            M[j,i]=M[i,j]
    for i in range(0, A.dimensions()[0]):
        M[i,i] = 2*A[i,i]
    return QuadraticForm(ZZ, M)

def locally_block_split(A, q):
    r"""
    Given a symmetric matrix $A$ with entries in $\mathbb{Z}$, return $U, D$ st $U^{T}AU=D \pmod{q}$ and $D$ is in blocks.

    q must be a prime power.
    """
    if q%2==0:
        return locally_block_split_2adic(A, q)
    else:
        return locally_block_split_padic(A, q)

def locally_block_split_padic(A, q):
    Q = make_similar_quadratic_form(A);
    assert len(prime_factors(q))==1
    p = prime_factors(q)[0]
    [_,U]=local_normal_form_with_change_vars(Q, p)
    return [U.transpose(), U*A*U.transpose()]

def locally_block_split_2adic(A, q):
    Q = make_similar_quadratic_form(A);
    assert len(prime_factors(q))==1
    p = prime_factors(q)[0]
    assert p == 2
    [_,U]=twonf_with_change_vars_correct(Q)
    return [U.transpose(), U*A*U.transpose()]


# Block splitting is all we need. Next we identify pieces, find a solution, apply U and then CRT to solve completely. I don't quite get how this works exactly.
# We just need one, so try every single splitting.
# Understand paper more: we want a primitive splitting, and I'd like to know if one exists or not

# The key is to look at the p-adic symbols. We store as pairs of (ord, symb) where ord is an integer representing the power of p^k
# symb is 1,3,5,7 for p=2, and 1 or -1 if p odd

def power_mod_form(q):
    p = prime_factors(q)[0]
    return (p, valuation(q, p))

def symbol_of(t, q):
    t = Mod(t, q).lift()
    p = prime_factors(q)[0]
    if p==2:
        ordert=valuation(t, p)
        if t == 0:
            symbolt = 0
        else:
            symbolt = (t//(p**ordert)) % 8
        return (ordert, symbolt)
    else:
        ordert=valuation(t, p)
        if t == 0:
            symbolt = 0
        else:
            symbolt = kronecker_symbol(t//(p**ordert), p)
        return ( ordert, symbolt)

def splitting_works(s1, s2, s3, q):
    ''' Returns True if there is are t, a, b in s1, s2, s3 st t=a+b, false if not '''
    # Woops, what about 0? 
    # NB that s1 is \gamma, s2 is \gamma_1, s3 is \gamma_2 in paper notation
    (val1, sym1)=s1
    (val2, sym2)=s2
    (val3, sym3)=s3
    p = prime_factors(q)[0]
    # Several special cases with infinity need to be handled
    if val1 == Infinity:
        if val2 != val3:
            return False
        else:
            if val2 == Infinity:
                return val3==Infinity
            else:
                ## Mistake in paper: this requires paying attention to the symbol and paper doesn't.
                if p==2:
                    return s2 == symbol_of((-(2**val3)*sym3)%q,q)
                else:
                    return sym2 == kronecker_symbol(-1, p)*sym3
    if val2 == Infinity:
        if s1 == s3:
            return True
        else:
            return False
    if val3 == Infinity:
        return s1 == s2

    # Now we know that infinity isn't present
    if val1==val2 and val1==val3:
        if p == 2:
            return False
        else:
            sols=p - (p %4) - (( sym2+sym3)*(kronecker_symbol(-1, p)*sym2+sym1))
            return sols != 0
    if val1 != val3:
        #swap to simplify next case
        (val3, val2) = (val2, val3)
        (sym3, sym2) = (sym2, sym3)
    # Casework from Lemma 21
    gamma_3 = 0
    if p != 2:
        if val1<val2:
            gamma_3 = s1
        else:
            gamma_3 = (val2, kronecker_symbol(-1, p)*sym2)   
    else:
        if val1<val2:
            magic_s_1 = (sym1 - (2 ** (val2-val1))*sym2)% 8
            gamma_3 = (val1, magic_s_1)
        else:
            magic_s_2 = (sym1*(2**(val1-val2)) - sym2)% 8
            gamma_3 = (val2, magic_s_2)
    return gamma_3 == (val3, sym3)

def sample_from_squareclass(s1, q):
    p = prime_factors(q)[0]
    ordq = valuation(q, p)
    if s1[0] == Infinity:
        return 0
    if p == 2:
        return (p**s1[0]*s1[1])%q
    else:
        t = 0
        while (kronecker_symbol(t, p)!=s1[1]):
            t = randrange(1, p)
        return (p**s1[0]*t)%q
    
def find_split(t, s1, s2, q):
    ''' Returns (a,b) with a inhabiting symbol s1, b inhabiting symbol s2 and t=a+b mod q'''
    #Insert test assertions
    if not splitting_works(symbol_of(t,q), s1, s2, q):
        assert False
    # Now what to do? Paper is somewhat unclear, but when valuations different any work
    p = prime_factors(q)[0]
    valt = valuation(t, p)
    if valt != s1[0]:
        a = sample_from_squareclass(s1, q)
        return (a, (t-a)%q)
    elif valt != s2[0]:
        a = sample_from_squareclass(s2, q)
        return ((t-a)%q, a)
    # Valuations same we apply lemma 26
    if valt == Infinity:
        return (0, 0)
    valc = s1[0]
    valq = valuation(q, p)
    crpt = t/p**valt
    done = False
    while not done:
        tau_1 = 0
        while gcd(tau_1, p) != 1:
            tau_1 = randrange(1, p**(valq-valc))
        tau_2 = (crpt - tau_1) %q
        if symbol_of(tau_1*p**valc, q) == s1 and symbol_of(tau_2*p**valc, q) == s2:
            done = True
    return (tau_1*p**valc, tau_2*p**valc)

# How to think about blocks? Starting point+type? Or by coeffs?
# Here by coeffs. s is a type, t is an actual number
# Possible refactoring

# For type I the condition is simple. Solving also easy. So why break up types?
def solveable_type_I(s, block, q):
    p = prime_factors(q)[0]
    ordb = valuation(block, p)
    ords = s[0]
    if ordb == Infinity:
        if ords == Infinity:
            return True
        else:
            return False
    if ords == Infinity:
        return True #Check
    if ords - ordb >= 0 and (ordb-ords) %2 ==0:
        if p != 2:
            return s[1] ==symbol_of(block, q)[1]
        else:
            nmod = min(8, 2**(valuation(q, 2)-s[0]))
            return s[1]% nmod == symbol_of(block, q)[1]% nmod
    else:
        return False

def solveable_type_I_prim(s, block, q):
    return s == symbol_of(block, q)

def solution_type_I(t, block, q):
    if not solveable_type_I(symbol_of(t,q), block, q):
        assert False
    p = prime_factors(q)[0]
    ordb = valuation(block, p)
    ordt = valuation(t, p)
    if ordt == Infinity:
        return 0
    tprime = t/(p**(ordt))
    blockprime = block/(p**ordb)
    if p == 2:
        k = valuation(q, 2)
        littlemod = 2**(k-ordt)
        y=(Mod(tprime, littlemod)/Mod(blockprime, littlemod))
        y = sqrt_mod_q(y, littlemod).lift()
        sqrp = y
    else:
        sqrp=(Mod(tprime,q)/Mod(blockprime,q))
        sqrp = sqrt_mod_q(sqrp, q).lift()
    ans = sqrp * p**((ordt-ordb)/2)
    if Mod(ans, q)**2*Mod(block, q) != Mod(t, q):
        assert False
    return ans

# The paper has a complicated counting algorithm that determines both of these together recusively
# important: we write our b as the coefficient of the quadratic form, not matrix element i.e. 2x matrix element

# Need to make lots of tests
def count_prim_nonprim_type_II(s, a, b, c, q):
    lp1 = valuation(b, 2) #This line is confusing to me also. Results from choice above
    k = valuation(q, 2)
    if lp1 >= k:
        if s[0] >= k:
            nonprim = 4**(k-1)
            prim = 4**k - nonprim
            return (prim, nonprim)
        else:
            return (0, 0)
    elif s[0]<lp1:
        return (0,0)
    else:
        prim = 0
        nonprim = 0
        k = k -lp1
        s = (s[0]-lp1, s[1])
        a = a/(2**lp1)
        b = b/(2**lp1)
        c = c/(2**lp1)
        # Now need to reconstruct s mod 2 for this next step
        if s[0] == Infinity:
            sm = 0
        else:
            sm = (s[1]*2**s[0])
        for pair in [(0,1), (1, 0), (1, 1)]:
            (a1, a2) = pair
            if Mod(a*a1**2+b*a1*a2+c*a2**2, 2) == Mod(sm, 2):
                prim += 2**(k-1)
        if sm % 4 == 0:
            if k == 1:
                nonprim += 1
            else:
                (n1,n2) = count_prim_nonprim_type_II((s[0]-2, s[1]), a,b,c, 2**(k-2))
                nonprim += nonprim +4*(n1+n2)
        prim = 4**(lp1)*prim
        nonprim = 4**(lp1)*nonprim
        return (prim, nonprim)

def solveable_type_II(s, a, b, c, q):
    (n1, n2)=count_prim_nonprim_type_II(s, a, b, c, q)
    return n1+n2>0

def solveable_type_II_prim(s, a, b, c, q):
    (n1, n2)=count_prim_nonprim_type_II(s, a, b, c, q)
    return n1>0

def solution_type_II(t, a, b, c, q):
    solution=(0,0)
    if q == 1:
        return (0, 0)
    sym = symbol_of(t, q)
    assert solveable_type_II(symbol_of(t, q), a, b, c, q)
    if solveable_type_II_prim(symbol_of(t, q), a,b, c,q):
        return solution_type_II_prim(t, a, b, c, q)
    else:
        # It isn't primitive solveable because we need to scale the solution
        assert Mod(q, 4) == 0
        assert Mod(t, 4) == 0
        (a1, a2)= solution_type_II(t/4, a, b, c, q/4)
        solution = (a1*2, a2*2)
    z = solution[0]
    y = solution[1]
    assert Mod(a*z*z+b*z*y+c*y*y, q)==Mod(t, q)
    return solution

def solution_type_II_prim(t, a, b, c, q):
    assert solveable_type_II_prim(symbol_of(t, q), a, b, c,q)
    solution=(0,0)
    lp1 = valuation(b, 2)
    if lp1 >= valuation(q,2):
        if valuation(t, 2) >=lp1:
            solution = (0,1) #Will solve equation, is primative
        else:
            assert False #Unsolveable case
    else:
        #rescale
        solution=(0,0)
        tp = t/2**lp1
        ap = a/2**lp1
        bp = b/2**lp1
        cp = c/2**lp1
        exp = valuation(q,2)-lp1
        #Find starting solution
        for tmp in [(0,1), (1, 1), (1, 0)]:
            (a1, a2) = tmp
            if Mod(ap*a1**2+bp*a1*a2+cp*a2**2, 2) == Mod(tp,2):
                solution = tmp
        if solution==(0,0):
            assert False
        #Henselize
        for i in range(1, exp+1):
            (y1, y2) = solution
            b1 = 0
            b2 = 0
            if Mod(y1, 2) == 1:
                b1 = 0
                b2 = Mod(((tp - (ap*y1*y1+bp*y1*y2+cp*y2*y2))/(bp*2**i)),2).lift()
            else:
                b2 = 0
                b1 = Mod(((tp - (ap*y1*y1+bp*y1*y2+cp*y2*y2))/(bp*2**i)),2).lift()
            solution=(y1+b1*(2**i), y2+b2*(2**i))
            z1 = solution[0]
            z2 = solution[1]
    z = solution[0]
    y = solution[1]
    assert Mod(a*z*z+b*z*y+c*y*y, q) == Mod(t, q)
    return solution

def random_test_type_II(q,lbound, coeffsize):
    # Create random type II block
    # Random x, y
    # Then ask for solution
    l=randrange(0, lbound)
    a=randrange(0, coeffsize)
    b=randrange(0, coeffsize)
    c=randrange(0, coeffsize)
    b=2*b+1
    a=a*(2**l)
    b=b*(2**l)
    c=c*(2**l)
    x=randrange(0, coeffsize*2**l)
    y=randrange(0, coeffsize*2**l)
    t=a*x*x+b*x*y+c*y*y
    t=Mod(t, q).lift()
    assert solveable_type_II(symbol_of(t, q), a,b,c, q)
    solution=solution_type_II(t, a, b, c, q)
    z = solution[0]
    w = solution[1]
    assert Mod(a*z*z+b*w*z+c*w*w,q)==Mod(t, q)

def random_test_find_split(q, coeffsize):
    a = randrange(0, coeffsize)
    b = randrange(0, coeffsize)
    t = (a+b)% q
    s = symbol_of(t, q)
    s1 = symbol_of(a, q)
    s2 = symbol_of(b, q)
    assert splitting_works(s, s1, s2, q)
    (c, d)=find_split(t, s1, s2, q)
    assert symbol_of(c,q)==s1
    assert symbol_of(d,q)==s2
    assert Mod(c+d, q) == t

def random_test_block_I(q, coeffsize):
    block = randrange(0, coeffsize)
    x = randrange(0, coeffsize)
    val = Mod(block*x*x, q).lift()
    y = solution_type_I(val, block, q)
    assert Mod(block*y*y, q)==Mod(val, q)
    
# The algorithm is now clear: search for a splitting that works, then go through and carry out with actual t, accumulate results
# We take each block, do it mod q, lift, and then represent the form as a list of blocks
# Then we apply dynamic programing to find a possible set of splits
# Then we find splits and representations, and put together to get a vector
# Lastly we apply U in reverse and take mod q
# And then do a final check

def find_modq_representation(A, t,q):
    tmp=locally_block_split(A, q)
    U=tmp[0]
    D=tmp[1]
    summands = form_into_subblocks(D,q)
    typelist=find_modq_rep_types_recurse(summands, symbol_of(t,q), q)
    if typelist == False:
        assert False
    vec = find_modq_rep_sum(summands,typelist,t,q)
    assert Mod(vec.inner_product(D*vec),q) ==Mod(t, q)
    ans=U*vec
    assert Mod(ans.inner_product(A*ans), q)==Mod(t, q)
    return ans

def is_modq_representable(A, t, q):
    return is_modq_rep_type(A, symbol_of(t, q), q)

def is_modq_rep_type(A, sym, q):
    tmp=locally_block_split(A, q)
    U=tmp[0]
    D=tmp[1]
    summands = form_into_subblocks(D,q)
    typelist=find_modq_rep_types_recurse(summands, sym, q)
    if typelist == False:
        return False
    else:
        return True

def form_into_subblocks(A, q):
    n = A.dimensions()[0]
    blocklist = list() # Blocklist has list either of tuples of 1 number or 3
    i = 0
    while i<n:
        a=Mod(A[i,i],q).lift()
        if i<n-1:
            b=Mod(2*A[i,i+1],q).lift()
            if b != 0:
                c=Mod(A[i+1, i+1],q).lift()
                blocklist.append([a, b,c])
                i+=2
            else:
                blocklist.append([a])
                i+=1
        else:
            blocklist.append([a])
            i+=1
    return blocklist

def solveable_block_any(block, symb, q):
    if len(block)==1:
        return solveable_type_I(symb, block[0], q)
    else:
        return solveable_type_II(symb, block[0], block[1], block[2], q)

def solveable_block_prim(block, symb, q):
    if len(block)==1:
        return solveable_type_I_prim(symb, block[0], q)
    else:
        return solveable_type_II_prim(symb, block[0], block[1], block[2], q)
#Todo: make this greedy and solve with as few blocks as possible
# Figure out speed?

def solution_symbol_prim(Q, s, q):
    p=prime_factors(q)[0]
    exp = valuation(q, p)
    allowed_symbols=list()
    allowed_symbols.append((Infinity, 0))
    for i in range(0, exp):
        if p == 2:
            allowed_symbols.append((i, 1))
            allowed_symbols.append((i, 3))
            allowed_symbols.append((i, 5))
            allowed_symbols.append((i, 7))
        else:
            allowed_symbols.append((i, 1))
            allowed_symbols.append((i, -1))
    if len(Q)==1:
        if solveable_block_prim(Q[0], s, q):
            return [(s, "end")]
        else:
            return False
    else:
        for ext in allowed_symbols:
            if solveable_block_prim(Q[0], ext, q):
                for res in allowed_symbols:
                    if splitting_works(s, ext, res, q):
                        recurse = solution_symbol_any(Q[1:], res, q)
                        if recurse != False:
                            recurse.insert(0, (ext, res))
                            return recurse

        for ext in allowed_symbols:
            if solveable_block_any(Q[0], ext, q):
                for res in allowed_symbols:
                    if splitting_works(s, ext, res, q):
                        recurse = solution_symbol_prim(Q[1:], res, q)
                        if recurse != False:
                            recurse.insert(0, (ext, res))
                            return recurse
    return False

def solution_symbol_any(Q, s, q):
    p=prime_factors(q)[0]
    exp = valuation(q, p)
    allowed_symbols=list()
    allowed_symbols.append((Infinity, 0))
    for i in range(0, exp):
        if p == 2:
            allowed_symbols.append((i, 1))
            allowed_symbols.append((i, 3))
            allowed_symbols.append((i, 5))
            allowed_symbols.append((i, 7))
        else:
            allowed_symbols.append((i, 1))
            allowed_symbols.append((i, -1))
    if len(Q)==1:
        if solveable_block_any(Q[0], s, q):
            return [(s, "end")]
        else:
            return False
    else:
        for ext in allowed_symbols:
            if solveable_block_any(Q[0], ext, q):
                for res in allowed_symbols:
                    if splitting_works(s, ext, res, q):
                        recurse = solution_symbol_any(Q[1:], res, q)
                        if recurse != False:
                            recurse.insert(0, (ext, res))
                            return recurse
        return False

def find_modq_rep_types_recurse(Q, symb, q):
    prim = solution_symbol_prim(Q, symb, q)
    if prim == False:
        return prim
    prim.insert(0, ("empty symbol", symb))
    return prim

#Test these, adapt rest of program
def find_modq_rep_types(Q, t,q):
    p=prime_factors(q)[0]
    exp = valuation(q, p)
    allowed_symbols=list()
    allowed_symbols.append((Infinity, 0))
    for i in range(0, exp):
        if p == 2:
            allowed_symbols.append((i, 1))
            allowed_symbols.append((i, 3))
            allowed_symbols.append((i, 5))
            allowed_symbols.append((i, 7))
        else:
            allowed_symbols.append((i, 1))
            allowed_symbols.append((i, -1))
    summand_count = len(Q)
    primitive_solutions = list()
    imprimitive_solutions = [[("empty symbol",symbol_of(t, q))]]
    # Redo this as DFS
    # Just use a helper function for it that extends solutions either making them primative
    # or not
    # two helper functions
    # Simple recusion, just write it: recursively solving the problem for a split or not
    for i in range(0, summand_count):
        # Update the current solution lists by seeing if they can extend
        # (Think through definition and edge and invariants)
        # Each solution is a list of pairs. The first pair is special
        # The i+1 pair (ext, res) is something so that res[i] splits into ext and res
        # and ext satisfies the Q[i] block
        # final res is special
        new_prim_sol = list()
        new_imprim_sol = list()
        for solution in primitive_solutions:
            if i<summand_count-1:
                tail=solution[i][1]
                for ext in allowed_symbols:
                    for res in allowed_symbols:
                        if splitting_works(tail, ext, res, q) and solveable_block_any(Q[i], ext, q):
                            new_prim_sol.append(solution+[(ext, res)])
            else:
                tail = solution[i][1]
                if solveable_block_any(Q[i],tail, q):
                    new_prim_sol.append(solution+[(tail, "end")])
        for solution in imprimitive_solutions:
            if i< summand_count-1:
                tail=solution[i][1]
                for ext in allowed_symbols:
                    for res in allowed_symbols:
                        if splitting_works(tail, ext, res, q) and solveable_block_prim(Q[i], ext, q):
                            new_prim_sol.append(solution+[(ext, res)])
                        elif splitting_works(tail, ext, res, q) and solveable_block_any(Q[i], ext, q):
                            new_imprim_sol.append(solution+[(ext, res)])
            else:
                tail = solution[i][1]
                if solveable_block_prim(Q[i], tail, q):
                    new_prim_sol.append(solution+[(ext, res)])
        primitive_solutions = new_prim_sol
        imprimitive_solutions = new_imprim_sol
    return primitive_solutions[0]

def solve_block(block, t, q):
    if len(block)==1:
        return [solution_type_I(t, block[0], q)]
    else:
        return solution_type_II(t, block[0], block[1], block[2], q)

def find_modq_rep_sum(Q, typelist, t, q):
    n = len(Q)
    ans = list()
    curr = t
    assert symbol_of(curr, q) == typelist[0][1]
    for i in range(1, n+1):
        if i<n:
            here,curr = find_split(curr, typelist[i][0], typelist[i][1], q)
        else:
            here = curr
        ans+=solve_block(Q[i-1], here, q)
    return vector(ZZ, ans)

# This tests above through random generation again
def test_find_modq_representation(q, n,coeffsize):
    p = prime_factors(q)[0]
    # Form random A
    A = Matrix(ZZ, n, n,0)
    for i in range(0, n):
        A[i,i]=randrange(0, coeffsize)
        for j in range(i+1, n):
            A[i,j]=randrange(0, coeffsize)
            A[j, i]=A[i,j]
    # Random vector
    v = vector(ZZ, n)
    for i in range(0, n):
        v[i] = randrange(0, coeffsize)
    if v[0] % p == 0:
        v[0] = 1
    # Compute result
    t = v.inner_product(A*v)
    t = t %q
    ans = find_modq_representation(A, t, q)
    assert Mod(ans.inner_product(A*ans),q)==Mod(t, q)

# Now some CRT to deal with non prime power q
def p_part(m, p):
    return p**(valuation(m, p))

def find_modular_representation(A, t, m): #THis needs testing
    primes = prime_factors(m)
    residues = list()
    moduli = list()
    n=A.dimensions()[0]
    retval = vector(ZZ, n)
    for i in range(0,n):
        residues.insert(i, list())
        moduli.insert(i,list())
    for p in primes:
        vec = find_modq_representation(A, t, p_part(m, p))
        for i in range(0, n):
            residues[i].append(vec[i]%(p_part(m,p)))
            moduli[i].append(p_part(m,p))
    for i in range(0, n):
        retval[i]=crt(residues[i], moduli[i])
    assert retval.inner_product(A*retval)%m==t%m
    return retval

def test_modular_representation(n, modulus, coeffsize):
    # Form random A
    A = Matrix(ZZ, n, n,0)
    for i in range(0, n):
        A[i,i]=randrange(0, coeffsize)
        for j in range(i+1, n):
            A[i,j]=randrange(0, coeffsize)
            A[j, i]=A[i,j]
    # Random vector
    v = vector(ZZ, n)
    for i in range(0, n):
        v[i] = randrange(0, coeffsize)
    v[0]=1 #Give right gcd
    # Compute result
    t = v.inner_product(A*v)
    t = t %modulus
    ans = find_modular_representation(A, t, modulus)
    assert Mod(ans.inner_product(A*ans),modulus)==Mod(t, modulus)

# Normalization of quadratic forms
# Think about strategy here

# P-adic normalization is mostly diagonalization followed by sorting by p-power (might be done?)
# and then combining the nonresidue bits in each scale compartment (bit tricky)
# We need some helpers for local transformation

def distinguished_nonsquare(p):
    P = Primes()
    i=2
    while kronecker_symbol(i, p) == 1:
        i = P.next(i)
    return i

def odd_p_normalization_from_diagonal(D,q):
    ''' Returns U, T st. U.transpose()*D*U=T, and T is in the p-adic normal form'''
    #We can just sweep down and write in the U?
    # Or are they best as compositions
    # The sorting by scales is done
    p = prime_factors(q)[0]
    L = Matrix(IntegerModRing(q), D) #local copy to check at end
    n = D.dimensions()[0]
    U = Matrix(IntegerModRing(q), n,n, 1)
    #First normalize the diagonal
    Ut = Matrix(IntegerModRing(q), n, n, 0)
    for i in range(0, n):
        num=L[i,i].lift()
        scale=valuation(num, p)
        red = num/(p**scale)
        if kronecker_symbol(red,p) == -1:
            red = Mod(red,q)/Mod(distinguished_nonsquare(p), q)
        Ut[i,i]=sqrt_mod_q(Mod(red, q)**(-1),q)
    L=Ut.transpose()*L*Ut
    U=U*Ut
    #Not done yet: need to walk things down the scale
    for i in range(0, n-1):
        num = L[i,i].lift()
        scale = valuation(num, p)
        nnum = L[i+1, i+1].lift()
        nscale = valuation(nnum, p)
        Ut = Matrix(IntegerModRing(q), n, n, 1)
        if scale == nscale:
            rednum = num/(p**scale)
            rednnum = nnum/(p**scale)
            if rednum == distinguished_nonsquare(p):
                if rednnum == 1:
                    Ut[i, i] = 0
                    Ut[i+1, i+1] = 0
                    Ut[i, i+1] = 1
                    Ut[i+1, i] = 1
                else:
                    A = odd_q_normalization_helper(q)
                    Ut[i,i] = A[0,0]
                    Ut[i, i+1] = A[0, 1]
                    Ut[i+1, i] = A[1,0]
                    Ut[i+1, i+1] = A[1,1]
        U=U*Ut
        L = Ut.transpose()*L*Ut
        assert gcd(Ut.determinant().lift(),q)==1
    return U, L

def odd_q_normalization_helper(q):
    # solution: (a^2+c^2)*r=1
    # b = c
    # d = -a
    #standard 2x2 matrix equation
    p = prime_factors(q)[0]
    A = Matrix(ZZ, 2, 2, [[distinguished_nonsquare(p), 0], [0, distinguished_nonsquare(p)]])
    (a, c) = find_modq_representation(A, 1, q)
    return Matrix(IntegerModRing(q), 2, 2, [a, c, c, -a])

# 2-adic case is much more complex
# Helper functions:

def twoadic_adjacent_blocks(A, ind1, ind2):
    pass

def twoadic_blocktype(A, ind):
    pass

def twoadic_same_scale(A, ind1, ind2):
    pass

# Normalize Type I and Type II blocks
def normalization_factor_type_I(block,q):
    val = valuation(block, 2)
    res = block/2**val
    typus = res % 8
    fixer = typus/res
    return sqrt_mod_q(Mod(fixer,q),q)

#Nota bena, uses matrix entries not the form coefficients
def normalization_matrix_type_II(a,b,c,q):
    scale = valuation(b, 2)
    a = a/(2**scale)
    b = b/(2**scale)
    c = c/(2**scale)
    s = q
    q = 2*q
    ourmat = Matrix(IntegerModRing(q), 2, 2, [a, b, b, c])
    lamb = ourmat.determinant().lift()%8
    vec = find_modq_representation(ourmat.lift(), 2, q)
    if vec[0] %2 == 0:
        T = Matrix(IntegerModRing(q), 2, 2, [[0, 1], [1, 0]])
        x1 = vec[1]
        x2 = vec[0]
    else:
        T = Matrix(IntegerModRing(q), 2, 2, 1)
        x1 = vec[0]
        x2 = vec[1]
    U = Matrix(IntegerModRing(q), 2, 2, [x1, 0, x2, Mod(1/x1,q)])
    G=T*U
    ndet=(G.transpose()*ourmat*G).determinant()
    xfix = sqrt_mod_q(lamb/ndet,q)
    V = Matrix(IntegerModRing(q), 2, 2, [1, 0, 0, xfix])
    G = G*V
    S = G.transpose()*ourmat*G
    W = Matrix(IntegerModRing(q), 2, 2, [1, Mod((1-S[0,1].lift())/2, q),0, 1])
    G = G*W
    checkmat=G.transpose()*ourmat*G
    checkmat=Matrix(IntegerModRing(s), matlift(checkmat))
    assert checkmat.determinant().lift() %8 == ourmat.determinant().lift()%8
    assert isplusminus(checkmat)
    return G

def isplusminus(mat):
    if mat[0,0]==2 and mat[1,0]==1 and mat[0,1]==1:
        return (mat[1,1] == 4 or mat[1,1]==2)
    else:
        return False

def matlift(mat):
    M=Matrix(ZZ, mat.dimensions()[0], mat.dimensions()[1],0)
    for i in range(0,mat.dimensions()[0]):
        for j in range(0,mat.dimensions()[1]):
            M[i,j]=mat[i,j].lift()
    return M
#TODO: write the global part that normalizes all blocks

# Break scales to either by all type I or all type II

# First the local part where we break a type I+type II sum
def find_local_scalebreak(tau, a, b, c, q):
    scale = valuation(tau, 2)
    tau = tau/(2**scale)
    a = a/(2**scale)
    b = b/(2**scale)
    c = c/(2**scale)
    # we should check block formedness. Assume normalized
    ourmat = Matrix(IntegerModRing(q), 3,3, [tau, 0, 0, 0, a, b, 0, b, c])
    oldmat = ourmat
    U=Matrix(IntegerModRing(q), 3, 3, 0)
    if c==4:
        V=normalization_matrix_type_II(8,1,2, q)
        V=V.inverse()
        U=Matrix(IntegerModRing(q), 3, 3, [1,0,0, 0, V[0,0], V[0,1], 0, V[1,0], V[1,1]])
    else:
        U=Matrix(IntegerModRing(q), 3, 3, 1)
    ourmat = U.transpose()*ourmat*U
    r = -ourmat[0,0]/(ourmat[1,1]+1)
    V = Matrix(IntegerModRing(q), 3, 3, [[1,1,1], [r, 1, 0], [0,1,1]])
    G = U*V
    ourmat = V.transpose()*ourmat*V
    W = locally_block_split_2adic(matlift(ourmat),q)[0]
    W = Matrix(IntegerModRing(q), W)
    G = G*W
    ourmat = W.transpose()*ourmat*W
    # We don't need to normalize split pieces as compartment will take care of it
    return G

# And now the global part that applies above repeatedly to a scale


# Sign walk
# We only need to put sign into or out of a compartment at the ends as canonicalization step
# will canonicalize. We have to manually track the total sign of a compartment as we might
# be forceably negating it

# Is this enough? Eeep....
# Does Pall/Jones give better answer?

# We use the entries, not 1/2 times the entries
# Again, this subject

#Each of these pulls the sign downwards in scale towards the front of the train

# Note that we know the inputs and outputs, so should be able to have strong assertions
# of correctness
def signwalk_empty_couple(tau1, tau2, q):
    # Todo more tests
    tau1 = tau1.lift()
    tau2 = tau2.lift()
    ourform = Matrix(IntegerModRing(q), 2, 2, [tau1, 0, 0, tau2])
    oldform = ourform
    scale = valuation(tau1, 2)
    tau1r = tau1/(2**scale)
    tau2r = tau2/(2**(scale+2))
    U = Matrix(IntegerModRing(q),2,2, [1, 4, 1, Mod(-tau1r/tau2r, q/4).lift()])
    ourform = U.transpose()*ourform*U
    a= normalization_factor_type_I(ourform[0,0].lift(),q)
    b =normalization_factor_type_I(ourform[1,1].lift(),q)
    V = Matrix(IntegerModRing(q),2, 2, [a, 0, 0, b])
    ourform=V.transpose()*ourform*V
    return U*V

def signwalk_I_II(tau1, a, b, c, q): #Here the type II block must be one scale up
    ourform = Matrix(IntegerModRing(q), 3, 3, [tau1, 0, 0, 0, a, b, 0, b, c])
    V = Matrix(IntegerModRing(q), 3, 3, [1, 0, 0, 1, 1, 0, 0, 0, 1])
    newform = V.transpose()*ourform*V
    W=locally_block_split(matlift(newform),q)[0]
    W = Matrix(IntegerModRing(q), W)
    newform = W.transpose()*newform*W
    # and now normalize newform
    a = normalization_factor_type_I(newform[0,0].lift(),q)
    A = Matrix(IntegerModRing(q), 1, 1, a)
    B = normalization_matrix_type_II(newform[1,1].lift(), newform[1,2].lift(), newform[2,2].lift(), q)
    X = A.block_sum(B)
    G = V*W*X
    return G

def signwalk_II_I(a,b,c,tau1,q): #And here one scale back
    ourform = Matrix(IntegerModRing(q), 3, 3, [a, b, 0, b, c, 0, 0, 0, tau1])
    V = Matrix(IntegerModRing(q), 3, 3, [1, 0, 0, 0, 1, 0, 0, 1, 1])
    newform = V.transpose()*ourform*V
    W = locally_block_split(matlift(newform), q)[0]
    W = Matrix(IntegerModRing(q), W)
    newform = W.transpose()*newform*W
    A = normalization_matrix_type_II(newform[0, 0].lift(), newform[0,1].lift(), newform[1,1].lift(), q)
    b = normalization_factor_type_I(newform[2,2].lift(),q)
    B = Matrix(IntegerModRing(q), 1, 1, b)
    X = A.block_sum(B)
    G = V*W*X
    return G

def signwalk_II_II(a1, b1, c1, a2, b2,c2, q): #Only in same scale.
    ourform = Matrix(IntegerModRing(q), 4, 4, [[a1, b1, 0, 0],
                                               [b1, c1, 0, 0],
                                               [0,  0, a2, b2],
                                               [0,  0, b2, c2]])
    G = Matrix(IntegerModRing(q), 4, 4, 1)
    if c1 == 4:
        G = Matrix(IntegerModRing(q), 4, 4, [[0,0,1,0],
                                             [0,0,0,1],
                                             [1,0,0,0],
                                             [0,1,0,0]])
    else:
        V= Matrix(IntegerModRing(q), 4, 4, [[1,0, 0,0],
                                            [0,1, 0,0],
                                            [0,1, 1,0],
                                            [0,0, 0,1]])
        newform = V.transpose()*ourform*V
        W = locally_block_split(matlift(newform), q)[0]
        W = Matrix(IntegerModRing(q), W)
        newform = W.transpose()*newform*W
        X = Matrix(IntegerModRing(q), 2, 2, 1)
        B = normalization_matrix_type_II(newform[2,2].lift(), newform[2,3].lift(), newform[3,3].lift(), q)
        X = X.block_sum(B)
        G = V*W*X
    return G

def block_sign_II(a,b,c,q):
    a = a.lift()
    b = b.lift()
    c = c.lift()
    scale = valuation(b, 2)
    fac = 2**scale
    a = a/fac
    b = b/fac
    c = c/fac
    det = a*c-b**2
    return kronecker_symbol(det, 2)


def block_sign_I(block, q):
    block = block.lift()
    scale = valuation(block, 2)
    return kronecker_symbol(block/(2**scale), 2)

# And now the global application of above

# Canonicalize compartments (by pulling out representible numbers and diagonalizing)
# with tricky case in n=3
# Now, we didn't sign walk all the way: these can are at most even. But that's ok! We're
# injecting or removing signs in this step anyway, and everything works out if Conway is OK

def extend_vec(x, q):
    primes = prime_factors(q)
    n = x.length()
    vecs = list()
    moduli = list()
    for p in primes:
        A = extend_vec_power(x, p_part(q, p))
        prob = vector(ZZ,n*n)
        for i in range(0,n):
            for j in range(0,n):
                prob[n*i+j] = A[i,j]
        vecs.append(prob)
        moduli.append(p**(valuation(q, p)))
    retval = Matrix(IntegerModRing(q), n, n, 0)
    retvec = sage.arith.misc.CRT_vectors(vecs, moduli)
    for i in range(0, n):
        for j in range(0, n):
            retval[i,j] = retvec[n*i+j]
    return Matrix(IntegerModRing(q),retval)
    
def extend_vec_power(x, q):
    n = x.length()
    ind = False
    good = False
    U = Matrix(IntegerModRing(q), n, n, 0)
    for i in range(0, n):
        if gcd(x[i], q)==1:
            ind = i
            good = True
            break
    assert good
    U[0,0] = x[ind]
    for i in range(1, ind):
        U[1,0] = x[i]
    U[ind, 0] = x[0]
    for i in range(ind+1, n):
        U[i, 0] = x[i]
    xinv = Mod(x[ind], q)**(-1)
    for i in range(1, n):
        U[i,i]=xinv
    U.swap_rows(0, ind)
    assert U.determinant()!=0
    return U

def lemma_19_application(tau1, tau2, tau3, q):#All rendered odd
    # This gets very caseworky: we need to look up the oddity and sign in a table
    # extract the thing.
    ourform = Matrix(IntegerModRing(q), 3, 3, [[tau1, 0, 0], [0, tau2,0],[0,0, tau3]])
    oddity = (tau1+tau2+tau3)%8
    sign = kronecker_symbol(tau1*tau2*tau3,2)
    sign_index = (sign+1)/2
    oddity_index = (oddity-1)/2
    oddity_sign_table = [[3,1,1,1], [7,1,3,1]]
    target = oddity_sign_table[sign_index][oddity_index]
    x = find_modq_representation(matlift(ourform), target, q)
    # Need to extend x to something in GL_3
    U = extend_vec(x, q)
    newform = U.transpose()*ourform*U
    # Now we need to canonicalize the lower corner of the block
    X = locally_block_split(matlift(newform), q)[0]
    X = Matrix(IntegerModRing(q), X)
    newform = X.transpose()*newform*X
    T = Matrix(IntegerModRing(q), 1,1,1).block_sum(normal_form(newform[1:, 1:],q)[0])
    newform = T.transpose()*newform*T
    G = Matrix(IntegerModRing(q), 3, 3, 0)
    # Need to sort the diagonal with G
    diag = [newform[0,0].lift(), newform[1,1].lift(), newform[2,2].lift()]
    perm = numpy.argsort(diag)
    for i in range(0, 3):
        G[perm[i], i] = 1
    ret = U*X*T*G
    return ret
    
# And now some more bookkeeping to canonicalize a compartment
# TODO: optimization, profiling
def canonicalize_compartment(M, q):
    n = M.dimensions()[0]
    k = valuation(q, 2)
    signs = [1, 3, 5, 7]
    target = False
    scale = valuation(M[0,0].lift(), 2)
    if n == 1:
        return Matrix(IntegerModRing(q), 1, 1, normalization_factor_type_I(M[0,0].lift(),q))
    elif n==2:
        found = False
        for exp in range(0, k):
            for sign in signs:
                if is_modq_rep_type(matlift(M), (exp, sign), q):
                    target = (exp, sign)
                    found = True
                    break
            if found:
                break
        targnum = (2**target[0])*target[1]
        vec = find_modq_representation(matlift(M), targnum, q)
        U = extend_vec(vec, q)
        nform = U.transpose()*M*U
        V = locally_block_split(matlift(nform), q)[0]
        V = Matrix(IntegerModRing(q), V)
        G = U*V
        nform = V.transpose()*nform*V
        L = Matrix(IntegerModRing(q), 2, 2, [normalization_factor_type_I(nform[0,0].lift(),q), 0,
                                             0, normalization_factor_type_I(nform[1,1].lift(),q)])
        G = G*L
        return G
    elif n==3:
        # The paper says we extract smallest and check if we need lemma 19.
        # I don't like that: we know going in if all elements are same scale, in which
        # case applying lemma 19 will work.
        a1 = M[0,0].lift()
        a2 = M[1,1].lift()
        a3 = M[2,2].lift()
        scale1 = valuation(a1, 2)
        scale2 = valuation(a2, 2)
        scale3 = valuation(a3, 2)
        G = Matrix(IntegerModRing(q), 3, 3, 1)
        if scale1==scale2 and scale2 == scale3:
            G = lemma_19_application(a1/(2**scale1), a2/(2**scale2), a3/(2**scale3), q)
        else:
            target = False
            found = False
            for exp in range(scale, k):
                for sign in signs:
                    if is_modq_rep_type(matlift(M), (exp, sign), q):
                        target = (exp, sign)
                        found = True
                        break
                if found:
                    break
            targnum = (2**target[0])*target[1]
            vec = find_modq_representation(matlift(M), targnum, q)
            U = extend_vec(vec, q)
            nform = U.transpose()*M*U
            V = locally_block_split(matlift(nform), q)[0]
            V = Matrix(IntegerModRing(q), V)
            G = U*V
            nform = V.transpose()*nform*V
            summand = nform[1:, 1:]
            L = Matrix(IntegerModRing(q), 1, 1, normalization_factor_type_I(nform[0,0].lift(),q))
            L = L.block_sum(canonicalize_compartment(summand,q))
            G = G*L
        return G
    else:
        funny = False
        i1 = valuation(M[0,0].lift(),2)
        i2 = valuation(M[1,1].lift(),2)
        i3 = valuation(M[2,2].lift(),2)
        i4 = valuation(M[3,3].lift(),2)
        t1 = M[0,0].lift()/(2**i1)
        t2 = M[1,1].lift()/(2**i2)
        t3 = M[2,2].lift()/(2**i3)
        ourform = M
        G = Matrix(IntegerModRing(q), n, n, 1)
        if i2-i1 ==0  and i3-i1==0 and i4-i1==1:
            funny = true
        if funny and (t1+t2+t3)%8 not in [3, 7]:
            W = Matrix(IntegerModRing(q), 2, 2, 1).block_sum(Matrix(IntegerModRing(q), 2, 2,
                                                                    [1, 0, 1, 1]))
            W = W.block_sum(Matrix(IntegerModRing(q), n-4, n-4, 1))
            ourform = W.transpose()*M*W
            G = G*W
            Z = locally_block_split(matlift(ourform[2:4, 2:4]),q)[0]
            Z = Matrix(IntegerModRing(q), Z)
            Z = Matrix(IntegerModRing(q), 2, 2, 1).block_sum(Z).block_sum(Matrix(IntegerModRing(q), n-4, n-4, 1))
            G = G*Z
            ourform = Z.transpose()*ourform*Z
            W = global_block_normalize(ourform, q)
            G = G*W
            ourform = W.transpose()*ourform*W
        target = False
        for exp in range(scale, k):
            if is_modq_rep_type(matlift(M), (exp, 1), q):
                target = (exp, 1)
                break
        targnum = (2**target[0])
        if funny:
            veci = find_modq_representation(matlift(ourform)[0:3, 0:3], targnum, q)
            vec = vector(ZZ, n)
            for i in range(0, 3):
                vec[i] = veci[i]
        else:
            vec = find_modq_representation(matlift(ourform), targnum, q)
        U = extend_vec(vec, q)
        nform = U.transpose()*ourform*U
        G = G*U
        V = locally_block_split(matlift(nform), q)[0]
        V = Matrix(IntegerModRing(q), V)
        G = G*V
        nform = V.transpose()*nform*V
        summand = nform[1:, 1:]
        # Todo: remove Type II blocks that may appear in summand
        L = Matrix(IntegerModRing(q), 1, 1, normalization_factor_type_I(nform[0,0].lift(),q))
        L = L.block_sum(global_scalebreak(summand,q))
        G = G*L
        nform = L.transpose()*nform*L
        summand = nform[1:, 1:]
        L = Matrix(IntegerModRing(q), 1, 1, normalization_factor_type_I(nform[0,0].lift(),q))
        L = L.block_sum(canonicalize_compartment(summand,q))
        G = G*L
        return G
# TODO: make  a test for this that creates random compartment and what?
# Oh, we can use oddity and scales to make isomorphic ones
# And then see they get to same thing
# Now we put everything together

#Global in the whole matrix sense
def global_block_normalize(M,q):
    # Use a block matrix already
    n = M.dimensions()[0]
    G = Matrix(IntegerModRing(q), 0, 0, 1)
    i = 0
    while i<n:
        if i<n-1 and M[i, i+1] != 0:
            G=G.block_sum(normalization_matrix_type_II(M[i,i].lift(), M[i, i+1].lift(),
                                                       M[i+1, i+1].lift(),q))
            i+=2
        else:
            G=G.block_sum(Matrix(IntegerModRing(q), 1, 1, normalization_factor_type_I(M[i,i].lift(), q)))
            i+=1
    assert is_block_normalized(G.transpose()*M*G, q)
    return G

def global_scalebreak(M,q):
    n = M.dimensions()[0]
    O = M
    G = Matrix(IntegerModRing(q), n, n, 1)
    i = 0
    while i<n-2: #Why is this enough? Because at i=n-2 we either have two type I or one type II
        if O[i,i+1] != 0:
            i+=2
        elif O[i+1, i+2] != 0:
            scale1 = valuation(O[i,i].lift(),2)
            scale2 = valuation(O[i+1, i+2].lift(),2)
            if scale1==scale2:
                T = Matrix(IntegerModRing(q), i, i, 1)
                T = T.block_sum(find_local_scalebreak(O[i,i].lift(),
                                                      O[i+1, i+1].lift(),
                                                      O[i+1, i+2].lift(),
                                                      O[i+2, i+2].lift(), q))
                x = T.dimensions()[0]
                T = T.block_sum(Matrix(IntegerModRing(q), n-x, n-x,1))
                G = G*T
                O = T.transpose()*O*T
            i+=1
        else:
            i+=1
    assert is_scale_split(G.transpose()*M*G, q)
    return G

def scale_indices(M,q):
    scalelist = sum_of_scales(M,q)
    idx = 0
    indexlist = list()
    for scale in scalelist:
        indexlist.append((idx, scale[0], scale[2]))
        idx += scale[1]
    return indexlist

#We don't represent entire trains: indexlist gives us enough to work with here.

# Sum of scale form is a list e_i, n_i, type,
# Not a genus symbol yet: we need to do more work before
# we can say what it is... well, could avoid some of it
# also doesn't include empty scales initially
def sum_of_scales(M,q):
    M = matlift(M)
    blocklist = form_into_subblocks(M,q) #Difference by two
    unconsolidated = list()
    scaleform = list()
    for block in blocklist:
        if len(block) == 3:
            scale = valuation(block[1], 2)-1
            typus = 2
        else:
            scale = valuation(block[0], 2)
            typus = 1
        unconsolidated.append((scale,typus))
    curr = unconsolidated[0]
    scale = curr[0]
    dim = curr[1]
    typus = curr[1]
    unconsolidated = unconsolidated[1:]
    for piece in unconsolidated:
        if scale != piece[0]:
            scaleform.append((scale, dim, typus))
            dim = piece[1]
            typus = piece[1]
            scale = piece[0]
        else:
            dim += piece[1]
    scaleform.append((scale, dim, typus))
    return scaleform


def compartment_list(M,q):
    scales = scale_indices(M, q)
    compartments = list()
    n = len(scales)
    i = 0
    exp = -1
    in_compartment = False
    while i<n:
        if in_compartment == False:
            if scales[i][2] == 1:
                start = scales[i][0]
                exp = scales[i][1]
                in_compartment = True
            i+=1
        else:
            if scales[i][2] == 2 or exp != scales[i][1]-1:
                compartments.append((start, scales[i][0]))
                in_compartment = False
            else:
                exp = scales[i][1]
                i+=1
    if in_compartment:
        compartments.append((start, M.dimensions()[0]))
    return compartments

def compartment_sign(M, compartment, q):
    start = compartment[0]
    end = compartment[1]
    prod = 1
    for i in range(start, end):
        a = M[i,i].lift()
        val = valuation(a,2)
        aodd = a/(2**val)
        prod *= aodd
    return kronecker_symbol(prod, 2)


def compartment_inverse_lookup(n, complist):
    ret = vector(ZZ, n)
    ind = len(complist)
    for j in range(0, n):
        ret[j]=-1
    for i in range(0, ind):
        comp = complist[i]
        for j in range(comp[0], comp[1]):
            ret[j]=i
    return ret

def global_sign_walk(M,q):
    # This gets nastily complicated. We want to walk signs, but need to deal with
    # inside scales (only for type II)
    n = M.dimensions()[0]
    compartments = compartment_list(M,q)
    # Rev compartments lets us know when we are in compartments
    rev_compartments = compartment_inverse_lookup(n, compartments)
    inputform = M #We change M
    G = Matrix(IntegerModRing(q), n, n, 1)
    # Think about iterating scales vs walking blocks
    i=n-1
    while i>=0:
        # Determine if we have a minus sign that can be pushed back, push it back
        # then step back
        if rev_compartments[i] != -1:
            comp = compartments[rev_compartments[i]]
            i = comp[0]
            if compartment_sign(M, comp, q)==-1:
                    if i>1 and M[i-2, i-1] != 0:
                        #Do a II_I signwalk
                        scalei = valuation(M[i,i].lift(),2)
                        scaleii = valuation(M[i-2, i-1].lift(),2)
                        if scaleii==scalei-1:
                            T=Matrix(IntegerModRing(q), n, n, 1)
                            T[i-2:i+1, i-2:i+1]=signwalk_II_I(M[i-2,i-2], M[i-2,i-1], M[i-1,i-1],
                                                           M[i,i], q)
                            G=G*T
                            M = T.transpose()*M*T
                    elif i>0:
                        scalei=valuation(M[i,i].lift(),2)
                        scaleim = valuation(M[i-1, i-1].lift(),2)
                        if scaleim==scalei-2:
                            T=Matrix(IntegerModRing(q), n,n,1)
                            T[i-1:i+1, i-1:i+1]=signwalk_empty_couple(M[i-1,i-1], M[i,i],q)
                            G=G*T
                            M = T.transpose()*M*T
            i-=1
        else:
            #We are in a type II block, approaching from the end
            # Below shouldn't underflow
            i-=1
            #Determine sign, and walk it back if necessary/possible
            if block_sign_II(M[i,i], M[i, i+1], M[i+1, i+1], q)==-1:
                if i>1 and M[i-1, i-2] != 0:
                    scalei = valuation(M[i,i+1].lift(),2)
                    scaleim = valuation(M[i-2, i-1].lift(),2)
                    if scalei==scaleim:
                        T=Matrix(IntegerModRing(q), n, n, 1)
                        T[i-2:i+2, i-2:i+2] = signwalk_II_II(M[i-2, i-2], M[i-2, i-1], M[i-1,i-1], M[i, i], M[i, i+1], M[i+1, i+1], q)
                        G = G*T
                        M = T.transpose()*M*T
                elif i>0:
                    # Maybe a type I is behind us
                    scalei=valuation(M[i,i+1].lift(),2)
                    scaleim = valuation(M[i-1, i-1].lift(),2)
                    if scaleim == scalei-1:
                        T=Matrix(IntegerModRing(q), n, n, 1)
                        T[i-1:i+2, i-1:i+2] = signwalk_I_II(M[i-1, i-1], M[i,i],
                                                            M[i,i+1], M[i+1,i+1],q)
                        G=G*T
                        M=T.transpose()*M*T
            i-=1
    return G

def global_canonicalize_compartments(M,q):
    compartments = compartment_list(M,q)
    n = M.dimensions()[0]
    G = Matrix(IntegerModRing(q), n, n, 1)
    for compartment in compartments:
        start, end = compartment
        G[start:end, start:end]=canonicalize_compartment(M[start:end, start:end],q)
    return G

def even_q_normal_form(M,q):
    G1=locally_block_split(matlift(M), q)[0]
    G1 = Matrix(IntegerModRing(q), G1)
    M = G1.transpose()*M*G1
    assert is_block_diagonalized(M,q)
    G2 = global_block_normalize(M, q)
    M = G2.transpose()*M*G2
    assert is_block_normalized(M,q)
    G3 = global_scalebreak(M,q)
    M = G3.transpose()*M*G3
    assert is_scale_split(M,q)
    G4 = global_sign_walk(M, q)
    M = G4.transpose()*M*G4
    G5 = global_canonicalize_compartments(M,q)
    M = G5.transpose()*M*G5
    assert is_block_normalized(M,q)
    return G1*G2*G3*G4*G5, M

def run_counterexamples():
    M=Matrix(ZZ, 5, 5, [[  2,  11,  14,   3,  13],
                        [ 11,  74,  80,  36,  94],
                        [ 14,  80,  84, 106,  42],
                        [  3,  36, 106,  47,  75],
                        [ 13,  94,  42,  75,  73]])
    q=2^12
    even_q_normal_form(Matrix(IntegerModRing(q), M), q)
    q = 8
    M=Matrix(ZZ, 4, 4, [[  3, 0, 0 ,0],
                        [ 0, 7, 0, 0],
                        [ 0, 0, 7, 0],
                        [ 0, 0, 0, 1]])
    even_q_normal_form(Matrix(IntegerModRing(q), M), q)
def odd_q_normal_form(M,q):
    G1=locally_block_split(matlift(M), q)[0]
    G1= Matrix(IntegerModRing(q), G1)
    M = G1.transpose()*M*G1
    G2 = odd_p_normalization_from_diagonal(M, q)[0]
    M = G2.transpose()*M*G2
    return G1*G2, M

def normal_form(M, q): # Only works modulo prime powers
    if q %2 == 0:
        return even_q_normal_form(M,q)
    else:
        return odd_q_normal_form(M,q)

#The random testing facility
def diagonal_form(q, entries):
    n = len(entries)
    M = Matrix(IntegerModRing(q), n, n, 1)
    for i in range(0, n):
        M[i, i] = entries[i]
    return M

def is_sign_walked(M,q):
    pass

def is_block_diagonalized(M,q):
    n = M.dimensions()[0]
    for i in range(0, n):
        for j in range(i+2, n):
            if M[i,j] != 0:
                return False
    return True

def random_transformation_test(M, coeffsize, q):
    # Apply a random element of SL_n(q) and see if we get same out
    # Issue with determinant is why SL_n(q)
    n = M.dimensions()[0]
    L = Matrix(IntegerModRing(q), n, n, 1)
    U = Matrix(IntegerModRing(q), n, n, 1)
    for i in range(0, n):
        for j in range(i+1, n):
            L[j, i] = randrange(0, coeffsize)
            U[i, j] = randrange(0, coeffsize)
    P = Matrix(IntegerModRing(q), n, n,1)
    for i in range(0, n):
        j = randrange(0, n)
        P.swap_rows(i,j)
    G = L*U*P
    N=G.transpose()*M*G
    Mpr=even_q_normal_form(M,q)[1]
    Npr=even_q_normal_form(N, q)[1]
    if Mpr != Npr:
        return N, M
    else:
        return True

def is_block_normalized(M,q):
    # Should also assert scales increase
    if not is_block_diagonalized(M,q):
        return False
    n = M.dimensions()[0]
    i = 0
    signs = [1,3,5,7]
    normalized = True
    while i<n:
        if i<n-1 and M[i,i+1]!=0:
            b = M[i,i+1].lift()
            scale = valuation(b, 2)
            power = 2**scale
            at = (M[i, i].lift())/power
            bt = b/power
            ct = (M[i+1, i+1].lift())/power
            if at != 2:
                 normalized = False
            if bt != 1:
                normalized = False
            if ct != 2 and ct !=4:
                 normalized = False
            i+=2
        else:
            a = M[i,i].lift()
            scale = valuation(a, 2)
            anorm = a/(2**scale)
            if not (anorm in signs):
                normalized = False
            i+=2
    return normalized

def is_scale_split(M,q):
    M = matlift(M)
    blocklist = form_into_subblocks(M,q) #Difference by two
    unconsolidated = list()
    for block in blocklist:
        if len(block) == 3:
            scale = valuation(block[1], 2)-1
            typus = 2
        else:
            scale = valuation(block[0], 2)
            typus = 1
        unconsolidated.append((scale,typus))
    scale = -1
    typus = 1
    for form in unconsolidated:
        if form[0] != scale:
            scale = form[0]
            typus = form[1]
        else:
            if form[1] != typus:
                return False
    return True

# And now for the main course: computing forms from genus

# How do we represent a genus? Answer: a set of canonical form forms. We are only dealing with positive definite (Check)
# The issue is with expansion to other primes, so we include the determinant to make that expansion possible


# Some auxiliary helpers

def kp(p):
    if p == 2:
        return 3
    else:
        return 1

def complete_integer(q):
    prod = q
    primes = prime_factors(2*q)
    for p in primes:
        prod *= p**(kp(p))
    return prod,primes


# A genus is a list of prime, form pairs where the form is normalized
# The determinant can be computed from where the form is not unimodular
# We only do positive definite forms

def compute_genus(M):
    det = M.determinant()
    q, primelist = complete_integer(det)
    rel_p = primelist
    genus=list()
    for p in rel_p:
        mod = p_part(q, p)
        genus.append((p, normal_form(Matrix(IntegerModRing(mod), M), mod)[1]))
    return genus

# Computation of the representable t
def maj(a, b, c):
    tot = a+b+c
    if tot>=2:
        return 1
    else:
        return 0

def compute_t_4dim_helper(form):
    R= form.base_ring()
    I = R.defining_ideal()
    q = I.gens_reduced()[0]
    scales = sum_of_scales(form, q)
    # Check if we have a type 2 block, and if so what the scale is
    for scale in scales:
        if scale[2] == 2:
            return scale[0]+1
    return valuation(form[3,3].lift(), 2)

def compute_t_4dim(genus):
    #first compute r
    r = 1
    t = 1
    for pair in genus:
        prime = pair[0]
        form = pair[1]
        exp = 0
        if prime == 2:
            exp = compute_t_4dim_helper(form)
            t *= 2**exp
        else:
            i1 = valuation(form[0,0].lift(), prime)
            i2  = valuation(form[1,1].lift(), prime)
            i3 = valuation(form[2,2].lift(), prime)
            exp = maj(i1%2, i2%2, i3%2)
        r *= prime**exp
    for pair in genus:
        prime = pair[0]
        form = pair[1]
        if prime == 2: #we already do 2 completely
            continue
        i1 = valuation(form[0,0].lift(), prime)
        i2 = valuation(form[1,1].lift(), prime)
        i3 = valuation(form[2,2].lift(), prime)
        ind = 0 
        if i1 %2 == i2 %2:
            ia = i1
            ib = i2
            ind = 0
        elif i1 %2 == i3 %2:
            ia = i1
            ib = i3
            ind = 0
        else:
            ia = i2
            ib = i3
            ind = 1
        taua = form[ind, ind].lift()/(prime**ia)
        rp = r/(prime**valuation(r, prime))
        if kronecker_symbol(r, prime) == kronecker_symbol(taua, prime):
            t*=prime**ia
        else:
            t*=prime**ib
    return t

def test_4_t_compute(coeffsize):
    M=Matrix(ZZ, 4, 4, 0)
    while M.determinant() == 0:
        for i in range(0, 4):
            v=vector(ZZ, 4)
            for j in range(i, 4):
                v[j]=randrange(coeffsize)
            M+=tensor_product(v,v)
    q,_ = complete_integer(M.determinant())
    genus = compute_genus(M)
    print(genus)
    t = compute_t_4dim(genus)
    vec = find_modular_representation(M, t, q)
    assert vec.inner_product(M*vec)%q == t%q

#n=3, n=2 cases remain.

def compute_t_3dim_helper(form):
    R= form.base_ring()
    I = R.defining_ideal()
    q = I.gens_reduced()[0]
    scales = sum_of_scales(form, q)
    # Check if we have a type 2 block, and if so what the scale is
    for scale in scales:
        if scale[2] == 2:
            return (scale[0]+1, False)

    lowscale = valuation(form[0,0].lift(), 2)
    tau=form[0,0].lift()/(2**lowscale)
    return (lowscale, tau)

def compute_t_3dim(genus): #similar to 4 case, but aux prime required
        #first compute r
    r = 1
    t = 1
    needaux = False
    maxprime = 2
    twoexp = 0
    for pair in genus:
        prime = pair[0]
        if prime > maxprime:
            maxprime = prime
        form = pair[1]
        exp = 0
        if prime == 2:
            exp,needaux = compute_t_3dim_helper(form)
            t *= 2**exp
            twoexp = exp
        else:
            i1 = valuation(form[0,0].lift(), prime)
            i2  = valuation(form[1,1].lift(), prime)
            i3 = valuation(form[2,2].lift(), prime)
            exp = maj(i1%2, i2%2, i3%2)
            r *= prime**exp
    if needaux != False:
        tau = needaux
        auxprime = next_prime(maxprime)
        while Mod(r, 8)*Mod(auxprime, 8) != Mod(tau, 8):
            auxprime = next_prime(auxprime)
        t *= auxprime
        r *= auxprime
    r *=2**twoexp
    for pair in genus:
        prime = pair[0]
        form = pair[1]
        if prime == 2: #we already do 2 completely
            continue
        i1 = valuation(form[0,0].lift(), prime)
        i2 = valuation(form[1,1].lift(), prime)
        i3 = valuation(form[2,2].lift(), prime)
        ind = 0 
        if i1 %2 == i2 %2:
            ia = i1
            ib = i2
            ind = 0
        elif i1 %2 == i3 %2:
            ia = i1
            ib = i3
            ind = 0
        else:
            ia = i2
            ib = i3
            ind = 1
        taua = form[ind, ind].lift()/(prime**ia)
        rp = r/(prime**valuation(r, prime))
        if kronecker_symbol(rp, prime) == kronecker_symbol(taua, prime):
            t*=prime**ia
        else:
            t*=prime**ib
    return t

def test_3_t_compute(coeffsize):
    M=Matrix(ZZ, 3,3, 0)
    for i in range(0, 3):
        for j in range(0, 3):
            M[i,j]=randrange(coeffsize)
    M=M.transpose()+M
    q = complete_integer(M.determinant())
    genus = compute_genus(M)
    t = compute_t_3dim(genus) #We keep getting 8. Hmm.
    vec = find_modular_representation(M, t, q)
    assert vec.inner_product(M*vec)%q == t%q

def test_3_t_compute_diag(coeffsize):
    M=Matrix(ZZ, 3,3, 0)
    while M.determinant()==0:
        for i in range(0, 3):
            M[i,i]=randrange(coeffsize)
    q,_ = complete_integer(M.determinant())
    genus = compute_genus(M)
    t = compute_t_3dim(genus) #We keep getting 8. Hmm.
    vec = find_modular_representation(M, t, q)
    assert vec.inner_product(M*vec)%q == t%q

def gcd_of_genus(genus):
    prod = 1
    for species in genus:
        p = species[0]
        form = matlift(species[1])
        if p==2:
            if form[0,1]==0:
                prod *= 2**(valuation(form[0,0], 2))
            else:
                prod *= 2**(valuation(form[0,1], 2))
        else:
            prod *= p**(valuation(form[0,0], p))
    return prod

def reduce_genus(genus): #Reduce the genus by dividing each form by minimal scale
    div = gcd_of_genus(genus)
    newgens = list()
    for species in genus:
        form = species[1]
        R= form.base_ring()
        I = R.defining_ideal()
        q = I.gens_reduced()[0]
        newgens.append((species[0], normal_form(Matrix(IntegerModRing(q),matlift(form)/div), q)[1]))
    return newgens

def compute_t_2dim(genus):
    for species in genus:
        if species[0]==2:
            if species[1][0,1]!=0:
                return compute_t_2dim_type_II(genus)
            else:
                if valuation(species[1][1,1].lift(),2)%2==0:
                    return compute_t_2dim_type_I_even(genus)
                else:
                    return compute_t_2dim_type_I_odd(genus)

def compute_S_set(d,S, apmap):
    ret = set()
    for p in S:
        ap = apmap[p]
        if p%8==d:
            ret.add(p)
    return ret

def compute_S_minus(S, apmap):
    ret = set()
    for p in S:
        ap = apmap[p]
        if kronecker_symbol(ap, p)==-1:
            ret.add(p)
    return ret
        
def compute_S(genus): #Only makes sense for 2
    S = set()
    apmap = {}
    for species in genus:
        p = species[0]
        if p!=2:
            form = species[1]
            b = form[1,1].lift()
            ip = valuation(b, p)
            if ip %2 == 1:
                S.add(p)
            apmap[p] = form[0,0].lift()
    return S, apmap

def p_in_progression(a, m):
    x = a
    while not is_prime(x):
        x += m
    return x


def compute_t_2dim_type_II(genus):
    # t is 2*auxprime*r^2 where r^2 is easy
    # auxprime is in a given arithmetic progression
    S, apmap = compute_S(genus)
    S3 = compute_S_set(3, S, apmap)
    S5 = compute_S_set(5, S, apmap)
    Sm = compute_S_minus(S, apmap)
    mod4 = 1
    if len(S3)+len(S5)+len(Sm)%2==1:
        mod4 = 1
    else:
        mod4 = -1
    # Only cases because we are always positive definite (better check)
    residues = list()
    moduli = list()
    for p in S:
        residues.append(2*apmap[p])
        moduli.append(p)
    residues.append(mod4)
    moduli.append(4)
    auxa = CRT(residues, moduli)
    auxm = lcm(moduli)
    auxprime = p_in_progression(auxa, auxm)
    r = 1
    for species in genus:
        p = species[0]
        if p==2 or p in S:
            continue
        else:
            if kronecker_symbol(apmap[p], p)!=kronecker_symbol(2*auxprime, p):
                r*= p**(valuation(species[1][1,1].lift(), p)/2)
    return 2*auxprime*r^2

def compute_t_2dim_type_I_even(genus):
    # Same, a little more complex at 2
    # Woops, looks like the paper splits into two cases based on i2, and I only
    # did the even one here
    # Weirdly this works? Or maybe my test is just failing
    S, apmap = compute_S(genus)
    Sm = compute_S_minus(S, apmap)
    S3 = compute_S_set(3, S, apmap)
    S7 = compute_S_set(7, S, apmap)
    a2 = 0
    p2b2 = 0
    for species in genus:
        if species[0]==2:
            a2 = species[1][0,0].lift()
            p2b2 = species[1][1,1].lift()
    b2 = p2b2/(2**(valuation(p2b2, 2)))
    X = set([a2 %8, b2 %8])
    mod8 = 0
    if len(Sm) %2 == 0:
        if 1 in X:
            mod8 = 1
        elif 5 in X:
            mod8 = 5
    else:
        if 3 in X:
            mod8 = 3
        elif 7 in X:
            mod8 = 7
    assert mod8 != 0
    residues = list()
    moduli = list()
    for p in S:
        residues.append(apmap[p])
        moduli.append(p)
    residues.append(mod8)
    moduli.append(8)
    auxa = CRT(residues, moduli)
    auxm = lcm(moduli)
    auxp = p_in_progression(auxa, auxm)
    # Now iterate over computing r
    r = 1
    for species in genus:
        p = species[0]
        if p == 2:
            if a2 %8 != auxp %8:
                r*= 2**(valuation(p2b2,2)/2)
        elif p not in S:
            if kronecker_symbol(auxp, p) != kronecker_symbol(apmap[p],p):
                form = species[1]
                r*=p**(valuation(form[1,1].lift(), p)/2)
    return auxp*(r)**2

def compute_t_2dim_type_I_odd(genus):
    S, apmap = compute_S(genus)
    Sm = compute_S_minus(S, apmap)
    S3 = compute_S_set(3, S, apmap)
    S7 = compute_S_set(7, S, apmap)
    mulS = 1
    mod8 = 0
    twiddle = 1
    a2 = 0
    p2b2 = 0
    for species in genus:
        if species[0]==2:
            a2 = species[1][0,0].lift()
            p2b2 = species[1][1,1].lift()
            b2 = p2b2/(2**(valuation(p2b2, 2)))
            #Nasty conditional
            if a2 %4 == 1:
                if (-1)**len(Sm) == kronecker_symbol(a2, 1):
                    mod8 = a2 %8
                else:
                    twiddle = 2**(valuation(p2b2, 2))
                    mulS = 2
                    mod8 = b2 % 8
            else:
                if (-1)**(len(Sm)+len(S3)+len(S7))==kronecker_symbol(a2, 1):
                    mod8 = a2 % 8
                else:
                    twiddle = 2**(valuation(p2b2, 2))
                    mulS = 2
                    mod8 = b2 % 8
    #Find auxprime
    residues = list()
    moduli = list()
    for p in S:
        residues.append(mulS*apmap[p])
        moduli.append(p)
    residues.append(mod8)
    moduli.append(8)
    auxa = CRT(residues, moduli)
    auxm = lcm(moduli)
    auxp = p_in_progression(auxa, auxm)
    #Find r
    r = 1
    for species in genus:
        p = species[0]
        if p!= 2 and p not in S:
            if kronecker_symbol(auxp*mulS, p) != kronecker_symbol(apmap[p],p):
                form = species[1]
                r*=p**(valuation(form[1,1].lift(), p))
    return auxp*r*twiddle

def test_2_t_compute(coeffsize):
    M=Matrix(ZZ, 2,2, 0)
    # Need to generate a random positive definite rank 2 quadratic form
    # add two squares together
    while M.determinant() == 0:
        a = randrange(coeffsize)
        b = randrange(coeffsize)
        c = randrange(coeffsize)
        d = randrange(coeffsize)
        M[0,0] = a**2+c**2
        M[0,1] = a*b+c*d
        M[1, 0] = M[0,1]
        M[1,1] = b**2+d**2
    q,_= complete_integer(abs(M.determinant()))
    genus = compute_genus(M)
    gcdl = gcd_of_genus(genus)
    M = M/gcdl
    genus = compute_genus(M)
    t = compute_t_2dim(genus)
    vec = find_modular_representation(M, t, q)
    assert vec.inner_product(M*vec)%q == t%q

def test_2_t_compute_diag_odd(coeffsize):
    M = Matrix(ZZ, 2, 2, 0)
    a = 2*randrange(coeffsize)+1
    i = randrange(3)*2+1
    b = (2*randrange(coeffsize)+1)*2**i
    M[0,0]=a
    M[1,1]=b
    genus = compute_genus(M)
    q = complete_integer(abs(M.determinant()))
    t = compute_t_2dim(genus)
    vec = find_modular_representation(M, t, q)
    assert vec.inner_product(M*vec)%q == t%q
#Some more auxilliary functions

def isomorphic_form_mod_q(genus, q, primelist):
    # Enough to match each one modulo q
    # as we carefully pick the q
    # NB: depends on having expanded the genus
    genuspart = dict()
    n = genus[0][1].dimensions()[0]
    moduli = list()
    vecs = list()
    for species in genus:
        p = species[0]
        genuspart[p]=species[1]
        form = matlift(species[1])
        prob = vector(ZZ,n*n)
        for i in range(0,n):
            for j in range(0,n):
                prob[n*i+j] = form[i,j]
        vecs.append(prob)
        moduli.append(p**(valuation(q, p)))
    retval = Matrix(IntegerModRing(q), n, n, 0)
    retvec = sage.arith.misc.CRT_vectors(vecs, moduli)
    for i in range(0, n):
        for j in range(0, n):
            retval[i,j] = retvec[n*i+j]
    return retval

def determinant_of_genus(genus):
    det = 1
    for species in genus:
        p = species[0]
        form = matlift(species[1])
        locdet = form.determinant()
        det *= p**valuation(locdet, p)
    return det

def expand_genus_all_primes(genus,q,primelist):
    newgen = list(genus)
    S = set()
    for species in genus:
        S.add(species[0])
    det = determinant_of_genus(genus)
    primeset = primelist
    n = genus[0][1].dimensions()[0]
    for prime in primeset:
        if prime not in S:
            mod = p_part(q, prime)
            M=Matrix(IntegerModRing(mod), n, n, 1)
            k = det%mod
            if kronecker_symbol(k, prime)==1:
                k = 1
            else:
                 k=distinguished_nonsquare(prime)
            M[n-1, n-1]=k
            newgen.append((prime, M))
    return newgen

def compute_t(genus):
    n = genus[0][1].dimensions()[1]
    if n>=4:
        return compute_t_4dim(genus)
    elif n==3:
        return compute_t_3dim(genus)
    else:
        return compute_t_2dim(genus)
#And now the final result

def transform_mod_pk(A, B, q):
    #Return U in GLn(Q) such that A=Q.transpose()*B*Q
    norm_A, classA = normal_form(A, q)
    norm_B, classB = normal_form(B, q)
    op = norm_B*norm_A.inverse()
    assert  op.transpose()*B*op == A
    return op

def transform_mod_q(A, B,q, primelist): #CRT on above
    n = A.dimensions()[0]
    primes =primelist
    vecs = list()
    moduli = list()
    for p in primes:
        mod = p_part(q, p)
        mat = transform_mod_pk(Matrix(IntegerModRing(mod), A), Matrix(IntegerModRing(mod),B),mod)
        vec = vector(ZZ, n**2)
        for i in range(0, n):
            for j in range(0, n):
                vec[n*i+j]=mat[i,j]
        vecs.append(vec)
        moduli.append(mod)
    retvec = sage.arith.misc.CRT_vectors(vecs, moduli)
    retmat = Matrix(ZZ, n, n, 0)
    for i in range(0,n):
        for j in range(0,n):
            retmat[i,j]=retvec[n*i+j]
    if matmod(retmat.transpose()*B*retmat,q)!=matmod(A,q):
        assert False
    return retmat


def tensor_product(v, w):
    return matrix(len(v),len(w),[a*b for a in v for b in w])

def compute_inner_genus(mat, mod, aux):
    primes = prime_divisors(mod)
    gen = list()
    for p in primes:
        loc = p_part(mod, p)
        gen.append((p, normal_form(Matrix(IntegerModRing(loc), mat), loc)[1]))
    return gen

def vecmod(vec,q):
    ret = vector(ZZ, len(vec))
    for i in range(0, len(vec)):
        ret[i]=vec[i] % q
    return ret

def matmod(mat, q):
    ret = Matrix(ZZ, mat.dimensions()[0], mat.dimensions()[1])
    for i in range(0, mat.dimensions()[0]):
        for j in range(0, mat.dimensions()[1]):
            ret[i,j]=mat[i,j]%q
    return ret

def form_from_genus(genus):
    #Todo: use LLL to shrink coeffsize
    dim = genus[0][1].dimensions()[1]
    if dim == 1:
        det = determinant_of_genus(genus)
        M = Matrix(ZZ, 1, 1, det)
        return M
    #Now induct
    gcdl=gcd_of_genus(genus)
    genus = reduce_genus(genus)
    t=compute_t(genus)
    det=determinant_of_genus(genus)
    q,primelist=complete_integer(det*t**(dim-1))
    genus = expand_genus_all_primes(genus,q,primelist)
    S = isomorphic_form_mod_q(genus, q,primelist)
    S = matlift(S)
    x = find_modular_representation(S, t, q) #There are now some ring problems
    # Doesn't always work: sign of a bug in the t calculation
    R = extend_vec(x, q)
    A = R[:, 1:]
    A = matlift(A)
    d = x*S*A
    d = vecmod(d,q)
    H = t*A.transpose()*S*A-tensor_product(d,d)
    H = matmod(H, q)
    assert H.determinant() % det*t**(dim-1)==0
    inner = compute_inner_genus(H, q, det*t**(dim-2))
    Htilde = form_from_genus(inner) #Important that this is integral
    if not genus_loose_equal(inner, compute_genus(Htilde)): #Not quite right
        assert False
    #Now need to equate H and Htilde modulo q and do cleanup
    U = transform_mod_q(Htilde, H, q, primelist)
    Q = Matrix(ZZ, dim, dim, 0)
    Q[0,0]=t
    v=U.transpose()*d
    bottom = Htilde+U.transpose()*tensor_product(d,d)*U
    for i in range(1, dim):
        Q[0,i]=v[i-1]
        Q[i,0]=v[i-1]
    for i in range(1, dim):
        for j in range(1, dim):
            Q[i,j] = bottom[i-1,j-1]/t
    return reduce_LLL(Q*gcdl)

def unmake_quadratic_form(Q):
    return Matrix(ZZ, Q.Hessian_matrix()*1/2)

def reduce_LLL(M):
    Q = make_similar_quadratic_form(M)
    return unmake_quadratic_form(Q.lll())

def sage_genus_to_our_genus(genus):
    pass

def genus_equal(gen1, gen2):
    mapgen1=dict()
    mapgen2 = dict()
    for species in gen1:
        mapgen1[species[0]]=species[1]
    for species in gen2:
        mapgen2[species[0]]=species[1]
    for p in mapgen1:
        if p not in mapgen2 or  mapgen2[p] != mapgen1[p]:
            return False
    for p in mapgen2:
        if p not in mapgen1 or mapgen1[p] != mapgen2[p]:
            return False
    return True

def genus_loose_equal(gen1, gen2):
    primelist = list()
    for species in gen1:
        primelist.append(species[0])
    for species in gen2:
        primelist.append(species[0])
    q = prod(primelist)
    return genus_equal(expand_genus_all_primes(gen1, q,primelist),
                       expand_genus_all_primes(gen2, q,primelist))
    

def form_from_genus_test(dim, coeffsize):
    M=Matrix(ZZ, dim, dim, 0)
    while M.determinant() == 0:
        for i in range(0, dim):
            v=vector(ZZ, dim)
            for j in range(i, dim):
                v[j]=randrange(coeffsize)
            M+=tensor_product(v,v)
    print("Testing with the genus of")
    print(M)
    gen1 = compute_genus(M)
    form2 = form_from_genus(gen1)
    print("Right det")
    assert form2.determinant()==M.determinant()
    gen2 = compute_genus(form2)
    print("Right genus")
    assert genus_equal(gen1, gen2)
