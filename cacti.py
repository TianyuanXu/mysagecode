# To be used only when hecke.py has been loaded.


# Digression: computatzion of T_w0 * c_w.
def ts_times_cw(type,s,w):
    """ Compute T_s*c_w in the Hecke algebra.

    INPUT:
    - "type" -- the Cartan type of a Coxeter group W
    - "s" -- a number i, representing the simple reflection s=s_i
    - "w" -- a tuple, considered as an expression(not necessarily reduced) for
             an element of W.   

    OUTPUT:
    - the product T_s*c_w in the Hecke algebra of W, where T_s is the standard
      basis element for s. The output is a linear combination in the kl-basis,
      and the indexing elements are already in normal forms.

    ALGORITHM:
    Recall that T_s=c_s-v^(-1), so compute c_s*c_w-v^(-1)*c_w.

    """
    W = CoxeterGroup(type, implementation = 'coxeter3')
    w_reduced = W(list(w)).reduced_word()

    d = defaultdict()
    d = s_times_w(type,s,w)
    
    if tuple(w_reduced) in d:# Todo: avoid this; += should do but has KeyError
        d[tuple(w_reduced)] += -v**(-1)
    else:       # must be the case that sw>w and c_sc_w doesn't contain c_w
        d[tuple(w_reduced)] = -v**(-1)

    return d

def tproduct_times_cw(type,t,w):
    """ Multiply a product of c_{s}'s with c_w.

    INPUT:
    - 'type' -- the Cartan type of the Coxeter group W
    - 't' -- a tuple of simple reflections (s1,s2,...,sn)
    - 'w' -- a tuple, representing an element of W

    OUTPUT:
    - the product (T_s1 * T_s2 * ... * T_sn) * c_w in the Hecke algebra of W.
    
    """
    W = CoxeterGroup(type, implementation='coxeter3')
    d = defaultdict()
    w_reduced = W(list(w)).reduced_word()

    d[tuple(w_reduced)] = 1
    for s in reversed(t): 
        ts_times_d = defaultdict(int)
        for x in d:
            dd = ts_times_cw(type,s,x)
            for term in dd:
                ts_times_d[term] += dd[term] * d[x]
        d = remove_zero(ts_times_d)  # see below for definition of remove_zero
    return d

def tw0_times_cw(type,w):
    """ Compute T_{w0}*c_w for a finite Coxeter group.

    INPUT:
    - "type" -- the Cartan type of a finite Coxeter group W 
    - "w" -- a tuple, considered as an expression(not necessarily reduced) for 
             an element of W

    OUTPUT:
    - the product T_{w0}*c_w in the Hecke algebra of W, where w0 is the longest
      element of W

    ALGORITHM:
    Recall that if (s1,...,sn) is a reduced expression for w0, then 
        T_{w0} = T_{s1} * \cdots * T_{sn},
    so use ts_times_cw(type,t,w) with t given by any reduced word for w0.
    """

    w0 = CoxeterGroup(type).w0.reduced_word()
    d = tproduct_times_cw(type,tuple(w0),w)
    return compress_key(d)

def compress_tuple(t):
    r"""
    Compress a tuple of digits into a number, e.g., (2,3,1) --> 231
    """
    return int(''.join(map(str,t)))

def compress_key(d):
    r""" 
    Compress the keys in a dictionary $d$, returning a new dictionary.
    """
    dd = defaultdict(int)
    for k in d:
        kk = compress_tuple(k)
        dd[kk] = d[k]
    return dd

def inv(type,l,w):
    d = tw0_times_cw(type,w)
    for x in l:
        if x in d:
            return x 
        else: 
            continue

def print_inv(type,d):
    print 'The involution sigma behaves as follows in type {}:'.format(type)
    print ''
    i = 1
    for k in d:
        print 'Left cell #{}: {}'.format(i,d[k])
        for elt in d[k]:
            print '{} --> {}'.format(elt,inv(type,d[k],tuple(int(i) for i in str(elt))))
        print ''
        i = i+1

def convert_expression_BDFH(rank,w):
    return compress_tuple((rank+1-int(i) for i in str(w)))

def convert_lcells_BDFH(rank,d):
    return {k:[convert_expression_BDFH(rank,w) for w in d[k]] for k in d}


