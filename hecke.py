import itertools
from collections import defaultdict


"""
Background:
    There are two sets of conventions for Hecke algebra.
    Ours: 
    (H_s-v)(H_s+v^{-1})=0
    Theirs (the one used by Sage and Coxeter3):
    (T_s-q)(T_s+1)=0, 
    implying (T_s/q^{1/2}-q^{1/2})(T_s/q^{1/2}+q^{-1/2})=0

    So there is an algebra homomorphism \phi sending our H_s to their
    T_s/q^{1/2}, and sending our v to their q^{1/2}. Meanwhile, on each version
    of the Hecke algebras there is a "bar involution" of algebras, with ours 
    sending H_s to its inverse and v to v^{-1} and theirs sending T_s to its 
    inverse and sending q^{1/2} to q^{-1/2}. Thus, the bar involution simply 
    sends each term to its inverse in the Hecke algebra. Being an algebra
    homomorphism, \phi commutes with the bar involution. So the \phi-image of 
    our Kazhdan-Lusztig basis elements are their Kazhdan-Lusztig basis elements 
    (for Sage, the images are denoted by Cp, with C used for another basis),
    which are defined to be special bar-invariant elements in each case.

    Conversion between the two versions of Kazhdan-Lusztig(KL) polynomials:
    For us, the KL polynomials are characterized as follows:
    c_w=\sum_{y\le w} p_{y,w}(v)H_y;
    For Sage and Coxeter3, the KL polynomials are characterized as follows:
    c_w=q^{-l(w)/2}\sum_{y\le w} p_{y,w}(q)T_y.
    (Again, in Sage, this holds with Cp in place of C.)
    Using the fact that \phi(c_w)=c_w, it is then straightforward to figure out 
    how to convert the two versions of KL polynomials.
    
    From ours to theirs: 
    Mulitply p_{y,w}(v) by v^{l(w)-l(y)}, then replace each v with q^{1/2}. 
    This results in p_{y,w}(q).

    From theirs to ours:
    Multiply p_{y,w}(q) by q^{(l(y)-l(w))/2}, then replace each q with v^2.
"""


""" Things I'm not using at the moment


# variable setup
R.<q,v>=LaurentPolynomialRing(QQ,2)


# Computations using Sage's internal algorithms
# A Coxeter group can be defined using its Cartan type, or its Coxeter matrix.

m=[[1,4,3,4],[4,1,3,2],[3,3,1,2],[4,2,2,1]] 
W1=CoxeterGroup(m)


# Immediately declare s=W.simple_reflections() after defining W; 
# Indexing of s starts at 1. 
# Or declare something like [s1,s2,s3]=W.simple_reflections().



# Computations using Gap3-Cheive:
# attach pycox.py and define Coxeter groups using Cartan type or Cartan matrix.
# Computattion is faster than Sage but mostly only deals with finite groups.

"""

# Computation using coxeter3, for faster computation of KL polynomials:

var('v')
# var('q')

def coxeter3(l,rank):
    return CoxeterGroup([l,rank],implementation='coxeter3')

def kl(type,y,w):
    """ Compute a Kazhdan-Lusztig polynomial p_{y,w}.

    INPUT:
    - "type" -- the Cartan type of a Coxeter group W
    - "y", "w" --tuples, considered as expressions(not necessarily reduced) for
                 elements in W
   
    OUTPUT:
    - the Kazhdan-Lusztig polynomial p_{y,w}, returned in our convention. 


    .. TODO::
    - use Coxeter matrix instead of Cartan type to specify Coxeter groups.
    """
   
    W = CoxeterGroup(type,implementation='coxeter3')
    f = W.kazhdan_lusztig_polynomial(y,w)
    g = f.substitute(q=v**2)*v**(len(y)-len(w))
    return g.expand()

def mu(type,y,w):
    """ Compute the mu_{y,w}, the mu-coefficient for a pair of elements y,w.
    
    .. TODO::
        Use the coxeter3 package to directly compute mu_{y,w} instead of 
        reading it off as a coefficient of p_{y,w}.
        The direct computation can be---and is---done in s_times_w, so this
        function is actually not used.
    """

    return kl(type,y,w).coefficient(v,-1)
    
def s_times_w(type,s,w):
    """ Compute the product c_s*c_w, ere s is a simple reflection.
    
    INPUT:
    - "type" -- the Cartan type of a Coxeter group W
    - "s" -- a number i, representing the simple reflection s=s_i
    - "w" -- a tuple, considered as an expression(not necessarily reduced) for
             an element of W.    

    OUTPUT:
    - the product c_s*c_w of the Kazhdan Lusztig basis elements c_s and c_w in
      the Hecke algebra of W, implemented as a dictionary, where the keys are
      tuples encoding the basis elements, and the keys are their coefficients
      in ZZ[v,v^(-1)].
    
    ALGORITHM:
    Recall that if sw<w in the Bruhat order, then 
        c_s * c_w = (v+v^(-1)) c_w;
    otherwise, 
        c_s * c_w = c_{sw} + \sum_{z:sz<z<w} \mu_{z,w} * c_z.

    .. NOTE: 
        We are intentionally implementing expressions in W as tuples instead of
        lists here, in accordance with the use of dictionaries to encode linear
        combinations, for the keys in a dictionary can be tuples but not list.
        Note, however, that the Sage functions such as W.bruhat_interval() all 
        use expressions implemented as lists.
    """
    
    W = CoxeterGroup(type, implementation = 'coxeter3')
    d = defaultdict()

    w_reduced = W(list(w)).reduced_word()   # a reduced word for w, as a list 
    
    if W(w_reduced).has_left_descent(s): 
        d[tuple(w_reduced)] = v + v**(-1)
    else:
        # Use the following three lines instead of d[(s,)+tuple(w_reduced)] = 1
        # to ensure that all tuples come from the normal form of elements
        t = (s,) + tuple(w_reduced)
        sw = W(list(t)).reduced_word()
        d[tuple(sw)] = 1 
        l = W.bruhat_interval([],w_reduced)
        # BIG assumpution here: we are assuming that elements in l are already
        # in normal forms (this appears to be true).
        for x in l:
            if s in x.left_descents() and x.mu_coefficient(W(w_reduced))!=0:
                d[tuple(x)] = x.mu_coefficient(W(w_reduced))
    return d



def remove_zero(d):
    dd = defaultdict(int)
    for k in d:
        if d[k] != 0:
            dd[k] += d[k]
    return dd


# computation of c_w * c_v. 
def sproduct_times_w(type,t,w):
    """ Multiply a product of c_{s}'s with c_w.

    INPUT:
    - 'type' -- the Cartan type of the Coxeter group W
    - 't' -- a tuple of simple reflections (s1,s2,...,sn)
    - 'w' -- a tuple, representing an element of W

    OUTPUT:
    - the product (c_s1 * c_s2 * ... * c_sn) * c_w in the Hecke algebra of W.
    
    """

    d = defaultdict()
    d[w] = 1
    for s in reversed(t):
        s_times_d = defaultdict(int)
        for w in d:
            dd = s_times_w(type,s,w)
            for term in dd:
                s_times_d[term] += dd[term] * d[w]
        d = s_times_d
    return d

def sproducts_times_w(type,d,w):
    result = defaultdict()
    for k in d:
        sproduct = tuple(itertools.chain.from_iterable(k))  # flatten d
        dd = sproduct_times_w(type,sproduct,w)
        for term in dd:
            if term in result:
                result[term] += dd[term] * d[k]
            else:
                result[term] = dd[term] * d[k]
    return remove_zero(result)



def break_elt_once(type,w):
    """
    INPUT:
    - 'w': an element of the Coxeter group, implemented as a tuple.

    OUTPUT:
    - a dictionary whose keys are tuples(products) of tuples(KL-basis
      elements), and whose values are the coefficients of the products in
      the expression for c_w.
      Note: all but one key in the resulting dictionary are just a plain
      KL-basis elements, considered as a "product of one element". 
    """

    d=defaultdict()
    s = w[0]
    w1= w[1:]
    
    if len(w) == 1:
        d[(w,)] = 1
    else:
        d = s_times_w(type,s,w1)
        del d[w]      
        d = {(k,):-d[k] for k in d}
        d[((s,),w1)] = 1

    return d

def done(t):
    """ Determine if a tuple represents a product of c_s's.

    INPUT:
    - t: a tuple of tuples, encoding the product of the KL-basis elements whose
      indices are its entries(which are tuples themselves).

    OUTPUT: 
    - a boolean value, describing if all the basis elements in the product are 
      already simple relfections
    """
    return len(t[-1]) == 1


def break_product(type,t):
    """ Break further a product c_s*...*c_w where w is not a simple reflection.
    INPUT:
    - 't': a tuple(product) of tuples(elements), of the form 
           t = ((s1,),...,(sn,),(w)),
           where w is not a simple reflection, so that t needs breaking.

    """
    result = defaultdict()
    head = t[:-1]
    d = break_elt_once(type,t[-1])
    for tail in d:
        result[head + tail] = d[tail]
    return result



def break_elt(type,w):
    result = defaultdict()
    to_be_checked = defaultdict()
    to_be_checked[(w,)] = 1
    while bool(to_be_checked): #equivalent to 'while to_be_checked is nonempty'
        new_to_be_checked = defaultdict()
        for t in to_be_checked:
            if done(t):
                if t in result:
                    result[t] += to_be_checked[t]
                else:
                    result[t] = to_be_checked[t]
            else:
                dd = break_product(type,t)
                for k in dd:
                    if k in new_to_be_checked:
                        new_to_be_checked[k] += dd[k] * to_be_checked[t]
                    else:
                        new_to_be_checked[k] = dd[k] * to_be_checked[t]
        to_be_checked = new_to_be_checked
    return result
                
def v_times_w(type,v,w):
    """ Compute c_v*c_w in the Hecke algebra.
    """

    d = break_elt(type,v)
    return sproducts_times_w(type,d,w)


def cbasis_into_tbasis(type,w):
    W3=CoxeterGroup(type,implementation='coxeter3')
    l=W3.bruhat_interval([],w)
    H=IwahoriHeckeAlgebra(type,v,-v**(-1))
    C=H.Cp();T=H.T()
    return sum([kl(type,y,w)*T[y] for y in l])
     

