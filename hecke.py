import itertools
from collections import defaultdict

r"""

This file contains code for fast computation of products of Kazhdan-Lusztig
basis elements in Hecke algebras.

-- Key Example

We illustrate our algorithm using a very simple example. 
(Please see the next section for more remarks on notations.)
In the Hecke algebra of type ['A',4] (i.e., that of the symmetric group 'S_5'), 
we have

\[
c_{121}c_{121} = (c_{1}c_{21} - c_1)c_{121} = (c_1c_2c_1 - c_1) c_121
\]

In this example, we expressed the left factor 'c_{121}' as a sum of products of



-- Notations and Conventions

We use the notations and conventions used by Lusztig in his book "Hecke
algebras with unequal parameters".  

- Coxeter groups are given by types(e.g., ['A',4]) and usually denoted by 'W'.
- Simple reflections are usually denoted by 's' or 't'.
- Elements of Coxeter groups are denoted by lower letters like the 'w' or 'y'.
- Elements of Coxeter groups are implemented as tuples, e.g., 's_1s_2=(1,2)'.
- For 'w\in W', 'l(w)' denotes the length of 'w'.
- The Bruhat order is denoted by '<'. 
- Standard basis elements of Hecke algebras are denoted by 'T_w'.
- Kazhdan-Lusztig basis elements of Hecke algebras are denoted by 'c_w'.
- 'kl' always stands for 'Kazhdan-Lusztig'.
- We use the normalization '(T_s-v)(T_s+v^{-1})=0' for the quadratic relation
  in the definition of the Hecke algebra. 
- Kazhdan-Lusztig polynomials are denoted by 'p_{y,w}'.
- We will use the following fact: 
     'c_s \times c_w = (v+v^{-1}) c_w' if 'sw<w';
     'c_s \times c_w = c_{sw} + \sum_{y: sy<y<w} \mu_{y,w} c_y' if 'sw>w'.
- We call the nonnegative integer coefficient '\mu_{y,w}' a mu-coefficient.
- Computation of mu-coefficients is at the center of our algorithm. We will
  compute them by using the package Coxeter3(wrapped in Sage), which seems to
  be the fastest way.  
- In the future we may also need to compute Kazhdan-Lusztig polynomials using
  Coxeter3. Coxeter3 uses a normalization for Hecke algebras that is different
  from ours, so a conversion would be necessary(see the function kl(type,y,w)).

"""

    
def v_times_w(type,v,w):
    r""" Compute c_v*c_w in the Hecke algebra.

    INPUT: 
 
    - 'type' -- the Cartan type of a Coxeter group 'W'
    - 'v,w' -- tuples, representing elements of 'W'

    OUTPUT:
    - the product 'c_v * c_w' in the Hecke algebra

    EXAMPLES:

        sage: v_times_w(['A',4],(2,1,3,2,1,4),(2,1))
        sage: {(2,1,3,2,1,4): (v + 1/v)^2}

        sage: v_times_w(['A',4],(1,2,1), (1,2,1))
        sage: {(1,2,1): (v + 1/v)^3 - v - 1/v} 
    """

    d = break_elt(type,v)
    return sproducts_times_w(type,d,w)



def s_times_w(type,s,w):
    r""" 
    Compute the product 'c_s*c_w', where 's' is a simple reflection.
    
    INPUT:
    - "type" -- the Cartan type of a Coxeter group 'W'
    - "s" -- a number 'i', representing the simple reflection 's=s_i'
    - "w" -- a tuple, considered as an expression(not necessarily reduced) for
             an element of 'W'.    

    OUTPUT:
    - the product 'c_s*c_w' of the kl basis elements in the Hecke algebra of W,
      implemented as a dictionary, where the keys are elements in 'W' indexing
      the kl basis elements, and the values are their coefficients in
      '\Z[v,v^(-1)]'.

    EXAMPLES:

        sage: s_times_w(['A',4],1,(1,2))
        sage: {(1,2): v + 1/v}
        
        sage: s_times_w(['A',4],1,(2,3,1,4))
        sage: {(1,2,1,3,4):1, (1,3,4): 1}
    
    ALGORITHM: 
    
    Recall that if 'sw<w', then 'c_s * c_w = (v+v^(-1)) c_w'; 
    otherwise, 'c_s * c_w = c_{sw} + \sum_{y:sy<y<w} \mu_{y,w} c_y'.

    .. NOTE: 
        - the elements returned in this function are in their canonical reduced
          form; in particular, the expression 'w' or 'sw' may not be one of the
          keys in the dictionary.

        - We are intentionally implementing expressions in W as tuples instead
          of lists here, in accordance with the use of dictionaries to encode
          linear combinations, for the keys in a dictionary can be tuples but
          not lists.  Note, however, that the Sage functions such as
          W.bruhat_interval() all use expressions implemented as lists, so we
          frequently have to convert tuples to lists or vice versa.  
    """
    
    W = CoxeterGroup(type, implementation = 'coxeter3') 
    d = defaultdict()
    
    # The canonical reduced word for w, produced by Sage as a list: 
    w_reduced = W(list(w)).reduced_word()  
    
    if W(w_reduced).has_left_descent(s): 
        d[tuple(w_reduced)] = v + v**(-1)
    else:
        # Use the following three lines instead of d[(s,)+tuple(w_reduced)] = 1
        # to ensure that all tuples come from the normal form of elements
        t = (s,) + w
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
    """ 
    Remove keys with value 0 from a dict/remove terms with coeff. 0 from a sum.
    """

    dd = defaultdict(int)
    for k in d:
        if d[k] != 0:
            dd[k] += d[k]
    return dd


def sproduct_times_w(type,t,w):
    r""" 
    Multiply a product of 'c_{s}'s with 'c_w'.

    INPUT:
    - 'type' -- the Cartan type of the Coxeter group W
    - 't' -- a tuple of simple reflections '(s1,s2,...,sn)'
    - 'w' -- a tuple, representing an element of W

    OUTPUT:
    - the product (c_s1 * c_s2 * ... * c_sn) * c_w in the Hecke algebra of W.
    
    EXAMPLES:

        sage: sproduct_times_w(['A',4],(1,),(1,2))
        sage: {(1,2): v + 1/v}

        sage: sproduct_times_w(['A',4],(1,1),(1,2))
        sage: {(1,2): (v + 1/v)^2}

        sage: sproduct_times_w(['A',4],(2,1),(2,3,1,4))
        sage: {(1,2,1,3,4): v + 1/v, (2,1,3,4): 1}
        # compare with
        sage: s_times_w(['A',4],1,(2,3,1,4))
        sage: {(1,2,1,3,4): 1, (1,3,4): 1}
        sage: s_times_w(['A',4],2,(1,2,1,3,4))
        sage: {(1,2,1,3,4): v + 1/v}
        sage: s_times_w(['A',4],2,(1,3,4))
        sage: {(2,1,3,4): 1}
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
    r"""
    Compute the product of a sum and 'c_w', where the sum is a linear
    combination of products of 'c_s's.

    INPUT: 
    - 'type' -- the Cartan type of the Coxeter group W
    - 'd' -- a dictionary representing a linear combinations of products of
            'c_s's, so the keys are like the input 't' in sproduct_times_w
    - 'w' -- a tuple, representing an element of 'W'

    EXAMPLE:

        sage: sproducts_times_w(['A',4],{(2,):1,(1,3):v},(2,3,1,4))
        sage: {(1,2,3,2,1,4): v, (2,1,3,4): v + 1/v, (1,3,4): (v+1/v)*v}
        
        # compare the above with the following

        sage: sproduct_times_w(['A',4],(2,),(2,3,1,4))
        sage: {(2,1,3,4): v + 1/v}
        sage: sproduct_times_w(['A',4],(1,3),(2,3,1,4))
        sage: {(1,2,3,2,1,4): 1, (1,3,4): v + 1/v}

    """

    result = defaultdict()
    for k in d:
        dd = sproduct_times_w(type,k,w)
        for term in dd:
            if term in result:
                result[term] += dd[term] * d[k]
            else:
                result[term] = dd[term] * d[k]
    return remove_zero(result)



def break_elt_once(type,w):
    r"""
    First step for breaking 'c_w' into linear comb. of products of 'c_s's.

    INPUT:
    - 'type' -- the Cartan type of a Coxeter group 'W'
    - 'w' -- a tuple, representing an element of 'W' 

    OUTPUT:
    - a dictionary whose keys are tuples(products) of tuples(kl basis
      elements), and whose values are the coefficients of the products in
      the expression for c_w. See the following note for more.


    .. NOTE:
      - Suppose 'w' is a canonical reduced expression starting with 's', and
        let 'w1' be the subword of 'w' following the 's'. Recall that

        'c_s * c_{w1} = c_w + \sum_{y: sy<y<w1} \mu_{y,w1} c_y',

        therefore

        'c_w = c_s * c_{w1} - \sum_{y: sy<y<w1} \mu_{y,w1} c_y'.

        The dictionary break_elt_once(type,w) returns is the following union:

        '{((s,),(w1)):1} \cup \cup_{y:sy<y<w1} {(y,):-\mu_{y,w1}}'

    EXAMPLES:

        sage: break_elt_once(['A',4],(1,))
        sage: {((1,),): 1}

        sage: break_elt_once(['A',4],(1,2))
        sage: {((1,),(2,)): 1}

        sage: break_elt_once(['A',4],(1,2,1))
        sage: {((1,),): -1, ((1,),(2,1)): 1}

        sage: break_elt_once(['A',4],(1,2,3))
        sage: {((1,),(2,3)): 1}
    """

    d = defaultdict()
    s = w[0]
    w1= w[1:]
    
    W = CoxeterGroup(type,implementation='coxeter3')
    w = tuple(W(list(w)).reduced_word())


    if len(w) == 1:
        d[(w,)] = 1
    else:
        d = s_times_w(type,s,w1)
        del d[w]      
        d = {(k,):-d[k] for k in d}
        d[((s,),w1)] = 1

    return d
    return len(t[-1]) == 1


def break_further(type,t):
    r""" 
    Break further a product 'c_s*...*c_w' if 'w' is not a simple reflection.

    INPUT:
    - 't': a tuple(product) of tuples(elements), of the form 
           't = ((s1,),...,(sn,),(w))'
           where 'w' is not a simple reflection.

    OUTPUT:
    - a dictionary obtained from the dictionary break_elt_once(type,w) by
      prepending the tuple ((s1,),...,(sn)) to all the keys.

    EXAMPLES:

        sage: break_further(['A',4],((4,),(3,),(1,2,1)))
        sage: {((4,),(3,),(1,))): -1, ((4,),(3,),(1,),(2,1)): 1}
        # compare with
        sage: break_elt_once(['A',4],(1,2,1))
        sage: {((1,),): 1, ((1,),(2,1)): 1}

    """
    result = defaultdict()
    head = t[:-1]
    d = break_elt_once(type,t[-1])
    for tail in d:
        result[head + tail] = d[tail]
    return result



def done(t):
    r""" 
    Determine if a tuple represents a product of c_s's.

    INPUT:
    - t: a tuple of tuples, encoding the product of the KL-basis elements whose
      indices are its entries(which are tuples themselves).

    OUTPUT: 
    - a boolean value, describing if all the basis elements in the product are 
      already simple relfections
    
    EXAMPLES:
    
        sage: done(((1,),))
        sage: True

        sage: done(((1,),(2,)))
        sage: True

        sage: done(((1,),(2,1)))
        sage: False

    """

    return len(t[-1]) == 1


def break_elt(type,w):
    r"""
    Write 'c_w' as a linear combination of product of 'c_s's.
    
    INPUT:
    - 'type' -- the Cartan type of a Coxeter group 'W'
    - 'w' -- a tuple, representing an element of 'W' 
    
    OUTPUT:
    - a dictionary encoding how to write 'c_w' as a linear combination of
      products of 'c_s's 
      
    EXAMPLE:

        # For type ['A',4], 'c_{121} = c_1c_{21} - c_1 = c_1c_2c_1 - c_1', so
     
        sage: break_elt(['A',4], (1,2,1))
        sage: {(1,): -1, (1,2,1):1}

        sage: break_elt(['A',4], (2,1,3,2,1,4))
        sage: {(2,1,3,2,1,4): 1, (2,1,3,4): -1, (1,4): 1, (1,2,1,4): -1}
        
        # The above result shows 'c_{213214}' can be written as follows:
            'c_2c_1c_3c_2c_1c_4 - c_2c_1c_3c_4 + c_1c_4 + c_1c_2c_1c_4'
    """
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
                dd = break_further(type,t)
                for k in dd:
                    if k in new_to_be_checked:
                        new_to_be_checked[k] += dd[k] * to_be_checked[t]
                    else:
                        new_to_be_checked[k] = dd[k] * to_be_checked[t]
        to_be_checked = new_to_be_checked

    flatten_result = defaultdict()
    for k in result:
        sproduct = tuple(itertools.chain.from_iterable(k))  # flatten d
        flatten_result[sproduct] = result[k]
    return flatten_result
               









########################### Code for future use ###############################

var('v')

def kl(type,y,w):
    """ Compute a Kazhdan-Lusztig polynomial 'p_{y,w}'.

    INPUT:
    - "type" -- the Cartan type of a Coxeter group 'W'
    - "y", "w" -- tuples, considered as expressions(not necessarily reduced) for
                 elements in 'W'
   
    OUTPUT:
    - the Kazhdan-Lusztig polynomial 'p_{y,w}', returned in our convention. 


    .. TODO::
    - use Coxeter matrix instead of Cartan type to specify Coxeter groups.
    """
    
    # compute the kl-polynomial (in 'q') using Coxeter3
    W = CoxeterGroup(type,implementation='coxeter3')
    f = W.kazhdan_lusztig_polynomial(y,w) 
    # convert the Coxeter3 result using our normalization ('v' instead of 'q')
    g = f.substitute(q=v**2)*v**(len(y)-len(w))
    return g.expand()

def mu(type,y,w):
    """ Compute the mu_{y,w}, the mu-coefficient for a pair of elements y,w.
    
    .. TODO::
        Use the coxeter3 package to directly compute mu_{y,w} instead of 
        reading it off as a coefficient of p_{y,w}.
        The direct computation can be---and is---done in s_times_w, so this
        function is actually not used.

    .. NOTE::
        It is known that '\mu_{y,w}' equals the coefficient of 'v^{-1}' in the
        kl-polynomials 'p_{y,w}'
    """

    return kl(type,y,w).coefficient(v,-1)
