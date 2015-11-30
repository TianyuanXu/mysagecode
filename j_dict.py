import itertools
from collections import defaultdict


def dihedral_string(s, t, l):
    """ Create the dihedral element sts... of length l as a tuple. """

    if is_even(l):
        return (s,t) * (l//2)
    elif is_odd(l):
        return (s,t) * (l//2) + (s,)

def numbers_from(a,n,step=-2):
    """ Return the tuple (a,a-2,a-4,...) of n numbers.  """

    return tuple(a + step * i for i in xrange(n))

def dihedral_product(u,w,M):
    """ Compute t_u*t_w in J for dihedral sequences u,w.  
    
    INPUT:
    -"u", "w" -- lists representing reduced expressions of subregular dihedral
                 elements. Both can have length only 1.

    -"M"      -- the Coxeter matrix, from which m-values are extracted if
                 necessary. If m(s,t)=\infty, the corresponding entry is
                 set to 0 instead.

    OUTPUT:
    - a list(representing the sum) of lists(representing elements). 
      The elements always have the same head as u and same tail as w.

    EXAMPLES:
    sage: M=[[1,4,7],[4,1,0],[7,0,1]]
    
    sage: dihedral_product((1,2,1),(2,3),M)
    sage: 0

    sage: dihedral_product((1,2,1),(1,),M)
    sage: [(1,2,1)]

    sage: dihedral_product((1,2,1),(1,3),M)
    sage: [(1,2,1,3)]

    sage: dihedral_product((2,3,2),(2,3,2,3,2),M)
    sage: [(2,3,2),(2,3,2,3,2),(2,3,2,3,2,3,2)]

    sage: dihedral_product((1,3,1,3),(3,1,3,1,3),M)
    sage: [(1,3,1,3),(1,3)]
    
    .. NOTE:
    The following is known: suppose u=st... and u,w collide, then
    \[ 
    t_u*t_w=\sum_{l\in S}t^{(n)}
    \]
    where t^{(n)} stands for t_x where x is the element sts... of length n, and
    S is the set that depends on k=len(u), l=len(w) and m.
        If m=\infty, then 
            S contains the \min(k,l) integers k+l-1, k+l-3, ... of the same
            parity.  
        If m is finite, then
            if k+l<m+1, then
                S still contains the \min(k,l) integers k+l-1, k+l-3, ... of
                the same parity.  
            if m-k+1\le l <m, then 
                S contains the m-(\max(k,l)) integers 2m-k-l-1, 2m-k-l-3, ...
                of the same parity.
    """
    d = defaultdict(int)
    if u[-1] != w[0]:
        d = {}
    elif len(u) == 1:
        d[w] += 1
    elif len(w) == 1:
        d[u] += 1
    elif u[-2] != w[1]:
        d[u[:-1] + w] += 1
    else: 
        k=len(u)
        l=len(w)
# Note the -1's below: we are assuming that the simple refelctions are labeled
# starting from 1, instead of 0 like the rows/columns of M.
        m=M[u[0]-1][u[1]-1]
        if m==0:
            S=numbers_from(k+l-1,min(k,l))
        elif k+l<m+1:
            S=numbers_from(k+l-1,min(k,l))
        elif k+l>= m+1:
            S=numbers_from(2*m-k-l-1,m-max(k,l))
        for l in S:
            d[dihedral_string(u[0],u[1],l)] += 1
    return d

def first_dihedral_segment(t):
    """ Return the dihedral segment at the beginning of the tuple t.
    
    INPUT:
    - any nonempty list t, allowed to have length 1

    OUTPUT:
    - the first dihedral segment in t, so t itself when it has length 1. 
      The output is never empty.


    EXAMPLES:
    
    sage: first_dihedral_segment((1,2,1,2,3,4))
    sage: (1,2,1,2)

    sage: first_dihedral_segment((1,2,1))
    sage: (1,2,1)

    sage: first_dihedral_segment((1))
    sage: (1)

    """
    return tuple(itertools.takewhile(lambda x: x == t[0] or x == t[1],t))

def is_dihedral(t):
    return all(i == t[0] or i == t[1] for i in t)


def remove_first_dihedral_segment(t):
    """ Return the empty tuple if t is dihedral; otherwise chop off all but the
        last number of the first dihedral segment of t.

    INPUT:
    - any nonempty list t

    OUTPUT:
    - empty if t is itself dihedral; 
      otherwise starts with the end of the first dihedral segment 

    Examples:
    
    sage: remove_first_dihedral_segment((1))
    sage: ()

    sage: remove_first_dihedral_segment((1,2,1,2))
    sage: ()

    sage: remove_first_dihedral_segment((1,2,1,3,4,2))
    sage: (1,3,4,2)
    
    """
    
    if is_dihedral(t):
        return ()
    else:
        return (first_dihedral_segment(t)[-1],) + tuple(itertools.dropwhile(lambda x: x == t[0] or x == t[1], t))

def dihedral_segments(t):
    """ Return the dihedral segments of t in a list.
    
    INPUT:
    - any nonempty tuple t

    OUTPUT:
    - a list of (nonempty) dihedral segments

    EXAMPLES:

    sage: dihedral_segments((1,))
    sage: [(1,)]

    sage: dihedral_segments((1,2,1,2))
    sage: [(1,2,1,2)]

    sage: dihedral_segments((1,2,1,3,2,3,4,3,4,1))
    sage: [(1,2,1),(1,3),(3,2,3),(3,4,3,4),(4,1)]
    """

    segment_list = []
    remainder = t
    while len(remainder) != 0:
        segment_list.append(first_dihedral_segment(remainder))
        remainder = remove_first_dihedral_segment(remainder)
    return segment_list


def left_mult_by_dihedral(u,w,M):
    """ Return the product t_u*t_w where u is dihedral.

    INPUT:
    - "u" -- a dihedral tuple, allowed to be of length only 1
    - "w" -- any nonempty tuple

    OUTPUT:
    - a list(sum) of tuples(elements)

    EXAMPLES:
    
    sage: M=[[1,4,7],[4,1,0],[7,0,1]]

    sage: left_mult_by_dihedral((1,2,1),(1,3,1),M)
    sage: [(1,2,1,3,1)]

    sage: left_mult_by_dihedral((1,2,1),(1,3,1,2,1),M)
    sage: [(1,2,1,3,1,2,1)]

    sage: left_mult_by_dihedral((1,3,1,3),(3,1,3,1,3,2,1),M)
    sage: [(1,3,1,3,2,1),(1,3,2,1)]

    """
    if u[-1] != w[0]:
        return {}
    elif is_dihedral(w):    # this includes the case len(w)=1
        return dihedral_product(u,w,M)
    else: 
        head_of_w=first_dihedral_segment(w)
        body_of_w=remove_first_dihedral_segment(w)  # guaranteed to be nonempty
        new_heads=dihedral_product(u,head_of_w,M) 
        d = defaultdict(int)
        for key in new_heads:
            d[key[:-1] + body_of_w] += new_heads[key]
        return d

def t_basis_mult(u,w,M):
    """ Return t_u*t_w for any subregular u and w.
    
    INPUTS:
    - "u", "w" -- any two tuples
    - "M"      -- the Coxeter matrix, providing the m-values when necessary

    OUTPUT:
    - a list of tuples (t_u*t_w as a sum of subregular elements)

    .. ALGORITHM:
    
    Break u into its dihedral segments d1, d2, ..., dn.
    Multiply each list(element) in [t_w] (list of list (t_w)) with t_dn 
    on the left using left_mult_by_dihedral, getting a list of lists.
    Repeat: multiply the resulting list of lists by t_d(n-1), ..., t_d1.

    """

    if u[-1]!=w[0]:
        return {}
    else:
        segments=dihedral_segments(u)
        n=len(segments)
        d = defaultdict(int)
        d[w] = 1
        while n>0:
            new_sum = defaultdict(int)
            for old_summand in d:
                new_summands = left_mult_by_dihedral(u,old_summand,M)
                for new_summand in new_summands:
                    new_sum[new_summand] += new_summands[new_summand] * d[old_summand]
            d = new_sum
            segments=segments[:-1]
            n=n-1
        return d

def print_basis_product(u,w,M):
    d = t_basis_mult(u,w,M)
    print 'The product is the sum of the following term(s):'
    for key in d:
        print d[key], '*' , key


def arbitrary_product(d1,d2,M):
    d = defaultdict(int)
    for u in d1:
        for w in d2:
            dd = t_basis_mult(u,w,M)
            for summand in dd:
                d[summand] += dd[summand] * d1[u] * d2[w]
    return d

def print_arbitrary_product(d1,d2,M):
    d = arbitrary_product(d1,d2,M)
    print 'The product is the sum of the following term(s):'
    for key in d:
        print d[key], '*', key


