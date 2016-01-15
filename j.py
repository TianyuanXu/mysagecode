import itertools
from collections import defaultdict



def dihedral_string(s,t,l):
    r""" 
    Return the alternating tuple $(s,t,s,...)$ of length $l$. 
    
    EXAMPLES: 

        sage: dihedral_string(1,2,1)
        sage: (1,)

        sage: dihedral_string(2,4,5)
        sage: (2,4,2,4,2)
    """

    if is_even(l):
        return (s,t) * (l//2)
    elif is_odd(l):
        return (s,t) * (l//2) + (s,)

def numbers_from(a,n,d=-2):
    r""" 
    Return the tuple $(a,a+d,a+2d,...,a+(n-1)d)$. 

    EXAMPLES:

        sage: numbers_from(3,1)
        sage: (3,)

        sage: numbers_from(4,3)
        sage: (4,2,0)

        sage: numbers_from(5,4,3)
        sage: (5,8,11,14)
    """

    return tuple(a + d*i for i in xrange(n))


def dihedral_product(u,w,M):
    r""" 
    Compute $t_u*t_w$ in $J$ for dihedral sequences $u,w$.  
    
    INPUT:
    
    - "M"     -- the Coxeter matrix for an ambient Coxeter system $(W,S)$. 
                 If $m(s,t)=\infty$ for $s,t$ in $S$, the corresponding entry
                 is set to 0 instead.

    - "u, w"  -- tuples representing reduced expressions of subregular dihedral
                 elements in $W$. Both can have length only 1.

    OUTPUT:

    - a dictionary where the keys are tuples representing elements of $W$, and
      whose values are the corresponding coefficients in the product $t_u*t_w$.
      The empty dictionary {} is returned if the product is zero.

    EXAMPLES:

        sage: M = [[1,4,7],[4,1,0],[7,0,1]]
    
        sage: dihedral_product((1,2,1),(2,3),M)
        sage: {}

        sage: dihedral_product((1,2,1),(1,),M)
        sage: {(1,2,1):1}

        sage: dihedral_product((1,2,1),(1,3),M)
        sage: {(1,2,1,3):1}

        sage: dihedral_product((2,3,2),(2,3,2,3,2),M)
        sage: {(2,3,2,3,2,3,2):1,(2,3,2,3,2):1,(2,3,2):1}

        sage: dihedral_product((1,3,1,3),(3,1,3,1,3),M)
        sage: {(1,3,1,3):1,(1,3):1}
    
    .. NOTE::

        Note that the keys in a dictionary are "unordered", that is, not
        necessarily appear in the order in which they are inserted. For
        example, the previous example may actually appear as follows:

        sage: dihedral_product((1,3,1,3),(3,1,3,1,3),M)
        sage: {(1,3):1,(1,3,1,3):1}

    
    ALGORITHM:

    The following fact is known: 
    
    Let $u=...sts, w=pqp...$ be two subregular dihedral elements and $t_u, t_w$
    be the corresponding basis elements in $J$, then:
    
        If $s\neq p$, then $t_u \times t_w=0$;
    
        If $s=p$ and $t\neq q$, then $t_u \times t_w=t_v$ where $v=...stsqs...$;

        If $s=p$ and $t=q$, then
    
            \[ 
            t_u*t_w=\sum_{l\in S} t^{(n)}
            \]
        where $t^{(n)}$ stands for $t_v$ where $v$ is the element $sts...$ of 
        length $n$, and $S$ is the set that depends on $k=len(u), l=len(w)$ and 
        $m=m(s,t)$. To be precise: 

            If $m=\infty$, then $S$ contains the $\min(k,l)$ integers $k+l-1, 
            k+l-3, ...$ of the same parity.  
        
            If $m$ is finite, then
                if $k+l<m+1$, $S$ still contains the $\min(k,l)$ integers
                    $k+l-1, k+l-3, ...$ of the same parity;  
                if $m-k+1\le l <m$, then $S$ contains the $m-(\max(k,l))$
                    integers $2m-k-l-1, 2m-k-l-3, ...$ of the same parity.

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
        if m == 0:
            S=numbers_from(k+l-1,min(k,l))
        elif k+l<m+1:
            S=numbers_from(k+l-1,min(k,l))
        elif k+l>= m+1:
            S=numbers_from(2*m-k-l-1,m-max(k,l))
        for l in S:
            d[dihedral_string(u[0],u[1],l)] += 1
    return d

def first_dihedral_segment(t):
    r""" 
    Return the dihedral segment at the beginning of the tuple $t$.
    
    INPUT:

    - "t" -- any nonempty tuple , allowed to have length 1

    OUTPUT:

    - the longest dihedral string at the beginning of $t$, so $t$
      itself when it has length 1.  The output is never empty. 

    EXAMPLES:
    
        sage: first_dihedral_segment((1,2,1,2,3,4))
        sage: (1,2,1,2)

        sage: first_dihedral_segment((1,2,1))
        sage: (1,2,1)

        sage: first_dihedral_segment((1,))
        sage: (1,)

    """

    return tuple(itertools.takewhile(lambda x: x==t[0] or x==t[1],t))

def is_dihedral(t):
    r"""
    Return if a tuple is dihedral.
    """

    return all(i==t[0] or i==t[1] for i in t)


def remove_first_dihedral_segment(t):
    r""" 
    Return the empty tuple if $t$ is dihedral; otherwise delete all but the
    last number of the first dihedral segment of $t$.

    INPUT:
   
    - "t" -- any nonempty tuple
    

    OUTPUT:
    
    - empty if $t$ is itself dihedral; otherwise, the part of $t$ starting with
      the end of the first dihedral segment 

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
    r""" 
    Return the dihedral segments of $t$ in a list.
    
    INPUT:
    
    - "t" -- any nonemtpy tuple

    OUTPUT:

    - the list of the dihedral segments of $t$

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
    r""" 
    Return the product $t_u*t_w$ where $u$ is dihedral.

    INPUT:

    - "M" -- the Coxeter matrix for an ambient Coxeter system $(W,S)$.
             If $m(s,t)=\infty$ for $s,t$ in $S$, the corresponding entry is
             set to 0 instead.

    - "u" -- a tuple representing a dihedral element in $W$, allowed to have
             length only 1
    - "w" -- any nonempty tuple representing an element in $W$

    OUTPUT:

    - a dictionary encoding the product $t_u*t_w$, with tuples representing
      elements in $W$ as keys and their coefficients in the sum as values.
      The empty dictionary {} is returned if the product is zero.
   
    EXAMPLES:
    
        sage: M=[[1,4,7],[4,1,0],[7,0,1]]

        sage: left_mult_by_dihedral((1,2,1),(1,3,1),M)
        sage: {(1,2,1,3,1):1}

        sage: left_mult_by_diheral((1,2,1),(1,3,1,2,1),M)
        sage: {(1,2,1,3,1,2,1):1}

        sage: left_mult_by_dihedral((1,3,1,3),(3,1,3,1,3,2,1),M)
        sage: {(1,3,1,3,2,1):1,(1,3,2,1):1}

    .. NOTE::
    
        Note that the keys in a dictionary are "unordered", that is, not
        necessarily appear in the order in which they are inserted. For
        example, the previous example may actually appear as follows:

        sage: left_mult_by_dihedral((1,3,1,3),(3,1,3,1,3,2,1),M)
        sage: {(1,3,2,1):1,(1,3,1,3,2,1):1}
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

def t_basis_product(u,w,M):
    r""" 
    Return $t_u*t_w$ for any subregular elements $u$ and $w$.
    

    INPUT:
    
    - "M" -- the Coxeter matrix for an ambient Coxeter system $(W,S)$
             If $m(s,t)=\infty$ for $s,t$ in $S$, the corresponding entry is
             set to 0 instead.

    - "u, w" --  any nonempty tuples representing elements in $W$

    OUTPUT:

    - a dictionary encoding the product $t_u*t_w$, with tuples representing
      elements in $W$ as keys and their coefficients in the sum as values.
      The empty dictionary {} is returned if the product is zero.
    
    ALGORITHM: 

    Break $u$ into its dihedral segments $d_n, ..., d_2, d_1$. Multiply $t_w$ 
    by $t_{d_1}$ on the left, then multiply the result on the left by
    $t_{d_2}$, then $t_{d_3}$, etc.


    EXAMPLES:

        sage: M=[[1,4,7],[4,1,0],[7,0,1]]

        sage: t_basis_product((1,),(1,3,1),M)
        sage: {(1,3,1):1}

        sage: t_basis_product((1,2,1),(1,3,1,2,1),M)
        sage: {(1,2,1,3,1,2,1):1}

        sage: t_basis_product((2,1,3,1,3),(3,1,3,1,3,2,1),M)
        sage: {(2,1,3,1,3,2,1):1,(2,1,3,2,1):1}

    """
    if u[-1] != w[0]:
        return {}
    else:
        segments = dihedral_segments(u)
        n=len(segments)
        d = defaultdict(int)
        d[w] = 1
        while n>0:
            x = segments[-1]
            new_d = defaultdict(int)
            for term in d:
                new_terms = left_mult_by_dihedral(x,term,M)
                for new_term in new_terms:
                    new_d[new_term] += new_terms[new_term] * d[term]
            d = new_d
            segments=segments[:-1]
            n=n-1
        return d

def print_t_basis_product(u,w,M):
    r""" 
    Print the summands of $t_u*t_w$ for any subregular elements $u$ and $w$.
    

    INPUT:
    
    - "M" -- the Coxeter matrix for an ambient Coxeter system $(W,S)$
             If $m(s,t)=\infty$ for $s,t$ in $S$, the corresponding entry is
             set to 0 instead.

    - "u, w" --  any nonempty tuples representing elements in $W$

    OUTPUT:

    - the summands of $t_u*t_w$, printed line by line

    EXAMPLES:

        sage: M=[[1,4,7],[4,1,0],[7,0,1]]

        sage: print_t_basis_product((1,),(1,3,1),M)
        sage: t_(1,) * t_(1,3,1) 
              equals the sum of the following term(s):
              1 * (1,3,1)


        sage: print_t_basis_product((1,2,1),(1,3,1,2,1),M)
        sage: t_(1,2,1)* t_(1,3,1,2,1) 
              equals the sum of the following term(s):
              1 * (1,2,1,3,1,2,1)

        sage: print_t_basis_product((2,1,3,1,3),(3,1,3,1,3,2,1),M)
        sage: t_(2,1,3,1,3) * t_(3,1,3,1,3,2,1)
              equals the sum of the following term(s):
              1 * (2,1,3,1,3,2,1)
              1 * (2,1,3,2,1)
    """

    d = t_basis_product(u,w,M)
    print 't_{} * t_{} \nequals the sum of the following term(s):'.format(u,w)
    for key in d:
        print d[key], '*' , key


def arbitrary_product(d1,d2,M):
    r"""
    Return the product of any two elements in the ring $J_C$.
    
    INPUT: 
    
    - "M"      -- the Coxeter matrix for an ambient Coxeter system $(W,S)$. 
                  If $m(s,t)=\infty$ for $s,t$ in $S$, the corresponding entry
                  is set to 0 instead.

    - "d1, d2" -- dictionaries representing elements in $J_C$, the subring of
                  the asymptotic Hecke algbera $J$ for $W$ corresponding to the 
                  subregular cell of $W$. 
                  The keys of the dictionaries are tuples representing the
                  $t$-basis elements and the values are the correponding
                  coefficients in the linear combination.
    
    OUTPUT:

    - the product in $J_C$ of the elements represented by "d1" and "d2". 
    
    EXAMPLES:
       
        sage: M=[[1,4,7],[4,1,0],[7,0,1]]

        sage: arbitrary_product({(1,2,1):1},{(1,3,1,2,1):1},M)
        sage: {(1,2,1,3,1,2,1):1}

        sage: arbitrary_product({(2,3):2, (1,3,1,3):4}, {(3,1,3,1,3,2):5}, M)
        sage: {(2,3,1,3,1,3,2): 10, (1,3,1,3,2): 20, (1,3,2): 20}


        sage: arbitrary_product({(2,3):2, (1,3,1,3,1}:3}, {(1,3,1,2):3, (3,):7
              },M)
        sage: {(2,3):14, (1,3,1,3,1,2):9, (1,3,1,2):9}

    """

    d = defaultdict(int)
    for u in d1:
        for w in d2:
            dd = t_basis_product(u,w,M)
            for summand in dd:
                d[summand] += dd[summand] * d1[u] * d2[w]
    return d

def print_arbitrary_product(d1,d2,M):
    r"""
    Print the summands in the product of any two elements in the ring $J_C$.
    
    INPUT: 
    
    - "M"      -- the Coxeter matrix for an ambient Coxeter system $(W,S)$. 
                  If $m(s,t)=\infty$ for $s,t$ in $S$, the corresponding entry
                  is set to 0 instead.

    - "d1, d2" -- dictionaries representing elements in $J_C$, the subring of
                  the asymptotic Hecke algbera $J$ for $W$ corresponding to the 
                  subregular cell of $W$. 
                  The keys of the dictionaries are tuples representing the
                  $t$-basis elements and the values are the correponding
                  coefficients in the linear combination.
    
    OUTPUT:

    - the summands of the product in $J_C$ of the elements represented by "d1"
      and "d2", printed line by line
    
    EXAMPLES:
       
        sage: M=[[1,4,7],[4,1,0],[7,0,1]]

        sage: arbitrary_product({(1,2,1):1},{(1,3,1,2,1):1},M)
        sage: The desired product is the sum of the following term(s):
              1 * (1,2,1,3,1,2,1)

        sage: arbitrary_product({(2,3):2, (1,3,1,3):4}, {(3,1,3,1,3,2):5}, M)
        sage: The desired product is the sum of the following term(s):
              10 * (2,3,1,3,1,3,2)
              20 * (1,3,1,3,2)
              20 * (1,3,2)


        sage: arbitrary_product({(2,3):2, (1,3,1,3,1}:3}, {(1,3,1,2):3, (3,):7
              },M)
        sage: The desired product is the sum of the following term(s):
              14 * (2,3) 
              9 * (1,3,1,3,1,2)
              9 * (1,3,1,2)

    """
    d = arbitrary_product(d1,d2,M)
    print 'The desired product is the sum of the following term(s):'
    for key in d:
        print d[key], '*', key


def path(x,y):
    m = matrix(3,[1,x,2,x,1,y,2,y,1])
    return m


