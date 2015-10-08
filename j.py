def dihedral_string(s, t, l):
    """ Create the dihedral list sts... of length l. """

    if is_even(l):
        return [s,t] * (l//2)
    elif is_odd(l):
        return [s,t] * (l//2) + [s]

def numbers_from(a,n,step=-2):
    """ Return the list [a,a-2,a-4,...] of n numbers.  """

    return [a + step * i for i in range(n)]

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
    
    sage: dihedral_product([1,2,1],[2,3],M)
    sage: 0

    sage: dihedral_product([1,2,1],[1],M)
    sage: [[1,2,1]]

    sage: dihedral_product([1,2,1],[1,3],M)
    sage: [[1,2,1,3]]

    sage: dihedral_product([2,3,2],[2,3,2,3,2],M)
    sage: [[2,3,2],[2,3,2,3,2],[2,3,2,3,2,3,2]]

    sage: dihedral_product([1,3,1,3],[3,1,3,1,3],M)
    sage: [[1,3,1,3],[1,3]]
    
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
    if u[-1] != w[0]:
        return 0
    elif len(u) == 1:
        return [w]
    elif len(w) == 1:
        return [u]
    elif u[-2] != w[1]:
        return [u[:-1] + w]
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
        return [dihedral_string(u[0],u[1],l) for l in S]

def first_dihedral_segment(l):
    """ Return the dihedral segment at the beginning of the list l.
    
    INPUT:
    - any nonempty list l, allowed to have length 1

    OUTPUT:
    - the first dihedral segment in l, so l itself when it has length 1. 
      The output is never empty.


    EXAMPLES:
    
    sage: first_dihedral_segment([1,2,1,2,3,4])
    sage: [1,2,1,2]

    sage: first_dihedral_segment([1,2,1])
    sage: [1,2,1]

    sage: first_dihedral_segment([1])
    sage: [1]

    """
    
    segment=[]
    for x in l:
        if x==l[0] or x==l[1]:
            segment=segment+[x]
        else: 
            break
    return segment

def remove_first_dihedral_segment(l):
    """ Return the empty list if l is dihedral; otherwise chop off all but the
        last number of the first dihedral segment of l.

    INPUT:
    - any nonempty list l

    OUTPUT:
    - empty if l is itself dihedral; 
      otherwise starts with the end of the first dihedral segment 

    Examples:
    
    sage: remove_first_dihedral_segment([1])
    sage: []

    sage: remove_first_dihedral_segment([1,2,1,2])
    sage: []

    sage: remove_first_dihedral_segment([1,2,1,3,4,2])
    sage: [1,3,4,2]
    
    """
    
    count = 0
    for x in l:
        if x==l[0] or x==l[1]:
            count = count + 1  
        else:
            break
    if count == len(l):
        return []
    else:
        return l[count - 1:]

def dihedral_segments(l):
    """ Return the dihedral segments of l in a list.

    INPUT:
    - any nonempty list l

    OUTPUT:
    - a list of (nonempty) dihedral segments

    EXAMPLES:
    
    sage: dihedral_segments([1])
    sage: [[1]]

    sage: dihedral_segments([1,2,1,2])
    sage: [[1,2,1,2]]

    sage: dihedral_segments([1,2,1,3,2,3,4,3,4,1])
    sage: [[1,2,1],[1,3],[3,2,3],[3,4,3,4],[4,1]]

    """
    segment_list=[first_dihedral_segment(l)]
    remainder=remove_first_dihedral_segment(l)
    while len(remainder) != 0:
        segment_list.append(first_dihedral_segment(remainder))
        remainder=remove_first_dihedral_segment(remainder)
    return segment_list

def left_mult_by_dihedral(u,w,M):
    """ Return the product t_u*t_w where u is dihedral.

    INPUT:
    - "u" -- a dihedral list, allowed to be of length only 1
    - "w" -- any nonempty list

    OUTPUT:
    - a list(sum) of lists(elements)

    EXAMPLES:
    
    sage: M=[[1,4,7],[4,1,0],[7,0,1]]

    sage: left_mult_by_dihedral([1,2,1],[1,3,1],M)
    sage: [[1,2,1,3,1]]

    sage: left_mult_by_diheral([1,2,1],[1,3,1,2,1],M)
    sage: [[1,2,1,3,1,2,1]]

    sage: left_mult_by_dihedral([1,3,1,3],[3,1,3,1,3,2,1],M)
    sage: [[1,3,1,3,2,1],[1,3,2,1]]

    """
    if u[-1] != w[0]:
        return 0
    elif w == first_dihedral_segment(w):    # this includes the case len(w)=1
        return dihedral_product(u,w,M)
    else: 
        head_of_w=first_dihedral_segment(w)
        body_of_w=remove_first_dihedral_segment(w)  # guaranteed to be nonempty
        new_heads=dihedral_product(u,head_of_w,M)
        return [(x[:-1]+body_of_w) for x in new_heads]

def t_mult(u,w,M):
    """ Return t_u*t_w for any subregular u and w.
    
    INPUTS:
    - "u", "w" -- any two lists
    - "M"      -- the Coxeter matrix, providing the m-values when necessary

    OUTPUT:
    - a list of lists (t_u*t_w as a sum of subregular elements)

    .. ALGORITHM:
    
    Break u into its dihedral segments d1, d2, ..., dn.
    Multiply each list(element) in [t_w] (list of list (t_w)) with t_dn 
    on the left using left_mult_by_dihedral, getting a list of lists.
    Repeat: multiply the resulting list of lists by t_d(n-1), ..., t_d1.

    """

    if u[-1]!=w[0]:
        return 0
    else:
        segments=dihedral_segments(u)
        n=len(segments)
        result=[w]
        while n>0:
            result=[y for x in result for y in
                left_mult_by_dihedral(segments[-1],x,M)]
            segments=segments[:-1]
            n=n-1
        return result


