def dihedral_string(s, t, l):
    """Create the dihedral sting sts... of length l."""

    if is_even(l):
        return [s,t]*(l//2)
    elif is_odd(l):
        return [s,t]*(l//2)+[s]

def numbers_from(a,n,step=-2):
    """
    Return the list [a,a-2,a-4,...] of n numbers.
    """

    return [a + step *i for i in range(n)]

def dihedral_product(u,w,M):
    r""" 
    Compute t_u*t_w in J for dihedral sequences u,w.

    INPUT:
    -"u", "w" -- lists representing reduced expressions of subregular elements
                 The expressions must 'collide', meaning they need to involve
                 the same two simple refelctions and u needs to end with what w
                 starts with.

    -"m" -- the number m(s,t) where s,t are the simple reflections in u and w

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
    if u[-1]!=w[0]:
        return 0
    elif u[-2]!=w[1]:
        return [u[:-1]+w]
    else: 
        k=len(u)
        l=len(w)
        m=M[u[0]-1][u[1]-1]
        if m==0:
            S=numbers_from(k+l-1,min(k,l))
        elif k+l<m+1:
            S=numbers_from(k+l-1,min(k,l))
        elif k+l>= m+1:
            S=numbers_from(2*m-k-l-1,m-max(k,l))
        return [dihedral_string(u[0],u[1],l) for l in S]

def first_dihedral_segment(l):
    segment=[]
    for x in l:
        if x==l[0] or x==l[1]:
            segment=segment+[x]
        else: 
            break
    return segment

def remove_first_dihedral_segment(l):
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
    segment_list=[first_dihedral_segment(l)]
    remainder=remove_first_dihedral_segment(l)
    while len(remainder)!=0:
        segment_list.append(first_dihedral_segment(remainder))
        remainder=remove_first_dihedral_segment(remainder)
    return segment_list

def left_mult_by_dihedral(u,w,M):
    if u[-1]!=w[0]:
        return 0
    elif w==first_dihedral_segment(w):
        return dihedral_product(u,w,M)
    else: 
        head_of_w=first_dihedral_segment(w)
        body_of_w=remove_first_dihedral_segment(w)
        new_heads=dihedral_product(u,head_of_w,M)
        return [(x[:-1]+body_of_w) for x in new_heads]

def t_mult(u,w,M):
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


