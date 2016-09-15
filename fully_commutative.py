def heap(m, w):
    r""" Return the heap of a word in a Coxeter group.

    INPUT:
    
    - 'm' -- the Coxeter matrix of a Coxeter group $W$
    - 'w' -- an expression in $W$

    OUTPUT:
    - a poset, the heap of $w$.

    """

    l = len(w)
    f = lambda i,j: i < j and m[w[i]-1][w[j]-1] != 2
    return Poset((range(l),f))

def a(m,w):
    l = heap(m,w).antichains()
    return max([len(x) for x in l])

