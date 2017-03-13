""" Tuple Operations """

def remove_first(t,y):
    """ Remove the first occurrence of t from a tuple y.
    
    EXAMPLE:
        sage: remove_first(3,(1,2,3,4,3,2)) 
        sage: (1,2,4,3,2)
    """

    a = y.index(t)
    return y[:a]+y[a+1:]

def neighbor_indices(s,y):
    """ Count the elements not commuting with s before the first s in y. 
    
    EXAMPLE:
        sage: count_neighbors(3,(2,1,4,2,3))
        sage: 3
    """

    i = 0
    l = []
    while len(l)<2 and i<y.index(s):
            if y[i] == s+1 or y[i] == s-1:
                l = l+[i]
            else:
                l = l
            i = i+1
    return l

def onetwo_on_left(y):
    i = 1
    if 1 in y and 2 in y:
        s = y[min(y.index(1),y.index(2))]
    elif 1 in y:
        s = 1
    else: 
        s = 2
    parabolic = tuple()
    minimal = y
    while i < 5 and s in minimal and neighbor_indices(s,minimal) == []:
        parabolic = parabolic + (s,)
        minimal = remove_first(s,minimal)    
        i = i + 1
        s = 3 - s
    return parabolic + minimal 


