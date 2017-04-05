import itertools
from collections import defaultdict
from collections import OrderedDict

""" conversions between strings, lists and tuples """

def word_to_list(w):
    """ Return the digits in a word (number) w as a list of integers.

    EXAMPLES:
        sage: word_to_list(231)
        sage: [2,3,1]
    """
    return  [int(i) for i in str(w)]

def word_to_tuple(w):
    """ Like word_to_list, but to tuples.

    EXAMPLES:
        sage: word_to_tuple(231)_
        sage: (2,3,1)
    """
    return tuple(int(i) for i in str(w))

def list_or_tuple_to_word(l):
    """ Concatenate the digits in a list or tuple to form a word. 

    EXAMPLES:
        sage: list_or_tuple_to_word([2,3,1])
        sage: 231

        sage: list_or_tuple_to_word((2,3,1)) 
        sage: 231
    """
    return int(''.join(map(str,l)))


""" number of fully-commutative but not subregular elements """

def fc_cardinality(n):
    """ 
    Return the number of fully-commutative elements, including the
    identity. 
    
        sage: fc_cardinality(2)
        sage: 9

        sage: fc_cardinality(3)
        sage: 44
    """
    return binomial(2*n+2,n+1)-2**(n+2)+n+3

def subregular_cardinality(n):
    """ 
    Return the number of subregular elements, including the
    identity. 
    
        sage: subregular_cardinality(2)
        sage: 9

        sage: subregular_cardinality(3)
        sage: 19
    """
    return 2*n*n+1

def fcnsr(n):
    """ 
    Return the number of elements that are fully-commutative but not
    subregular.
    
        sage: fcnsr(2)
        sage: 0

        sage: fcnsr(3)
        sage: 25
    """
    return fc_cardinality(n)-subregular_cardinality(n)


""" heaps and a-values, for type H """

def heap(w):
    r""" 
    Return the heap of a word in a Coxeter group of type H.
        
    EXAMPLE:
        sage: heap(132)
        sage: Finite poset containing 3 elements 
    """
    w = word_to_list(w)
    l = len(w)
    f = lambda i,j: i < j and (w[i] == w[j] + 1 or w[i] == w[j] - 1)
    return Poset((range(l),f))

def a(w):
    r"""
    Compute the a-value of w.
    
    EXAMPLE:
        sage: a(132)
        sage: 2 

        sage: a(135)
        sage: 3
    """
    l = heap(w).antichains()
    return max([len(x) for x in l])


""" canonical word from heap """

def leveled_heap(w):
    """
    Return a dictionary encoding the drawing of the heap of w, with the levels
    of letters as keys.

    EXAMPLE:
        sage: leveled_heap(1231)
        sage: {1:[1], 2:[2], 3:[3,1]}

        sage: leveled_heap(132143)
        sage: {1:[1,3], 2:[2,4], 3:[1,3]}
    """

    w = word_to_list(w)
    l = len(w)
    levels = defaultdict(list)
    height = defaultdict(int)
    for i in range(l):
        height[i] = 1
        for j in range(i):
            if w[j] == w[i] + 1 or w[j] == w[i] - 1:
                height[i] = max(height[i], height[j]+1)
        # print i, w[i], height[i]
        levels[height[i]] += [w[i]]
    return OrderedDict(sorted(levels.items(),key=lambda t: t[0]))

def canonical_word(w):
    """
    Return a canonical word of w obtained by reading the heap of w from bottom
    to top, from left to right.

    EXAMPLE:
        sage: canonical_word(1231)
        sage: 1213

        sage: canonical_word(132143)
        sage: 132413
    """
    l = []
    d = leveled_heap(w)
    for i in d:
        l = l + sorted(d[i])
    return list_or_tuple_to_word(l)






""" Tuple Operations """

def remove_first(t,y):
    """ 
    Remove the first occurrence of t from a tuple y.
    
    EXAMPLE:
        sage: remove_first(3,(1,2,3,4,3,2)) 
        sage: (1,2,4,3,2)
    """

    a = y.index(t)
    return y[:a]+y[a+1:]

def neighbors_before(s,y):
    """ 
    Find the indices of the neighbors of s before its first appearance in y.
    Stop if two neighbors are found.

    Note: If the resulting list is empty, sy < y. If two neighbors are found,
    sy is still fully-commutative. 
 
    EXAMPLE:
        sage: neighbors_before(3,(2,1,4,2,3))
        sage: [0,2] 
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

def first_12(y):
    """ 
    Find the index of the first occurrence of 1 or 2 in y. Return -1 if
    there's no 1 or 2 in y.
    
    EXAMPLE:
        sage: first_12((1,3,4,2))
        sage: 0

        sage: first_12((3,2,4,1))
        sage: 1

        sage: first_12((3,4))
        sage: -1 
    """
    i = 0
    index = -1
    while i < len(y):
        if y[i] > 2:
            i = i + 1
        else: 
            index = i
            break
    return index

def first_neighbor(s,y):
    """ 
    Find the index of the leftmost neighbor of s in y. 

    Note: We will only call this function when we know there is such a
    neighbor.

    EXAMPLE:
        sage: first_neighbor(2,(5,3,1,4))
        sage: 1

        sage: first_neighbor(2,(5,4,1,2))
        sage: 2
    """
    i = 0
    while i < len(y):
        if y[i] == s-1 or y[i] == s+1:
            return i
        else: 
            i = i+1
def my_index(s,y):
    """ 
    Return the index of s in y if s is in y; otherwise return -1. 
    
    EXAMPLE:
        sage: my_index(1,(2,1,3))
        sage: 1

        sage: my_index(1,(2,3,4))
        sage: -1
    """
    try: 
        return y.index(s)
    except ValueError:
        return -1

""" canonical and justified words, via repeated 12-coset decomposition """

def before_12(y):
    """
    Return the part of y before the first appearance of 1 or 2.

    EXAMPLE: 
        sage: before_12((3,4,5))
        sage: (3,4,5)

        sage: before_12((3,4,1,2,5))
        sage: (3,4)

        sage: before_12((3,4,2,5))
        sage: (3,4)
    """

    if first_12(y) == -1:
        return y
    else:
        return y[:first_12(y)]

def after_12(y):
    """
    Return the part of y starting from the first 1 or 2. 

    EXAMPLE: 
        sage: after_12((3,4,5))
        sage: ()

        sage: after_12((3,4,1,2,5))
        sage: (1,2,5)

        sage: after_12((3,4,2,5))
        sage: (2,5)
    """
    if first_12(y) == -1:
        return tuple()
    else: 
        return y[first_12(y):]

def onetwo_on_left(y):
    """ 
    Return the coset decomposition y_1 * y_2 of y, where y_1 is in the
    parabolic subgroup generated by 1 and 2 and neither 1 or 2 reduces y_2.  
    
    EXAMPLE:
        sage: onetwo_on_left((1,2,3)) 
        sage: ((1,2), (3,))
        
        sage: onetwo_on_left((1,2,3,1))
        sage: ((1,2,1,),(3,))

        sage: onetwo_on_left((1,2,3,2,1))
        sage: ((1,2,1),(3,1))
    """

    parabolic = tuple()
    minimal = y

    if first_12(y) > -1:
        s = y[first_12(y)]
        i = 1
        while i < 5 and my_index(s,minimal) > -1 and neighbors_before(s, minimal) == []:
            s = minimal[first_12(minimal)]
            parabolic = parabolic + (s,)
            minimal = remove_first(s,minimal)    
            i = i + 1
            s = 3 - s
    return parabolic, minimal 

def my_append(l,x):
    """
    Append x to l if x is not the empty tuple; return l as is otherwise.

    EXAMPLE:
        sage: my_append(((1,2)),(3,))
        sage: ((1,2),(3,))

        sage: my_append((1,2),())
        sage: ((1,2))
    """
    if x == tuple():
        return l
    else: 
        return l + [x]

def left_justify(y):
    """ 
    Return the segments of the left justified word of y.

    EXAMPLE:
        sage: left_justify((1,2,3,1,4,2,1,5))
        sage: [(1,2,1),(3,4),(2,1),(5,)]

        sage: left_justify((3,2,1,3,2,5,1))
        sage: [(3),(2,1),(3,),(2,1),(5,)]
    """


    l = []
    remain = y
    while remain != tuple():
        x = onetwo_on_left(remain)[0]
        y = onetwo_on_left(remain)[1]
        ll = my_append(l,x)
        l = my_append(ll,before_12(y))
        remain = after_12(y) 
    return l

def right_justify(y):
    """ 
    Return the segments of the right justified word of y.
    
    EXAMPLE:
        sage: right_justify((1,2,3,1,4,2,1,5))
        sage: [(1,2),(3,4,5),(1,2,1)]
    
    """
    z = y[::-1]
    l = left_justify(z)
    ll = [i[::-1] for i in l]
    return ll[::-1]


""" Green's f-basis from right-justified word """

# variables for the six types of 12-chuncks in a right justified word. They
# represent, in order: c1c2-1, c1c2c1-c1, c2c1c2-c2, c1c2c1c2-2c1c2,
# c2c1c2-2c2, c2c1c2c1-2c2c1
var('A')
var('B')
var('C')
var('D')
var('E')
var('F')


def ffactors(y):
    """
    Return the factors of Green's f-basis vector for y from a right-justified
    word.
    
    EXAMPLE:
        sage: factors([1,2,1,3,4,2,1,5])
        sage: [A,3,4,5,B]

        sage: factors((1,2,3,1,4,2,1,5))
        sage: [A,3,4,5,B]

        sage: factors((3,1,2,3,1,4,5,2))
        sage: [3,A,3,4,5,1,2]

    """
    l = right_justify(y)
    factors = []
    for i in range(len(l)):
        seg = l[i]
        if seg[0] == 1 or seg[0] == 2:
            if seg == (1,) or seg == (2,) or seg == (2,1):
                factors = factors + [i for i in seg]
            elif seg == (1,2):
                if i < len(l) - 2 and l[i+2][0] == 1:
                    factors = factors + [A]
                else:
                    factors = factors + [1,2]
            elif seg == (1,2,1):
                factors = factors + [B]
            elif seg == (2,1,2):
                if i < len(l) - 2 and l[i+2][0] == 1:
                    factors = factors + [E]
                else:
                    factors = factors + [C]
            elif seg == (1,2,1,2):
                factors = factors + [D]
            elif seg == (2,1,2,1):
                factors = factors + [F]
        else:
            factors = factors + [j for j in seg]
    return factors

def factor_to_tuple(l):
    """
    Recover the right-justified word of the f-basis vector with the given
    factors.
    
    EXAMPLES:
        sage: factor_to_tuple([4,A,6,3,1])
        sage: [4,1,2,6,3,1]

        sage: factor_to_tuple([4,A,6,3,D])
        sage: [4,1,2,6,6,3,1,2,1,2]
    """
    t = tuple()
    for x in l:
        if x == A:
            t = t + (1,2)
        elif x == B:
            t = t + (1,2,1)
        elif x == C:
            t = t + (2,1,2)
        elif x == D:
            t = t + (1,2,1,2)
        elif x == E:
            t = t + (2,1,2)
        elif x == F:
            t = t + (2,1,2,1)
        else: 
            t = t + (x,)
    return t

def factor_to_poly(l):
    """
    Multiply the factors of a f-basis element to return the element as a
    polynomial.

    EXAMPLES: 
        sage: factor_to_poly([4,3,A])
        sage: {(4,3):-1, (4,3,1,2):1}

        sage: factor_to_poly([A,4,B])
        sage: {(1,2,4,1):-1, (4,1,2,1):-1, (4,1):1, (1,2,4,1,2,1):1}
    """
    d = defaultdict(int)
    d[tuple()] = 1
    for i in l:
        if i == A:
            for w in list(d):
                d[w+(1,2)] += d[w]
                d[w] = -d[w]
        elif i == B:
            for w in list(d):
                d[w+(1,2,1)] += d[w]
                d[w+(1,)] += -d[w]
                del d[w]
        elif i == C:
            for w in list(d):
                d[w+(2,1,2)] += d[w]
                d[w+(2,)] += -d[w]
                del d[w]
        elif i == D:
            for w in list(d):
                d[w+(1,2,1,2)] += d[w]
                d[w+(1,2)] += -2*d[w]
                del d[w]
        elif i == E:
            for w in list(d):
                d[w+(2,1,2)] += d[w]
                d[w+(2,)] += -2*d[w]
                del d[w]
        elif i == F:
            for w in list(d):
                d[w+(2,1,2,1)] += d[w]
                d[w+(2,1)] += -2*d[w]
                del d[w]
        else:
            for w in list(d):
                d[w+(i,)] += d[w]
                del d[w]
    return d


""" multiplication of c_s * c_w or c_w * c_s """

var('v')

def s_once(s,y):
    """ 
    Return the result of the first step in computing c_s * c_w.

    INPUTS:
    -- 's': a simple reflection.
    -- 'w': the reduced word of an f.c. element

    OUTPUT:
    -- a dictionary of pairs (poly, w), where more computation is potentially
    needed to finish the computation of c_s * c_y. Here poly is a polynomial in
    the c_s' to be multiplied onto c_w with the pair's value as the
    coefficient. If poly is empty, it means it's 1 and the pair multiplies to
    c_w.

    EXAMPLES:
        sage: s_once(1,(2,3))
        sage: {((), (1,2,3)):1}

        sage: s_once(1,(5,4,1,2))
        sage: {((), (5,4,1,2)): v+1/v}

        sage: s_once(3,(5,2,4,1,3,2))
        sage: {((), (3,5,2,4,1,3,2)): 1}

        sage: s_once(6,(4,1,2,1,7,3,6))
        sage: {((4,),(3,6,1)): -1, ((4,1,2), (3,6,1)):1}

        sage: s_once(1,(4,2,1,3))
        sage: {((4,),(1,2,1,3):1, ((4,),(1,3)):1}
    """
    d = defaultdict(int)
    left = tuple()
    if my_index(s,y) == -1:             # if s does not appear in y
        d[(s,)+y] += 1 
    elif neighbors_before(s,y) == []:   # this is equivalent to sy<y 
        d[y] += (v + v**(-1))     
    elif len(neighbors_before(s,y)) == 2: # so sy is f.c.
        d[(s,)+y] += 1
    elif s > 3 or (s == 3 and y[first_neighbor(3,y)] == 4):
        factors = ffactors(y)
        neighbor = first_neighbor(s,factors)
        left = tuple(factors[:neighbor])
        d[factor_to_tuple(factors[neighbor+1:])] += 1 
    else:           # s appears in y, has only 1 neighbor before it, and sy>y
        if s == 1: 
            first2 = y.index(2)      # locate the first 2 in y
            left = y[:first2]
            z = y[first2:]
            parabolic = onetwo_on_left(z)[0] # the 12-star move
            if parabolic == (2,):
                d[(1,)+z] += 1 
            elif parabolic == (2,1,2,1):
                d[z[1:]] += 1 
            else:
                d[(1,)+z] += 1 
                d[z[1:]] += 1 
        elif s == 2 and y[first_neighbor(2,y)] == 1:
            first1 = y.index(1) 
            left = y[:first1]
            z = y[first1:]
            parabolic = onetwo_on_left(z)[0] # the 21-star move
            if parabolic == (1,):
                d[(2,)+z] += 1 
            elif parabolic == (1,2,1,2):
                d[z[1:]] += 1 
            else:
                d[(2,)+z] += 1 
                d[z[1:]] += 1 
        elif s == 2 and y[first_neighbor(2,y)] == 3: # the 23-star move
            first3 = y.index(3) 
            left = y[:first3]
            d[y[first3+1:]] += 1
        elif s == 3 and y[first_neighbor(3,y)] == 2:
            parabolic = onetwo_on_left(y)[0] 
            if parabolic == (1,2,1):   # the most subtle case
                pass             # since c_s * c_y = 0
            else:
                first2 = y.index(2)
                left = y[:first2]
                d[y[first2+1:]] += 1 
    dd = factor_to_poly(left)
    e = defaultdict(int)
    for k in dd:
        for j in d:
            e[(k,j)] += dd[k] * d[j] 
    return e

def s_times_w(s,w):
    """
    Return c_s * c_w as a linear combination of KL basis elements.

    EXAMPLES: 
        sage: s_times_w(1,(2,3))
        sage: {123:1}

        sage: s_times_w(1,(5,4,1,2))
        sage: {1254: v+1/v}

        sage: s_times_w(3,(5,4,2,1,3,2))
        sage: {3524132: 1}

        sage: s_once(6,(4,1,2,1,7,3,6))
        sage: {146231: 1} 
    """
    todo = defaultdict(int)
    todo[((s,),w)] = 1
    done = defaultdict(int)
    while todo != {}:
        for pair in list(todo):
            if pair[0] == tuple():
                done[pair[1]] += todo[pair]
                del todo[pair]
            else:
                monomial = pair[0]
                w = pair[1]
                c = todo[pair]
                del todo[pair]
                s = monomial[-1]
                d = s_once(s,w)
                for k in d:
                    todo[(monomial[:-1]+k[0],k[1])] += d[k] * c
    return clean_up(done)

def clean_up(d):
    """
    Delete the keys with value 0 from d and convert elements to canonical words. 

    EXAMPLE:
        sage: clean_up({(1,2,3):1, (4,5):0})
        sage: {123: 1}

        sage: clean_up({3,4,1}:1, (4,1):0)
        sage: {134: 1}
    """
    dd = defaultdict(int)
    for k in d:
        if d[k] != 0:
            dd[canonical_word(list_or_tuple_to_word(k))] += d[k]
    return dd

def inverse(w):
    """
    Return the reverse of a word w.
    
    EXAMPLE:
        sage: inverse(123)
        sage: 321
    """
    return list_or_tuple_to_word(reversed(word_to_tuple(w)))

def w_times_s(w,s):
    """
    Return c_w * c_s as a linear combination of KL basis elements.

    EXAMPLES: 
        sage: w_times_s((3,2),1)
        sage: {321:1}

        sage: w_times_s((2,1,4,5),1)
        sage: {2415: v+1/v}
    """
    y = w[::-1]
    d = s_times_w(s,y)
    dd = defaultdict(int)
    for k in d:
        dd[canonical_word(inverse(k))] += d[k]
    return dd



""" cells """

def left_descendents(x,n):
    """
    Return all elements lower than w in the left KL order in H_n.

    INPUT:
    -- 'w': the reduced word of an f.c. elements
    -- 'n': rank of the type-H Coxeter group

    Output: 
    -- A dictionary encoding a directed graph, whose vertices are the keys and
    the value of each key v is a list containing the vertices with an
    incoming edge coming from v. The vertices have all been converted to words
    from tuples for easier reading.
    
    NOTE: 
    The output graph is generated in the following way: start with w as
    the first vertex and multiply each vertex v by all simple reflections. If a
    new KL basis element appears in such a product, add the element to the
    vertex set and draw an edge from v to it. Continue until new elements arise
    this way. 
    
    EXAMPLE: 
        sage: left_descendents(1,2)
        sage: {1: [21], 21: [121,1], 2121: [121], 121: [2121, 21]} 
    """
    w = canonical_word(x)
    new_vertices = [w]
    current_vertices = [w]
    edges = defaultdict(list)
    while new_vertices != []:
        for z in new_vertices:
            # print z
            for s in range(1,n+1):
                d = s_times_w(s,word_to_tuple(z))
                for y in d:
                    if y != z and (y not in edges[z]):
                        edges[z] += [y]
                        # print s, y
                        if y not in current_vertices:
                            new_vertices = new_vertices + [y]
                            current_vertices = current_vertices + [y]
            new_vertices.remove(z)  
    d = defaultdict(list)
    for vert in edges:                # convert all tuples to words
        d[vert] = sorted(edges[vert])
    return OrderedDict(sorted(d.items(),key=lambda t: t[0]))

def descendents(x,n):
    """
    Return all elements lower than w in the left or right KL order in H_n.

    INPUT:
    -- 'w': the reduced word of an f.c. elements
    -- 'n': rank of the type-H Coxeter group

    Output: 
    -- A dictionary encoding a directed graph, whose vertices are the keys and
    the value of each key v is a list containing the vertices with an
    incoming edge coming from v.
    
    NOTE: 
    The output graph is generated in the following way: start with w as
    the first vertex and multiply each vertex v by all simple reflections. If a
    new KL basis element appears in such a product, add the element to the
    vertex set and draw an edge from v to it. Continue until new elements arise
    this way. 
    
    EXAMPLE: 
        sage: left_descendents(1,2)
        sage: {1: [21], 21: [121,1], 2121: [121], 121: [2121, 21]} 
    """
    w = canonical_word(x)
    new_vertices = [w]
    current_vertices = [w]
    edges = defaultdict(list)
    while new_vertices != []:
        for z in new_vertices:
            for s in range(1,n+1):
                d = s_times_w(s,word_to_tuple(z))
                for y in d:
                    if y != z and (y not in edges[z]):
                        edges[z] += [y]
                        # print s, y
                        if y not in current_vertices:
                            new_vertices = new_vertices + [y]
                            current_vertices = current_vertices + [y]
                dd = w_times_s(word_to_tuple(z),s)
                for y in dd:
                    if y != z and (y not in edges[z]):
                        edges[z] += [y]
                        # print s, y
                        if y not in current_vertices:
                            new_vertices = new_vertices + [y]
                            current_vertices = current_vertices + [y]
            new_vertices.remove(z)  
    d = defaultdict(list)
    for vert in edges:                # convert all tuples to words
        d[vert] = sorted(edges[vert])
    return OrderedDict(sorted(d.items(),key=lambda t: t[0]))


def left_graph(w,n):
    """
    Return the directed graph encoded by left_descendents(w,n).

    Note: 
    If w has a-value 2, the vertices of the strongly connected component of the
    output graph containing w is exactly the left cell of w.

    EXAMPLE:
        sage: left_graph(13,4)
        sage: Digraph on 18 vertices

        sage: left_graph(13,5)
        sage: Digraph on 66 vertices

    """
    return DiGraph(left_descendents(w,n))

def descendents_graph(w,n):
    """
    Return the directed graph encoded by descendents(w,n).

    Note: 
    If w has a-value 2, the vertices of the strongly connected component of the
    output graph containing w is exactly the 2-sided cell of w.

    EXAMPLE:
        sage: descendents_graph(13,4)
        sage: Digraph on 162 vertices

        sage: descendents_graph(13,5)
        sage: Digraph on 753 vertices
    """
    return DiGraph(descendents(w,n))


def left_cell(w,n):
    """
    Return the left cell in H_n of the element w (assuming a(w) = 2).

    EXAMPLE:
        sage: left_cell(13,4)
        sage: [13, 143, 213, 1213, 2143, 12143, 21213, 32143, 132143, 212143,
        321213, 2132143, 3212143, 4321213, 12132143, 212132143, 3212132143,
        43212132143]

        sage: left_cell(13,5)
        sage: [13, 143, 213, 1213, 1543, 2143, 12143, 21213, 21543, 32143,
        121543, 132143, 212143, 321213, 321543, 532143, 2121543, 2132143,
        3212143, 4321213, 4321543, 4532143, 12132143, 32121543, 53212143,
        54321213, 212132143, 432121543, 453212143, 3212132143, 43212132143,
        543212132143]

    """
    y = canonical_word(w)
    d = left_graph(y,n)
    return sorted(d.strongly_connected_component_containing_vertex(y))

def cell(w,n): 
    """
    Return the 2-sided cell in H_n of the element w (assuming a(w) = 2).

    EXAMPLE:
        sage: left_cell(13,4)
        sage: [13, 143, 213, 1213, 2143, 12143, 21213, 32143, 132143, 212143,
        321213, 2132143, 3212143, 4321213, 12132143, 212132143, 3212132143,
        43212132143]

        sage: left_cell(13,5)
        sage: [13, 143, 213, 1213, 1543, 2143, 12143, 21213, 21543, 32143,
        121543, 132143, 212143, 321213, 321543, 532143, 2121543, 2132143,
        3212143, 4321213, 4321543, 4532143, 12132143, 32121543, 53212143,
        54321213, 212132143, 432121543, 453212143, 3212132143, 43212132143,
        543212132143]
    """
    y = canonical_word(w)
    d = descendents_graph(y,n)
    return sorted(d.strongly_connected_component_containing_vertex(y))

def right_cell(w,n):
    """
    Return the right cell in H_n of the element w (assuming a(w) = 2).  
    """
    y = inverse(w)
    l = left_cell(y,n)
    ll = []
    for x in l:
        yy = inverse(x)
        tt = canonical_word(yy)
        ll = ll + [tt]
    return ll

def cell_inverse(C):
    """
    Return the inverse of a left cell, with elements in canonical form.
    """
    return [canonical_word(inverse(w)) for w in C]

    

def cell_intersection(w,n):
    """
    Return the intersection of the left cell containing w and its inverse in
    H_n, assuming a(w) = 2.
    
    EXAMPLE:

        sage: cell_intersection(13,4)
        sage:
    """
    A = left_cell(w,n)
    B = cell_inverse(A) 
    return Set(A).intersection(Set(B))


""" products c_x * c_y, where x, y have a-value at most 2 """

def monomial_times_y(t,y):
    """
    Compute the product of a monomial in c_s and c_y.
    
    INPUT:
    -- 't': a tuple encoding a monomial in c_s
    -- 'y': the reduced word of an element of a-value at most 2

    EXAMPLE:
        sage: monomial_times_y((1,2),3)
        sage: {123: 1}

        sage: monomial_times_y((1,2),31)
        sage: {1213: 1, 13:1}
    """
    d = defaultdict(int)
    d[y] = 1 
    for s in reversed(t):
        s_times_d = defaultdict(int)
        for w in d:
            dd = s_times_w(s,word_to_tuple(w))
            for term in dd:
                s_times_d[term] += dd[term] * d[w]
        d = s_times_d
    return d



def x_times_y(x,y):
    """
    Compute c_x * c_y modulo linear combinations of elements c_z where a(z)>2.

    INPUT:
    -- 'x', 'y': reduced words of elements of a-value at most 2

    EXAMPLE:
        sage: x_times_y(121,1212)
        sage: {12: v + 1/v}
        
        sage: x_times_y(132143,132143)
        sage: {132413: (v+1/v)^2, 13: (v+1/v)^2}
    """
    poly = factor_to_poly(ffactors(word_to_tuple(x)))
    d = defaultdict(int)
    for m in poly:
        dd = monomial_times_y(m,y)
        for term in dd:
            d[term] += poly[m] * dd[term]
    for k in list(d):
        if d[k] == 0:
            del d[k]
    return d




""" products of 2 commuting letters """

def a2_pairs(n):
    """
    Return all products ij (i<j) of 2 commuting letters i and j.

    EXAMPLE:
        sage: a2_pairs(5)
        sage: [13,14,15,24,25,35]
    """
    d = []
    for i in range(1,n+1):
        d += [list_or_tuple_to_word((i,j)) for j in range(i+2,n+1)]
    return d

def a2_cells(n):
    """
    Partition the 2-sided cell of a-value 2 in H_n into left cells.

    NOTE:
    The cardinality of the 2-sided cell with a-value 2 is 25, 162 = 2 * 9^2,
    392 = 2 * 14^2, 800 = 2 * 20^2, 1458 = 2 * 27^2 for H3, H4, H5, H6, H7,
    respectively.
    """
    C = cell(13,n)
    l = len(C)
    lcells = left_cells_in(C,n)
    k = len(lcells)
    intersections = [Set(c).intersection(Set(cell_inverse(c))) for c in lcells]
    with open('mysagecode/a2cells.txt','a+') as f:
        f.write("************** H" + str(n) + " *****************\n\n")
        f.write("The 2-sided cell of a-value 2 in H" + str(n) + " has " +
                str(l) + " elements.\n\n")
        f.write("The cell consists of " + str(k) + " left cells. For each cell,\nwe list its elements, its intersection with its inverse,\nand its unique distinguished involution below. \n\n")
        for i in range(k):
            f.write("Cell " + str(i+1) + ": ")
            f.write("%s" % lcells[i])
            f.write(".\n")
            f.write("The intersection: ")
            f.write("%s" % intersections[i])
            f.write(".\n")
            f.write("The distinguished involution: ")
            f.write("%s" % dist_inv(list(intersections[i])))
            f.write(".\n\n")


def left_cells_in(C,n):
    """
    Return the list of left cells in a 2-sided cell C.

    EXAMPLE:
        sage: left_cells_in(cell(13,3),3)
        sage: [[13, 213, 1213, 21213, 321213],
 			   [132, 2132, 12132, 212132, 3212132],
 			   [1321, 21321, 121321, 2121321, 32121321],
 			   [13212, 213212, 1213212, 21213212, 321213212],
 			   [132123, 2132123, 12132123, 212132123, 3212132123]]
    """
    remain = C
    l = []
    while remain != []:
        w = remain[0]
        lcell = left_cell(w,n)
        remain = [i for i in remain if i not in lcell]
        l = l + [lcell]
    return l

""" distinguished involutions """

def dist_inv(l):
    """
    Find the distinguished involution in the intersection l of a left cell and
    its intersection.

    EXAMPLE:
        sage: dist_inv(13,4)
        sage: {13}
        
        sage: dist_inv(1213,4)
        sage: {13}

        sage: dist_inv(13,5)
        sage: {13}
    """ 
    S = Set(x_times_y(inverse(l[0]),l[0]).keys())
    i = 1
    while len(S) > 1:
        SS = Set(x_times_y(inverse(l[i]),l[i]).keys())
        S = S.intersection(SS)
        i = i+1
    return list(S)[0]


