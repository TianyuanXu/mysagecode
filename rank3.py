import pprint
import itertools


def f(b):
    return 2 * b + 3

def all_length(n):
    l = list(itertools.product([0,1],repeat=n))
    ll = [int(''.join(map(str,map(f,i)))) for i in l]
    return ll 
        
def all_less(n):
    l = [all_length(i) for i in range(1,n)]
    return reduce(lambda x,y: x+y, l)

def table(rows,cols,b):
    ll = [[str(i)]+[path_prod(i,j) for j in cols] for i in rows]
    top_row = ['0'] + [str(i) for i in cols]
    lll = [top_row] + ll
    for i in lll:
        i = [j.rjust(b) for j in i]
        print(''.join(i))

# for i in list_of_lists:
    # i = [str(j).rjust(2n) for j in i]
    # print(''.join(i))

def triangle(m12,m13,m23):
    m = matrix(3,[1,m12,m13,m12,1,m23,m13,m23,1])
    return m

def path(m12,m23):
    return triangle(m12,2,m23)

### code for the Coxeter group with m12=4, m23=6, m13=2

def code_to_word(s):
    t = ((1,)+dihedral_string(2,3,int(l))+(1,) for l in tuple(str(s)))
    word = reduce(lambda x,y: x+y[1:], t)
    return word

def path_prod(s,t):
    word_dic = t_basis_product(code_to_word(s),code_to_word(t),path(4,6))
    l = [str(word_to_code(key))+'*'+str(word_dic[key]) for key in word_dic]
    return '+'.join(l)

def word_to_code(t):
    x = dihedral_segments(t)
    if x[-1] == (1,):
        return 0
    elif x[-1] == (1,2,1):
        return 1
    else:
        codes = [len(x[2*i+1]) for i in range(len(x)//2)]
        n = int(''.join(map(str,codes)))
        return n 

