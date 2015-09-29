R.<q,v>=LaurentPolynomialRing(QQ,2)

a4=['A',4]

S5_sage=CoxeterGroup(a4)

card=S5_sage.order()

s=S5_sage.simple_reflections()

w_word=[1,2,3,4,2,1]

w_sage=S5_sage.from_reduced_word([1,2,3,4,2,1])

S5_coxeter3=CoxeterGroup(a4,implementation='coxeter3')

def kl_coxeter3(y,w):
    f=S5_sage.kazhdan_lusztig_polynomial(y,w)
    g=f.substitute(q=v**2)*v**(len(y)-len(w))
    return g

H=IwahoriHeckeAlgebra(a4,v**2,-1)

T=H.T(); C=H.Cp()

# This isn't working so far. I believe there is a bug. The error message is "... has no attribute 'truncate'"
def kl_sage(y,w):
    KL=KazhdanLusztigPolynomial(S5_sage,q)
    f=KL.P(S5_sage.from_reduced_word(y),S5_sage.from_reduced_word(w))
    g=f.substitute(q=v**2)*v**(len(y)-len(w))
    return g

# This isn't working either, for some silly reason.
def C_to_T(w):
    l=S5_coxeter3.bruhat_interval([],w)
    exp=0
    for y in l:
        coeff=kl_coxeter3(y,w)
        exp=exp+coeff*T[y]
    return exp
