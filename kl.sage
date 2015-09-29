"""
Background:
    There are two sets of conventions for Hecke algebra.
    Ours: 
    (H_s-v)(H_s+v^{-1})=0
    Theirs (the one used by Sage and Coxeter3):
    (T_s-q)(T_s+1)=0, 
    implying (T_s/q^{1/2}-q^{1/2})(T_s/q^{1/2}+q^{-1/2})=0

    So there is an algebra homomorphism \phi sending our H_s to their
    T_s/q^{1/2}, and sending our v to their q^{1/2}. Meanwhile, on each version
    of the Hecke algebras there is a "bar involution" of algebras, with ours 
    sending H_s to its inverse and v to v^{-1} and theirs sending T_s to its 
    inverse and sending q^{1/2} to q^{-1/2}. Thus, the bar involution simply 
    sends each term to its inverse in the Hecke algebra. Being an algebra
    homomorphism, \phi commutes with the bar involution. So the \phi-image of 
    our Kazhdan-Lusztig basis elements are their Kazhdan-Lusztig basis elements 
    (for Sage, the images are denoted by Cp, with C used for another basis),
    which are defined to be special bar-invariant elements in each case.

    Conversion between the two versions of Kazhdan-Lusztig(KL) polynomials:
    For us, the KL polynomials are characterized as follows:
    c_w=\sum_{y\le w} p_{y,w}(v)H_y;
    For Sage and Coxeter3, the KL polynomials are characterized as follows:
    c_w=q^{-l(w)/2}\sum_{y\le w} p_{y,w}(q)T_y.
    (Again, in Sage, this holds with Cp in place of C.)
    Using the fact that \phi(c_w)=c_w, it is then straightforward to figure out 
    how to convert the two versions of KL polynomials.
    
    From ours to theirs: 
    Mulitply p_{y,w}(v) by v^{l(w)-l(y)}, then replace each v with q^{1/2}. 
    This results in p_{y,w}(q).

    From theirs to ours:
    Multiply p_{y,w}(q) by q^{(l(y)-l(w))/2}, then replace each q with v^2.
"""

# variable setup

R.<q,v>=LaurentPolynomialRing(QQ,2)


# Computations using Sage's internal algorithms
# A Coxeter group can be defined using its Cartan type, or its Coxeter matrix.

m=[[1,4,3,4],[4,1,3,2],[3,3,1,2],[4,2,2,1]] 
W1=CoxeterGroup(m)


# Immediately declare s=W.simple_reflections() after defining W; 
# Indexing of s starts at 1. 
# Or declare something like [s1,s2,s3]=W.simple_reflections().



# Computations using Gap3-Cheive:
# attach pycox.py and define Coxeter groups using Cartan type or Cartan matrix.
# Computattion is faster than Sage but mostly only deals with finite groups.


# Computation using coxeter3, for faster computation of KL polynomials:

W2=CoxeterGroup(m,implementation='coxeter3')

def cycle_type(n):
    # Cartan type of the type A affine Coxeter group whose graph is an n-cycle.
    return ['A',n-1,1]

W_triangle=CoxeterGroup(cycle_type(3),implementation='coxeter3')

def cycle(l,n): 
    return l*n+[l[0]]

def kl(type,y,w):    
    # input: the Coxeter type of W, 
    # y,w-two elements given using reduced words, such as [1,2] for s1*s2   
    # output: the KL polynomial in our convention

    W=CoxeterGroup(type,implementation='coxeter3')
    f=W.kazhdan_lusztig_polynomial(y,w)
    g=f.substitute(q=v**2)*v**(len(y)-len(w))
    return g


# The following is not working so far.
def cbasis_into_tbasis(type,w):
    W3=CoxeterGroup(type,implementation='coxeter3')
    l=W3.bruhat_interval([],w)
    H=IwahoriHeckeAlgebra(type,v,-v**(-1))
    C=H.Cp();T=H.T()
    return sum([kl(type,y,w)*T[y] for y in l])
     
 
