import hecke
import cacti



############################ Dihedral Case ####################################

def dihedral_ts_times_cw(rank,s,w):
    """ 
    Compute $T_s*c_w$ in the dihedral group(with symbols 1,2) with given rank.
    """

    d = defaultdict()
    if len(w) == rank:
        d[dihedral_string(1,2,rank)] = v
    elif s == w[0]:
        d[w] = v
    else:
        d[(s,)+w]=1
        d[w]=-v**(-1)
        if len(w) > 1: 
            d[w[1:]]=1
    return d

def dihedral_tproduct_times_cw(rank,t,w):
    d = defaultdict()

    d[w] = 1
    for s in reversed(t): 
        ts_times_d = defaultdict(int)
        for x in d:
            dd = dihedral_ts_times_cw(rank,s,x)
            for term in dd:
                ts_times_d[term] += dd[term] * d[x]
        d = remove_zero(ts_times_d)  # see below for definition of remove_zero
    return d

def dihedral_tw0_times_cw(rank,w):
    w0=dihedral_string(1,2,rank)
    d = dihedral_tproduct_times_cw(rank,w0,w)
    return compress_key(d)

def dihedral_inv(rank,l,w):
    d = dihedral_tw0_times_cw(rank,w)
    for x in l:
        if x in d:
            return x
        else: 
            continue

def dihedral_print_inv(rank):
    print 'The involution sigma behaves as follows in the dihedral group I_2({}):'.format(rank)
    print ''
    i = 1
    d = dihedral_lcells(rank)
    for k in d:
        print 'Left cell #{}: {}'.format(i,d[k])
        for elt in d[k]:
            print '{} --> {}'.format(elt,dihedral_inv(rank,d[k],tuple(int(i)
                for i in str(elt)))) 
        print ''
        i = i+1



