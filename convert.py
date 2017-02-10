def convert_exp(rank,w):
    t=tuple(rank+Integer(1)-int(i) for i in str(w))
    return int(''.join(map(str,t)))

def convert_cell(rank,l):
    return [convert_exp(rank,x) for x in l]

def number_to_tuple(n):
    return tuple(int(z) for z in str(n))

