def row(n,i):
    return (2,) * (i-2) + (3,1,3) + (2,) * (n-i-1)

def fan(n):
    row1 = (1,) + (0,) * n
    row2 = (0,) + (1,3) + (2,) * (n-2)
    l = tuple((0,)+row(n,i) for i in range(2,n))
    rown = (0,) + (2,) * (n-2) + (3,1)
    return (row1, row2)  + l + (rown,)


