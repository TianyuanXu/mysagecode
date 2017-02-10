r"""
This file contains various sorts of conversions that are useful for dealing
with words in Coxeter groups.
"""

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

def convert_element(n,w):
    """ Change a word digit-wise by the 1-n, 2-(n-1), ... .
    
    EXAMPLES:
        sage: reverse_diagram(4,123)
        sage: 432

    .. NOTE:
      -- This operation is useful to switch between the conventions of Sage and
         and Coxeter3, since they label the vertices in some types of Coxeter
         groups in the opposite order. In particular, for type H, Coxeter3 puts
         the strong bond between 1 and 2, while Sage puts it between n and n-1.
         The element 124 in H4 in  Coxeter3  is therefore 431. 
      -- We will get cell data from Coxeter, convert them using this map to
         words for Sage, then compute minimal coset representatives using Sage.
    """
    t=tuple(n+Integer(1)-int(i) for i in str(w))
    return int(''.join(map(str,t)))

def convert_cell(n,cell):
    """ Apply convert_element to every element in a cell. 

    EXAMPLES:
        sage: convert_cell(4,[123,13,24])
        sage: [432, 42, 31]
    """

    return [convert_exp(n,w) for w in cell]


