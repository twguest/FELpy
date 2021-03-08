#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.0"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""


import numpy as np    

def lu_decomp3(a):
    """
    c,d,e = lu_decomp3(a).
    LU decomposition of tridiagonal matrix a = [c\d\e]. On output
    {c},{d} and {e} are the diagonals of the decomposed matrix a.
    """
    n = np.diagonal(a).size
    assert(np.all(a.shape ==(n,n))) # check if square matrix

    d = np.copy(np.diagonal(a)) # without copy (assignment destination is read-only) error is raised 
    e = np.copy(np.diagonal(a, 1))
    c = np.copy(np.diagonal(a, -1)) 

    for k in range(1,n):
        lam = c[k-1]/d[k-1]
        d[k] = d[k] - lam*e[k-1]
        c[k-1] = lam
    return c,d,e

def lu_solve3(c,d,e,b):
    """
    x = lu_solve(c,d,e,b).
    Solves [c\d\e]{x} = {b}, where {c}, {d} and {e} are the
    vectors returned from lu_decomp3.
    """
    n = len(d)
    y = np.zeros_like(b)

    y[0] = b[0]
    for k in range(1,n): 
        y[k] = b[k] - c[k-1]*y[k-1]

    x = np.zeros_like(b)
    x[n-1] = y[n-1]/d[n-1] # there is no x[n] out of range
    for k in range(n-2,-1,-1):
        x[k] = (y[k] - e[k]*x[k+1])/d[k]
    return x

from scipy.sparse import diags
def create_tridiagonal(size = 4):
    diag = np.random.randn(size)*100
    diag_pos1 = np.random.randn(size-1)*10
    diag_neg1 = np.random.randn(size-1)*10

    a = diags([diag_neg1, diag, diag_pos1], offsets=[-1, 0, 1],shape=(size,size)).todense()
    return a

# =============================================================================
# a = create_tridiagonal(4)
# b = np.random.randn(4)*10
# 
# print('matrix a is\n = {} \n\n and vector b is \n {}'.format(a, b))
# 
# c, d, e = lu_decomp3(a)
# x = lu_solve3(c, d, e, b)
# 
# print("x from our function is {}".format(x))
# 
# print("check is answer correct ({})".format(np.allclose(np.dot(a, x), b)))
# 
# 
# =============================================================================
## Test Scipy
from scipy.linalg import solve_banded

def diagonal_form(a, upper = 1, lower= 1):
    """
    a is a numpy square matrix
    this function converts a square matrix to diagonal ordered form
    returned matrix in ab shape which can be used directly for scipy.linalg.solve_banded
    """
    n = a.shape[1]
    assert(np.all(a.shape ==(n,n)))

    ab = np.zeros((2*n-1, n))

    for i in range(n):
        ab[i,(n-1)-i:] = np.diagonal(a,(n-1)-i)

    for i in range(n-1): 
        ab[(2*n-2)-i,:i+1] = np.diagonal(a,i-(n-1))


    mid_row_inx = int(ab.shape[0]/2)
    upper_rows = [mid_row_inx - i for i in range(1, upper+1)]
    upper_rows.reverse()
    upper_rows.append(mid_row_inx)
    lower_rows = [mid_row_inx + i for i in range(1, lower+1)]
    keep_rows = upper_rows+lower_rows
    ab = ab[keep_rows,:]


    return ab


def make_diags(diags):
    # Make a linear array for the whole matrix
    n = len(diags[0])
    a = np.zeros(n * n, dtype=diags[0].dtype)
    # Assign each diagonal to the right stride
    step = n + 1
    for i, diag in enumerate(diags):
        a[i:(n - i) * n:step] = diag
    # Reshape
    return a.reshape(n, n)
    

# =============================================================================
# ab = diagonal_form(a, upper=1, lower=1) # for tridiagonal matrix upper and lower = 1
# 
# x_sp = solve_banded((1,1), ab, b)
# print("is our answer the same as scipy answer ({})".format(np.allclose(x, x_sp)))
# =============================================================================


if __name__ == '__main__':

    a = make_diags([np.arange(1,100)*np.pi, np.arange(1,99)*np.pi])
    print("shape a", a.shape)
    
    ab = diagonal_form(a)
    print("shape ab", ab.shape)
    
    b = np.random.rand(99,1)
    b[8] = 0
    print("shape b", b.shape)
    
    x = solve_banded((1,1), ab, b)
    print(x.shape)
    
    rb = np.matmul(a,x) ## resulting b
    
    from matplotlib import pyplot as plt
    plt.plot(b-rb)