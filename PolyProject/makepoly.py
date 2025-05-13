#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 15:19:25 2025

@author: Apol Medrano
"""

import numpy as np

'''
makepoly: computes a Matrix that contains the divided difference coefficents given x and y-values
=================================================================================================
input:  xn : the x-values of x_0, x_1, x_2,..., x_n
      fxn : the y-values of f(x_0), f(x_1), f(x_2),..,f(x_n)

output: returns matrix F that contains the divided difference coefficents
        along the diagonal (entries F[i,i])
'''

def makepoly(xn, fxn):

    # Sets N to the size of x-values, and N2 to the size of y-values 
    N = np.size(xn)
    N2 = np.size(fxn)

    # Raises number of x-values and y-values are of different sizes
    if N != N2:
        raise RuntimeError(f"Number of x and y values should agree.  There are {N} x-values and {N2} y-values.") 

    # Initializes matrix F
    F = np.zeros([N,N])

    # Sets the first column of F to the y-values, that is, Q[k, 0] = fxn[k] for k = 0, 1, ... n
    for k in range(N2):
        F[k, 0] = fxn[k] 
     
    # Computes the matrix F
    for i in range(1, N):  # iterates through the columns i of Q  
        for j in range(1, i+1): # iterates through the rows j of Q

            # computes F[i,j]
            F[i,j] = (F[i, j -1] - F[i -1, j -1])/(xn[i] - xn[i-j])

    # returns F after successful completion
    return F
 
'''
evalpoly: perform polynomial evaulation given nodes and coefficients 
====================================================================
input:  x : the value that we are approximating 
        xn : the x-values of x_0, x_1, x_2,..., x_n
        an : the coefficients a0, a1, ..., an

output: px: the output the polynomial
'''
def evalpoly(x, xn, an):

    # Sets N to the size of x-values, and N2 to the size of y-values 
    N = np.size(xn)
    N2 = np.size(an)

    # Raises number of x-values and y-values are of different sizes
    if N != N2:
        raise RuntimeError(f"Number of x and a values should agree.  There are {N} x-values and {N2} a-values.") 

    # Initializes the polynomail with a_n, that is, an[-1]  
    px = an[N-1]

    # Computes the nested form Pn(x)
    for k in range(N-2, -1, -1):
        xd = x - xn[k]
        px = px * xd + an[k]

    # returns the evaulated polynomial after successful completion
    return px

'''
Following are tests for makepoly() and evalpoly()
'''
if __name__ == "__main__":
    '''
    ++++++++++++++++++++
    Test 1:
    Example 1 from the section 3.3 notes
    ++++++++++++++++++++
    '''

    xn = np.array([1, 2, 3, 4]) # x-values
    yn = np.array([3,2, 3, 2])  # f(x-values)

    try:
        # Compute the divided difference Matrix F 
        F = makepoly(xn,yn) 

        # Done by hand computation of Matrix F 
        Fhand = np.array([[3, 0, 0, 0], [2, -1, 0, 0], [3, 1, 1, 0], [2, -1, -1, -2/3]])

        # Check if F is computed accurately
        assert(np.all(np.abs(Fhand- F))<1e-12)

    # prints error message 
    except RuntimeError as e:
        print(e)
    
    '''
    ++++++++++++++++++++
    Test 2:
    Example 1 from the section 3.3 notes
    ++++++++++++++++++++
    '''
    # x value to approximate
    x = 4
    xn = np.array([1, 2, 3, 4]) # x-values 
    yn = np.array([3,2, 3, 2])  # f(x-values)

    try:
        # Compute divided difference Matrix F
        F = makepoly(xn,yn)

        # Extract co-efficents for polynomial approximate
        A = np.diagonal(F) 

        # Evaluate polynomial to interpolate at x
        P = evalpoly(x, xn, A)

        # The interpolated value we computed by hand
        yHand = 2
        
        # Check if it successful computed reasonable values
        assert(np.all(np.abs(yHand - P)<1e-12))

    # prints error message 
    except RuntimeError as e:
        print(e)

    '''
    ++++++++++++++++++++
    Test 3: 
    Different sizes for x and y-values for makepoly
    ++++++++++++++++++++
    '''
    # x value to approximate
    x = 4
    xn = np.array([1, 2, 3]) # x-values 
    yn = np.array([3,2, 3, 2])  # f(x-values)

    try:
        # Compute divided difference Matrix F
        F = makepoly(xn,yn)

        # Extract co-efficents for polynomial approximate
        A = np.diagonal(F) 

        # Evaluate polynomial to interpolate at x
        P = evalpoly(x, xn, A)
        
        # Check if it successful computed reasonable values
        assert(np.all(np.abs(yn - P)<1e-12))

    # prints error message 
    except RuntimeError as e:
        msg = f"Number of x and y values should agree.  There are 3 x-values and 4 y-values."
        assert str(e) == msg

    '''
    ++++++++++++++++++++
    Test 4: 
    Different sizes for x and coefficent-values for evalpoly
    ++++++++++++++++++++
    '''
    # x value to approximate
    x = 4
    xn = np.array([1, 2, 3,4]) # x-values 
    yn = np.array([3,2, 3, 2])  # f(x-values)

    try:
        # Compute divided difference Matrix F
        F = makepoly(xn,yn)

        # Extract co-efficents for polynomial approximate
        A = np.array([1,2]) 

        # Evaluate polynomial to interpolate at x
        P = evalpoly(x, xn, A)
        
        # Check if it successful computed reasonable values
        assert(np.all(np.abs(yn - P)<1e-12))

    # prints error message 
    except RuntimeError as e:
        msg = f"Number of x and a values should agree.  There are 4 x-values and 2 a-values."
        assert str(e) == msg



    '''
    ++++++++++++++++++++
    Test 5:
    Exercise 2 from the Section 3.3 notes
    ++++++++++++++++++++
    '''
    # x value to approximate
    x = -1
    xn = np.array([0, 1, 2, -1])    # x-values
    yn = np.array([1, 0, -3, -6])   # f(x-values) 

    try:
        # Compute divided difference Matrix F
        F = makepoly(xn, yn)

        # Done by hand computation of Matrix F 
        Fhand = np.array([[1, 0, 0, 0], [0, -1, 0, 0], [-3, -3, -1, 0], [-6, 1, -2, 1]])

        # Check if F is computed accurately
        assert(np.all(np.abs(Fhand- F))<1e-12)

        # Extract co-efficents for polynomial approximate
        A = np.diagonal(F) 
        
        # Evaluate polynomial to interpolate at x
        P = evalpoly(x, xn, A)

        # The interpolated value we computed by hand
        yHand = -6

        # Check if it successful computed reasonable values
        assert(np.all(np.abs(yHand - P)<1e-12))

    # prints error message 
    except RuntimeError as e:
        print(e)

    '''
    ++++++++++++++++++++
    Test 6: Degree 3
    P3(x) = 2x^3 - 3x^2 + x - 5 
    ++++++++++++++++++++
    '''
    # define polynomial of degree shown above
    p3 = lambda x: 2*x**3 - 3*x**2 + x - 5 
    xn = np.linspace(-10, 10, 4)        # 4 x-values for degree 3
    yn = p3(xn)                         # 4 y-values for degree 3     
    x_test = np.linspace(-10, 10, 100)  # Points to test

    try:
        # Create the co-efficents
        F = makepoly(xn, yn)
        an = np.diagonal(np.array(F))
        
        # Interpolate testing points 
        y_inter = evalpoly(x_test, xn, an)
        
        # Compute exact testing points values 
        y_exact = p3(x_test)
        
        # Check if F is computed accurately
        assert(np.all(np.abs(y_inter - y_exact))<1e-8)
        
    # prints error message 
    except RuntimeError as e:
        print(e)

    '''
    ++++++++++++++++++++
    Test 7: Degree 5
    P5(x) = x^5 - 2x^3 + 3x - 7
    ++++++++++++++++++++
    '''
    # define polynomial of degree shown above
    p5 = lambda x: x**5 - 2*x**3 + 3*x - 7
    xn = np.linspace(-10, 10, 6)        # 6 x-values for degree 5
    yn = p5(xn)                         # 6 y-values for degree 5
    x_test = np.linspace(-10, 10, 100)  # Points to test

    try:
        # Create the co-efficents
        F = makepoly(xn, yn)
        an = np.diagonal(F)
        
        # Interpolate testing points 
        y_inter = evalpoly(x_test, xn, an)
        
        # Compute exact testing points values 
        y_exact = p5(x_test)
        
        # Check if F is computed accurately
        assert(np.all(np.abs(y_inter - y_exact))<1e-8)
        
    # prints error message 
    except RuntimeError as e:
        print(e)

    '''
    ++++++++++++++++++++
    Test 8: Degree 7
    P7(x) = -(1/2)x^7 + x^6 - x^5 + x^4 - 2x^2 + 1 
    ++++++++++++++++++++
    '''
    # define polynomial of degree shown above
    p7 = lambda x : -0.5*x**7 + x**6 - x**5 + x**4 - 2*x**2 + 1 
    xn = np.linspace(-10, 10, 8)        # 8 x-values for degree 7
    yn = p7(xn)                         # 8 y-values for degree 7
    x_test = np.linspace(-10, 10, 100)  # Points to test

    try:
        # Create the co-efficents
        F = makepoly(xn, yn)
        an = np.diagonal(np.array(F))
        
        # Interpolate testing points 
        y_inter = np.array([evalpoly(x, xn, an) for x in x_test])
        
        # Compute exact testing points values 
        y_exact = p7(x_test)
        
        # Check if F is computed accurately
        assert(np.all(np.abs(y_inter - y_exact))<1e-8)
        
    # prints error message 
    except RuntimeError as e:
        print(e)
