#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:47:56 2025

@author: Apol Medrano
"""

import makepoly as m
import matplotlib.pyplot as plt
import numpy as np

'''
hermite: Computes the divided difference coefficents for a Hermite polynomial
=================================================================================================
input: xn : the x-values of x_0, x_1, x_2,..., x_n
       fxn : the y-values of f(x_0), f(x_1), f(x_2),..,f(x_n)
       fpxn : the f' values of f'(x_0), f'(x_1), ..., f'(x_n)

output: returns matrix Q that contains the divided difference coefficents
        for the Hermite interploating polynomial along the diagonal (entries Q[i,i])
'''
def hermite(xn, fxn, fpxn):

    # Grab sizes of argument arrays
    N = np.size(xn)
    N2 = np.size(fxn)
    N3 = np.size(fpxn)

    # raise error of sizes don't match
    if N != N2 or N != N3 or N2 != N3:
        raise RuntimeError(f"Number of x, f(x), or f'(x) values should agree.  There are {N} x-values, {N2} f(x)-values and {N3} f'(x)-values")
    
    # Create z-array to help compute the Q-matrix
    z = np.zeros([2*N])
    
    # Create Q-matrix 
    Q = np.zeros([2*N, 2*N])

    # Compute the first 2 two columns of Q
    for i in range(N):
        
        # Grab values needed for computation
        z[2*i] = xn[i]
        z[(2*i) + 1] = xn[i]
        Q[2*i, 0] = fxn[i]
        Q[(2*i)+1, 0] = fxn[i]
        Q[(2*i)+1, 1] = fpxn[i]

        # if i > 0, then compute Q[2i, 1]
        if i != 0:
            Q[2*i, 1] = (Q[2*i, 0] - Q[(2*i)-1, 0])/(z[2*i] - z[(2*i) -1])

    # Compute the values for the rest of the Q matrix
    for i in range(2, 2*N):
        for j in range(2, i+1):

            # Compute Q[i,j]
            Q[i,j] = (Q[i, j-1] - Q[i-1, j-1])/(z[i] - z[i-j])
    
    # After successful computation return Q
    return Q


'''
Following are tests for hermite()
'''
if __name__ == "__main__":

    '''
    ++++++++++++++++++++
    Test 1: Checking if Q matrix matches from
    Exercise 1 from the section 3.4 notes
    ++++++++++++++++++++
    ''' 

    xn = np.array([-1, 0, 1])   # x-values` 
    fxn = np.array([2, 29/20, 1]) # f(x)-values
    fpxn = np.array([-31/60, -29/60, -8/15]) # f'(x)-values
    n = np.size(xn)  # size of x-values

    try:
        # Compute Q matrix
        Q = hermite(xn, fxn, fpxn) 

        # Expected values of the Q matrix done by hand
        Qhand = np.array([[2,0,0,0,0,0],[2,-31/60,0,0,0,0],
        [29/20, -11/20, -1/30,0,0,0], [29/20, -29/60, 1/15, 1/10,0,0], 
        [1, -9/20, 1/30, -1/60,-7/120, 0 ], [1,-8/15, -1/12, -7/60, -1/20, 1/240]])
        
        # Check if Q and Qhand have the same values
        assert(np.all(np.abs(Q - Qhand)<1e-12))

    # print error message
    except RuntimeError as e:
        print(e)

    '''
    ++++++++++++++++++++
    Test 2:
    number of nodes, function values and derivative values are not the same size
    ++++++++++++++++++++
    '''
    xn = np.array([-1, 1])   # x-values` 
    fxn = np.array([2, 29/20, 1]) # f(x)-values
    fpxn = np.array([-31/60, -29/60, -8/15]) # f'(x)-values

    try:
        # Compute Q matrix
        Q = hermite(xn, fxn, fpxn) 

        Qhand = np.array([[2,0,0,0,0,0],[2,-31/60,0,0,0,0],
        [29/20, -11/20, -1/30,0,0,0], [29/20, -29/60, 1/15, 1/10,0,0], 
        [1, -9/20, 1/30, -1/60,-7/120, 0 ], [1,-8/15, -1/12, -7/60, -1/20, 1/240]])
        
    # print error message
    except RuntimeError as e:
        msg = f"Number of x, f(x), or f'(x) values should agree.  There are 2 x-values, 3 f(x)-values and 3 f'(x)-values"
        assert str(e) == msg

    '''
    ++++++++++++++++++++
    Test 3: Check if interpolates the nodes correctly from
    Exercise 1 from the section 3.4 notes
    ++++++++++++++++++++
    ''' 

    xn = np.array([-1, 0, 1])   # x-values` 
    fxn = np.array([2, 29/20, 1]) # f(x)-values
    fpxn = np.array([-31/60, -29/60, -8/15]) # f'(x)-values
    n = np.size(xn)  # size of x-values

    try:
        # Compute Q matrix
        Q = hermite(xn, fxn, fpxn)

        # Grab zn values
        zn = np.zeros(2*n)
        for j in range(n):
            zn[2*j] = xn[j]
            zn[(2*j) + 1] = xn[j]

        # Compute hermite polynomial
        yH = m.evalpoly(xn, zn, np.diagonal(Q))

        # Check hermite polynomial interpolates nodes correctly
        assert(np.all(np.abs(yH - fxn)<1e-12))
        
    except RuntimeError as e:
        print(e)

    '''
    ++++++++++++++++++++
    Test 4: 
    Exercise 2 from the section 3.4 notes
    ++++++++++++++++++++
    '''
    
    xn = np.array([-1, 0, 1])   # x-values` 
    fxn = np.array([1, -4, 9]) # f(x)-values
    fpxn = np.array([-9.5, 6, 17.5]) # f'(x)-values
    n = np.size(xn)  # size of x-values

    try:
        # Compute Q matrix
        Q = hermite(xn, fxn, fpxn)

        # Grab zn values
        zn = np.zeros(2*n)
        for j in range(n):
            zn[2*j] = xn[j]
            zn[(2*j) + 1] = xn[j]

        # Compute hermite polynomial
        yH = m.evalpoly(xn, zn, np.diagonal(Q))

        # Expected values of the Q matrix done by hand
        Qhand = np.array([[-1, 0,0,0,0,0], [-1, -9.5, 0,0,0,0],
                          [-4, -3, -6.5,0,0,0], [-4, 6, 9, 2.5, 0, 0],
                         [9, 13, 7, -1, -1.75, 0], [9, 17.5, 4.5, -2.5, -0.75, 1/2]])

        # Check if Q and Qhand have the same values
        # assert(np.all(np.abs(Q - Qhand)<1e-12))

        # Check hermite polynomial interpolates nodes correctly
        assert(np.all(np.abs(yH - fxn)<1e-12))
        
    except RuntimeError as e:
        print(e)

    '''
    ++++++++++++++++++++
    Test 5: Degree 5
    P5(x) = x^5 - 2x^3 + 3x - 7
    ++++++++++++++++++++
    '''
    # define polynomial of degree 5 and its derivative shown above
    p5 = lambda x: x**5 - 2*x**3 + 3*x - 7
    dp5 = lambda x: 5*x**4 - 6*x**2 + 3 

    xn = np.linspace(-10, 10, 3)        # 3 x-values for degree 5
    fxn = p5(xn)                        # 3 f(x)-values for degree 5
    fpxn = dp5(xn)                      # 3 f'(x)-values for degree 5
    x_test = np.linspace(-10, 10, 100)  # Points to test
    n = np.size(xn)  # size of x-values


    try:
        # Create the co-efficents
        Q = hermite(xn, fxn, fpxn)
        an = np.diagonal(Q)

        # Grab zn values
        zn = np.zeros(2*n)
        for j in range(n):
            zn[2*j] = xn[j]
            zn[(2*j) + 1] = xn[j]

        # Interpolate testing points 
        yH = np.array(m.evalpoly(x_test, zn, an))
        
        # Compute exact testing points values 
        y_exact = p5(x_test)
        
        # Check if F is computed accurately
        assert(np.all(np.abs(yH - y_exact))<1e-8)
        
    # prints error message 
    except RuntimeError as e:
        print(e)

    '''
    ++++++++++++++++++++
    Test 6:     
    Plot interpolated Runge with 5 nodes vs exact Runge
    ++++++++++++++++++++
    '''

    # Define runge function and its derivative
    Runge = lambda x: (1 + 25*x**2)**(-1)
    Rungep = lambda x: -50*x*(1+ 25*x**2)**(-2)

    # Grab x, f(x), and f'(x) values 
    xn = np.linspace(-1, 1, 5)
    f_exact = Runge(xn)
    fp_exact = Rungep(xn)
    n = np.size(xn)  # size of x-values

    # Interval to graph
    x_test = np.linspace(-1, 1, 200)
    
    try:
        # Create the co-efficents
        Q = hermite(xn, f_exact, fp_exact)
        an = np.diagonal(Q)

        # Grab zn values
        zn = np.zeros(2*n)
        for j in range(n):
            zn[2*j] = xn[j]
            zn[(2*j) + 1] = xn[j]

        # Interpolate testing points 
        yH = np.array(m.evalpoly(x_test, zn, an))

        # Actual graphing points
        yR = Runge(x_test)

        # Plot
        plt.figure()
        plt.plot(x_test,yR, 'blue', label='Runge function')
        plt.plot(x_test, yH, 'r', label="Hermite poly")
        plt.legend()
        plt.title("Hermite with 5 nodes vs exact Runge")
        plt.xlabel("x-values")
        plt.ylabel("y-values")
        #plt.savefig(“RungeHermite.png”, dpi = 200)
        plt.show()
        plt.close()

    except RuntimeError as e:
        print(e) 
    
    '''
    ++++++++++++++++++++
    Test 6:     
    Plot interpolated Runge with 9 nodes vs exact Runge
    ++++++++++++++++++++
    '''
    # Define runge function and its derivative
    Runge = lambda x: (1 + 25*x**2)**(-1)
    Rungep = lambda x: -50*x*(1+ 25*x**2)**(-2)

    # Grab x, f(x), and f'(x) values 
    xn = np.linspace(-1, 1, 9)
    f_exact = Runge(xn)
    fp_exact = Rungep(xn)
    n = np.size(xn)  # size of x-values

    # Interval to graph
    x_test = np.linspace(-1, 1, 200)
    
    try:
        # Create the co-efficents
        Q = hermite(xn, f_exact, fp_exact)
        an = np.diagonal(Q)

        # Grab zn values
        zn = np.zeros(2*n)
        for j in range(n):
            zn[2*j] = xn[j]
            zn[(2*j) + 1] = xn[j]


        # Interpolate testing points 
        yH = np.array(m.evalpoly(x_test, zn, an))

        # Actual graphing points
        yR = Runge(x_test)

        # Plot
        plt.figure()
        plt.plot(x_test,yR, 'blue', label='Runge function')
        plt.plot(x_test, yH, 'r', label="Hermite poly")
        plt.legend()
        plt.title("Hermite with 9 nodes vs exact Runge")
        plt.xlabel("x-values")
        plt.ylabel("y-values")
        #plt.savefig(“RungeHermite.png”, dpi = 200)
        plt.show()
        plt.close()

    except RuntimeError as e:
        print(e)
