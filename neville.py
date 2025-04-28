#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 10:58:10 2025

@author: Apol Medrano
"""

import numpy as np

'''
Neville: To evaluate the interpolating polynomial P given x-values and y-values
at the number x which approximates f(x)   =====================================
input x : the value of x, we are approximating for f(x)
      xn : the x-values of x_0, x_1, x_2,..., x_n
      fxn : the y-values of f(x_0), f(x_1), f(x_2),..,f(x_n)

output returns matrix Q that holds P_{n}(x) = Q_{n,n}
'''

def Neville(x, xn, fxn):

    # Sets N to the size of x-values, and N2 to the size of y-values 
    N = np.size(xn)
    N2 = np.size(fxn)

    try:
        # Raises number of x-values and y-values are of different sizes
        if N != N2:
            raise RuntimeError(f"Number of x and y values should agree.  There are {N} x-values and {N2} y-values.") 
    
            
        # Initializes matrix Q which will hold our P_{n,n}(x) at Q[n,n]
        Q = np.zeros([N,N])

        # Sets the first column of Q to the y-values
        # That is: Q[k, 0] = f(x_k) for k = 0 .. n
        for k in range(0, N2):
            Q[k, 0] = fxn[k] 
         
        # Computes the matrix Q
        for i in range(1, N):  # iterates through the columns i of Q  
            for j in range(1, i+1): # iterates through the rows j of Q

                # computes Q[i,j]
                numerator = (x - xn[i-j])*Q[i, j-1] - (x - xn[i])*Q[i-1, j-1]
                denomator = xn[i] - xn[i-j]
                Q[i,j] = numerator/denomator

        # returns Q after successful computation
        return Q

    # prints error message
    except RuntimeError as e:
        print(e)

'''
++++++++++++++++++++
Test 1: 
Exercise 4 from the notes from class on section 3.2
++++++++++++++++++++
'''

xn = np.array([-1,1,2,4])
yn = np.array([1,2,3,2])
x = 0
P = np.array([[1, 0, 0, 0 ],[2, 3/2, 0,0], [3, 1,4/3, 0], [2, 4, 0, 16/15]])
Q = Neville(x,xn,yn)
assert(np.any(np.abs(P-Q))<1e-12)

'''
++++++++++++++++++++
Test 2:
This example test problem was made up by Apol.
Approximate 4^2 given f(x) = x^2, x_0 = 2, x_1 = 3, and x_2 = 5
++++++++++++++++++++
'''

xn = np.array([2,3,5])
yn = np.array([4,9,25])
x = 4
P = np.array([[4, 0, 0 ],[9, 14, 0], [25, 17, 16]])
Q = Neville(x,xn,yn)
assert(np.any(np.abs(P-Q))<1e-12)

'''
++++++++++++++++++++
Test 3:
Test to raise error if sizes of xn and yn differ: 
++++++++++++++++++++
'''
xn = np.array([np.pi,np.e,np.sqrt(2)])
yn = np.array([3,3])
x = 3
Q = Neville(x,xn,yn)

'''
++++++++++++++++++++
HW 3.2: 3a Pg. 121
Use Neville's Method to approximate sqrt(3) with the following function and values:
f(x) = 3^x and the values x0 = -2, x1 = -1, x2 = 0, x3 = 1, and x4 = 2
++++++++++++++++++++
'''
xn = np.array([-2, -1, 0, 1, 2])
yn = np.array([(3 ** -2), (3 ** -1), (3 ** 0), (3 ** 1), (3 ** 2)])
x = 0.5
Q = Neville(x, xn, yn)
print(f"Q = {Q}")
