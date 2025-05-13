#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 13 2:05 PM

@author: Apol Medrano
"""

import numpy as np
import matplotlib.pyplot as plt

'''
NaturalCubicSpline: creates nx4 matrix where Row j contains coefficients a_j, b_j, c_j, and d_j for the spline S_j(x) defined between the nodes x_j and x_j+1 
=====================================================================
INPUT: x_0,x_1, ..., x_n as the nodes and y_0, y_1, ..., y_n as the y-values for the nodes

OUTPUT: Coeffs, an nx4 matrix of coefficients. 
'''
def NaturalCubicSpline(xn, yn):

    # Grab sizes of argument arrays
    N = np.size(xn)
    N2 = np.size(yn)

    # raise error if nodes and y-values do not match
    if N != N2:
        raise RuntimeError(f"Number of nodes and y-values do not agree.There are {N} nodes and {N2} y-values.")

    # Create array h of length N and set h[0] = x_1 - x_0
    h = np.zeros(N-1)
    h[0] = xn[1] - xn[0]

    # Create array alpha of lenght n + 1 and set alpha[0] = alpha[n] = 0
    alpha = np.zeros(N)
    alpha[0] = 0
    alpha[N-1] = 0
    
    # Create Coeffs, a nx4 matrix
    Coeffs = np.zeros([N-1, 4])

    # Create a (n+1)x(n+1) matrix M and set M[0,0] = M[N,N] = 1
    M = np.zeros([N, N])
    M[0,0] = 1
    M[N-1,N-1] = 1

    # Compute entries of the matrix M 
    for k in range(1, N-1):
        h[k] =  xn[k+1] - xn[k]
        alpha[k] = (3/h[k])*(yn[k+1] - yn[k]) - (3/h[k-1]) * (yn[k] - yn[k-1])
        M[k, k-1] = h[k-1]
        M[k, k] = 2*(h[k-1] + h[k])
        M[k, k+1] = h[k]
    
    # We solve the system of equations to grab our coefficients c_j for j = 0, ..., n
    c = np.linalg.solve(M, alpha)

    # Compute the coefficients
    for j in range(0, N-1):
        a = yn[j]
        b = ((yn[j+1] - yn[j])/(h[j])) - ((h[j] * (c[j+1] + 2 * c[j]))/3)
        d = (c[j+1] - c[j])/(3 * h[j])
        Coeffs[j:] = np.array([a,b, c[j], d])

    # upon successful completion, it will return our Coeffs matrix
    return Coeffs


'''
evalpoly: evaluate interpolating polynomial given the divided differences
=====================================================================
input x, input x-value(x), and F the matrix of divided differences

output returns fx the output polynomial value(s) corresponding to input x
'''
def evalpoly(x,xn, A):
    N = np.size(A)
    N2 = np.size(xn)
    if N != N2:
        raise RuntimeError(f"Number of nodes and number of coefficients does not agree  There are {N2} nodes  and {N} coefficients.")
    px = A[N-1]
    for k in range(N-2,-1,-1):
        xd = x - xn[k]
        px = A[k]+px*xd
    return px

def evalCubicSpline(x, xn, Coeffs):
    N = np.size(xn)
    K = np.size(x)
    y = np.zeros(K)
    for i in range(K):
        A = Coeffs[N-2]
        y[i] = evalpoly(x[i], xn[N-2]*np.ones(4),A)
        for j in range(N-2, 0, -1):
            if x[i]<xn[j]:
                A = Coeffs[j-1]
                y[i] = evalpoly(x[i], xn[j-1]*np.ones(4),A)
    return y

if __name__ == "__main__":
    '''
    Test Data set 1
    '''
    xn = np.array([-4, -2, 0, 2, 4, 7])
    yn = np.array([-128, -16, 0, -40, 16, 51])

    # Compute the Coefficients for the spline
    Coeffs = NaturalCubicSpline(xn, yn)
    x = np.linspace(-4.5, 7.5, 100)
    y = evalCubicSpline(x, xn, Coeffs)

    plt.figure()
    plt.scatter(xn, yn,  color='r', label='nodes' )

    # Graph our cublic spline interpolation
    plt.plot(x, y, color="blue", label="cubic spline")

    plt.legend()
    plt.title("Test Data 1: Natural Cubic Spline Interpolation")
    # plt.savefig("myfile.png", dpi =200)
    plt.show()
    plt.close()

    ''' 
    Duck Data
    '''
    xd = np.array([0.9, 1.3, 1.9, 2.1, 2.6, 3.0, 3.9, 4.4, 4.7, 5.0, 6.0, 7.0, 8.0, 9.2, 10.5, 11.3, 11.6, 12.0, 12.6, 13.0, 13.3])
    yd = np.array([1.3, 1.5, 1.85, 2.1, 2.6, 2.7, 2.4, 2.15, 2.05, 2.1, 2.25, 2.3, 2.25, 1.95, 1.4, 0.9, 0.7, 0.6, 0.5, 0.4, 0.25])
    x = np.linspace(0.9, 13.3, 100)

    # Compute the Coefficients for the spline
    Coeffs = NaturalCubicSpline(xd, yd)
    y = evalCubicSpline(x, xd, Coeffs)

    plt.figure()
    plt.scatter(xd, yd,  color='r', label='nodes' )

    # Graph our cublic spline interpolation
    plt.plot(x, y, color="blue", label="cubic spline")

    plt.legend()
    plt.title("Duck Data and Natural Cubic Spline")
    plt.axis('equal')
    plt.grid()
    # plt.savefig("myfile.png", dpi =200)
    plt.show()
    plt.close()


    
