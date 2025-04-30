#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:14:27 2025

@author: Apol Medrano
"""

import makepoly as m
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    # Given x-values and y-values
    x_i = np.array([1, 2, 3, 4, 5])
    y_i = np.array([1, 1/2, 1/3, 1/4, 1/5])

    # Construct the function f(x) = 1/x
    f = lambda x: 1/x

    # interval 0.25 <= x <= 7
    x_test =  np.linspace(1/4, 7, 200)



    try:
        
        # Create the co-efficients
        F = m.makepoly(x_i, y_i) 
        an = np.diagonal(np.array(F))
        
        # Interpolate testing points 
        inter_val = np.array([m.evalpoly(x, x_i, an) for x in x_test])
        
        # Exact values of testing points
        exact_val = f(x_test)
        
        # Plot degree 2 and exact polynomial on the interval 0.25 <= x <= 7 w/ f(1/x)
        plt.figure()
        plt.plot(x_test, inter_val, 'orange', label='degree 4')
        plt.plot(x_test, exact_val, 'blue', label='f(x) = 1/x')
        plt.title("degree 4 vs exact values")
        plt.xlabel("x-values")
        plt.ylabel("y-values")
        #plt.savefig("", dpi = 200)
        plt.show()
        plt.close()


        # Put a scatter plot of the nodes on the previous graph
        plt.figure()
        plt.plot(x_test, inter_val, 'orange', label='degree 4')
        plt.plot(x_test, exact_val, 'blue', label='f(x) = 1/x')
        plt.scatter(x_i, y_i, color='black', marker='o', label="exact")
        plt.title("degree 4 vs exact values")
        plt.xlabel("x-values")
        plt.ylabel("y-values")
        #plt.savefig("", dpi = 200)
        plt.show()
        plt.close()

    # Print out error messages
    except RuntimeError as e:
        print(e)     
