#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 29 14:27:48 2025

@author: Apol Medrano
"""

import makepoly as m
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    # Define the Runge Function 
    runge = lambda x: (1 + 25*(x**2))**(-1) 

    # Initialize x_0 and x_1
    x0 = -1
    x1 = 1

    # x-values to feed our functions
    x_test = np.linspace(x0, x1, 200)

    # Construct the values of Runge Function
    y_runge = runge(x_test)

    '''
    ++++++++++++++++++++
    Task 1: 
    Construct P4 and Plot against Runge Function
    ++++++++++++++++++++
    '''
    
    # x-values and y-values for degree 4 polynomial
    x_4 = np.linspace(x0, x1, 5)
    y_4 = runge(x_4)
    
    try:
        # Construct P4
        F = m.makepoly(x_4, y_4)
        an = np.diagonal(F)
        inter_4 = np.array([m.evalpoly(x, x_4, an) for x in x_test])

        # Plot P4 against Runge function
        
    # Print error message
    except RuntimeError as e:
        print(e)

    '''
    ++++++++++++++++++++
    Task 2: 
    Construct P8 and Plot against Runge Function
    ++++++++++++++++++++
    '''

    # x-values and y-values for degree 4 polynomial
    x_8 = np.linspace(x0, x1, 9)
    y_8 = runge(x_8)

    try:
        # Construct P8
        F = m.makepoly(x_8, y_8)
        an = np.diagonal(F)
        inter_8 = np.array([m.evalpoly(x, x_8, an) for x in x_test])

        # Plot P8 against Runge function
        
    # Print error message
    except RuntimeError as e:
        print(e)


    '''
    ++++++++++++++++++++
    Task 3: 
    Construct P16 and Plot against Runge Function
    ++++++++++++++++++++
    '''

    # x-values and y-values for degree 4 polynomial
    x_16 = np.linspace(x0, x1, 17 )
    y_16 = runge(x_16)

    try:
        # Construct P4
        F = m.makepoly(x_16, y_16)
        an = np.diagonal(F)
        inter_16 = np.array([m.evalpoly(x, x_16, an) for x in x_test])

        # Plot P16 against Runge function
        
    # Print error message
    except RuntimeError as e:
        print(e)


    '''
    ++++++++++++++++++++
    Task 4: 
    Plot the relative error between P16 and Runge Function
    ++++++++++++++++++++
    '''

    # Plot

