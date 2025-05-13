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
    inter_4 = []

    print(f"Degree 4 x_4 = {x_4} ")
    print(f"Degree 4 y_4 = {y_4} ")


    try:
        # Construct P4
        F = m.makepoly(x_4, y_4)
        an = np.diagonal(F)
        print(f"Degree 4 a_n = {an}")
        inter_4 = np.array([m.evalpoly(x, x_4, an) for x in x_test])

        # Plot P4 against Runge function
        plt.figure()
        plt.plot(x_test, inter_4, 'orange', label='degree 4')
        plt.plot(x_test, y_runge, linestyle=":", color='blue', label="runge function")
        plt.scatter(x_4, y_4, color='blue', marker="o", label='data nodes')
        plt.legend()
        plt.title("Degree 4 vs Runge function")
        plt.xlabel("x-values")
        plt.ylabel("y-values")
        # plt.savefig("degree 4", dpi = 200)
        plt.show()
        plt.close()
        
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
    inter_8 = []

    print(f"Degree 8 x_8 = {x_8} ")
    print(f"Degree 8 y_8 = {y_8} ")


    try:
        # Construct P8
        F = m.makepoly(x_8, y_8)
        an = np.diagonal(F)
        inter_8 = np.array([m.evalpoly(x, x_8, an) for x in x_test])

        # Plot P8 against Runge function
        plt.figure()
        plt.plot(x_test, inter_8, 'orange', label='degree 8')
        plt.plot(x_test, y_runge, linestyle=":", color='blue', label="runge function")
        plt.scatter(x_8, y_8, color='blue', marker="o", label='data nodes')
        plt.legend()
        plt.title("Degree 8 vs Runge function")
        plt.xlabel("x-values")
        plt.ylabel("y-values")
        # plt.savefig("degree 8", dpi = 200)
        plt.show()
        plt.close()

        
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
    inter_16 = []
    print(f"Degree 16 x_16 = {x_16} ")
    print(f"Degree 16 y_16 = {y_16} ")

    
    try:
        # Construct P4
        F = m.makepoly(x_16, y_16)
        an = np.diagonal(F)
        inter_16 = m.evalpoly(x_test, x_16, an)

        # Plot P16 against Runge function
        plt.figure()
        plt.plot(x_test, inter_16, 'orange', label='degree 16')
        plt.plot(x_test, y_runge, linestyle=":", color='blue', label="runge function")
        plt.scatter(x_16, y_16, color='blue', marker="o", label='data nodes')
        plt.legend()
        plt.title("Degree 16 vs Runge function")
        plt.xlabel("x-values")
        plt.ylabel("y-values")
        # plt.savefig("degree 16", dpi = 200)
        plt.show()
        plt.close()

        '''
        ++++++++++++++++++++
        Task 4: 
        Plot the relative error between P16 and Runge Function
        ++++++++++++++++++++
        '''

        plt.figure()
        plt.plot(x_test, np.abs(inter_16 - y_runge)/(y_runge), 'red', label='relative error')
        plt.title("Relative Error Between Degree 16 and Runge function")
        plt.xlabel("x-values")
        plt.ylabel("y-values")
        # plt.savefig("relative error", dpi = 200)
        plt.show()
        plt.close()
 
    # Print error message
    except RuntimeError as e:
        print(e)    

    '''
    Compund Graph of Degree 4, 8, and 16 polyinomials
    '''
    fig = plt.figure()
    ax_array= fig.subplots(2,2,squeeze=True)
    ax_array[0,0].scatter(x_4,y_4, label="nodes")
    ax_array[0,0].plot(x_test,y_runge, color='k', label=r"$(1+25x^2)^{-1}$")
    ax_array[0,0].plot(x_test, inter_4, "--", color = 'b', label="poly interp d=4")
    ax_array[0,0].set_title("Degree 4")
    ax_array[0,1].scatter(x_8,y_8, color="r", label="nodes")
    ax_array[0,1].plot(x_test,y_runge, color='k', label=r"$(1+25x^2)^{-1}$")
    ax_array[0,1].plot(x_test, inter_8, "--", color='r', label="poly interp d=8")
    ax_array[0,1].set_title("Degree 8")
    ax_array[1,0].scatter(x_16,y_16, color="g", label="nodes")
    ax_array[1,0].plot(x_test,y_runge, color='k', label=r"$(1+25x^2)^{-1}$")
    ax_array[1,0].plot(x_test, inter_16, "--", color='g', label="poly interp d=16")
    ax_array[1,0].set_title("Degree 16")
    ax_array[1,1].plot(x_test,np.abs(inter_4 - y_runge)/(y_runge), color='b', label=r"absolute error")
    ax_array[1,1].plot(x_test,np.abs(inter_8 - y_runge)/(y_runge), color='r', label=r"absolute error")
    ax_array[1,1].plot(x_test,np.abs(inter_16 - y_runge)/(y_runge), color='g', label=r"absolute error")
    ax_array[1,1].set_title("Relative Error")
    # fig.savefig("runge_multi.png", dpi=200)
    plt.show()
    plt.close()

    '''
    Compund Graph of Relative Error
    '''
