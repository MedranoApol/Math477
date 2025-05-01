#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 14:01:37 2025

@author: Apol Medrano
"""

import makepoly as m
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":

    # Dataset given
    T_F = np.array([220, 230, 240, 250, 260, 270, 280, 290, 300])
    P_psi = np.array([17.188, 20.78, 24.97, 29.82, 35.42, 41.85, 49.18, 57.53, 66.98])

    # Values to test for the polynomial
    x_test = np.linspace(180, 350, 200)

    
    '''
    ++++++++++++++++++++
    Task a:
    Construct P2(t) with data nodes T0=220, T1=260, T2=300 with their respective values. 
    Plot interpolating polynomial and the data in the given table
    ++++++++++++++++++++
    '''
    # Nodes from T_F to create for our x-values
    T_a = np.array([220, 260, 300]) 

    # Values corresponding nodes for our y-values
    Psi_a =  np.array([17.188, 35.42, 66.98])

    try:
        # Create the co-efficients
        F = m.makepoly(T_a, Psi_a) 
        an = np.diagonal(np.array(F))
        
        # Interpolate testing points 
        y_inter_a = np.array([m.evalpoly(x, T_a, an) for x in x_test])

        # Plot the polynomal vs the data in the table
        plt.figure()
        plt.plot(x_test, y_inter_a, 'orange', label='degree 2')
        plt.scatter(T_F, P_psi, color='blue', marker='o', label="exact")
        plt.legend()
        plt.title("degree 2")
        plt.xlabel("Temperature (F)")
        plt.ylabel("Steam Pressure (psi)")
        #plt.savefig("", dpi = 200)
        plt.show()
        plt.close()

    # print error message
    except RuntimeError as e:
        print(e)

    '''
    ++++++++++++++++++++
    Task b:
    Construct P4(t) with nodes T0=220, T1=260, T2=300 with their respective values.
    Plot the interplotated P2(t) and P4(t) with the data in the given table
    ++++++++++++++++++++
    '''
    # Nodes from T_F to create for our x-values
    T_b = np.array([220, 240, 260, 280, 300]) 

    # Values corresponding nodes for our y-values
    Psi_b =  np.array([17.188, 24.97, 35.42, 49.18, 66.98])
    
    try:
        # Create the co-efficients
        F = m.makepoly(T_b, Psi_b) 
        an = np.diagonal(np.array(F))
        
        # Interpolate testing points 
        y_inter_b = np.array([m.evalpoly(x, T_b, an) for x in x_test])

        # Plot P2 vs P4 vs the data on the same graph
        plt.figure()
        plt.plot(x_test, y_inter_a, 'orange', label='degree 2')
        plt.plot(x_test, y_inter_b, 'green', label="degree 4")
        plt.scatter(T_F, P_psi, color='blue', marker='o', label="exact")
        plt.legend()
        plt.title("degree 2 vs degree 4")
        plt.xlabel("Temperature (F)")
        plt.ylabel("Steam Pressure (psi)")
        #plt.savefig("", dpi = 200)
        plt.show()
        plt.close()

    # print error message
    except RuntimeError as e:
        print(e)

    '''
    ++++++++++++++++++++
    Task c:
    Find and scatter plot the relative error between data and interpolated values.
    ++++++++++++++++++++
    '''
    
    try:
        # Create the co-efficients for degree 2
        F_2 = m.makepoly(T_a, Psi_a) 
        an_2 = np.diagonal(np.array(F_2))

        # Create the co-efficients for degree 4
        F_4 = m.makepoly(T_b, Psi_b) 
        an_4 = np.diagonal(np.array(F_4))
        
        # Interpolate the data set for degree 2 and 4
        y_2 = np.array(m.evalpoly(T_F, T_a, an_2))
        y_4 = np.array(m.evalpoly(T_F, T_b, an_4))

        # Plot the relative error for P2 vs P4 against exact values
        plt.figure()
        plt.scatter(T_F, np.abs(y_2 - P_psi)/(P_psi), color='red' , label="degree 2")
        plt.scatter(T_F, np.abs(y_4 - P_psi)/(P_psi), color='black', marker='o', label="degree 4")
        plt.legend()
        plt.title("Relative error comparision between degree 2 and degree 4")
        plt.xlabel("Temperature (F)")
        plt.ylabel("Difference from Steam Pressure (psi)")
        #plt.savefig("", dpi = 200)
        plt.show()
        plt.close()

    # print error message
    except RuntimeError as e:
        print(e)

    '''
    ++++++++++++++++++++
    Task d:
    Construct P8(t) with the given table values.
    Plot P2(t), P4(t), P8(t) and the data given.
    ++++++++++++++++++++
    ''' 
    
    try:
        # Create the co-efficients
        F = m.makepoly(T_F, P_psi)
        an = np.diagonal(np.array(F))
        
        # Interpolate testing points 
        y_inter_d = np.array([m.evalpoly(x, T_F, an) for x in x_test])
        
        # Plot P2 vs P4 vs the data on the same graph
        plt.figure()
        plt.plot(x_test, y_inter_a, 'orange', label='degree 2')
        plt.plot(x_test, y_inter_b, 'green', label="degree 4")
        plt.plot(x_test, y_inter_d, 'red', label="degree 8")
        plt.scatter(T_F, P_psi, color='blue', marker='o', label="exact")
        plt.legend()
        plt.title("degree 2 vs degree 4 vs degree 8")
        plt.xlabel("Temperature (F)")
        plt.ylabel("Steam Pressure (psi)")
        #plt.savefig("", dpi = 200)
        plt.show()
        plt.close()
    
    except RuntimeError as e:
        print(e)
    

