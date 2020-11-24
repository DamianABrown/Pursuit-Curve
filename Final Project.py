# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:43:39 2020

@author: damia
"""
"""
'''
This is the final project for Mathematical Programming Course.

Given a function (or curve) create a pursuit curve (using Forward Euler method) that will chase the 
function. Plot the curve, chasing curve, and its distance from curve and the time 
needed to catch the curve.
'''
"""


import numpy as np
import math
import matplotlib.pyplot as plt

## Input Values
k_value = 1
N = 50000
P_0 = [0, 0]
time_range = [0, 24*np.pi]
a = lambda a: math.sqrt((-13*np.sin(a)+(9*6)/13*np.sin(13*a/6))**2 + (13*np.cos(a) - (9*6)/13*np.cos(13*a/6))**2)
T = lambda t: [(13*np.cos(t)-(9*np.cos(13*t/6))), (13*np.sin(t)-(9*np.sin(13*t/6)))]
### Lambda Function a
### Copy equation and paste is spot of Lambda a. Match with same d() for Function T
# d(1) 1
# d(2) 1
# d(3) math.sqrt((-np.sin(a))**2 + (np.cos(a))**2)
# d(4) math.sqrt((-np.sin(a))**2 + (np.cos(a))**2)
# d(5) math.sqrt((-13*np.sin(a)+(9*6)/13*np.sin(13*a/6))**2 + (13*np.cos(a) - (9*6)/13*np.cos(13*a/6))**2)
### Lambda Function T
### Copy equation and paste is spot of Lambda T. Match with same d() for Function a
# d(1) [0, t]
# d(2) [0, t]   
# d(3) [np.cos(t), np.sin(t)]   
# d(4) [np.cos(t), np.sin(t)]   
# d(5) [(13*np.cos(t)-(9*np.cos(13*t/6))), (13*np.sin(t)-(9*np.sin(13*t/6)))]
## Creating matrix vector for Pursuit Curve
P_vector = np.zeros((2, N));T_vector = np.zeros((2, N))
def PursuitCurve():
    # Putting initial values into P matrix
    P_vector[0, 0] = P_0[0];P_vector[1, 0] = P_0[1]
    t, h = np.linspace(time_range[0], time_range[1], N, retstep=True)
    # Creating new variable for initial time value
    tt = t[0]
    time_Vector = T(tt)
    T_vector[0, 0] = time_Vector[0]
    T_vector[1, 0] = time_Vector[1]
    a_a = a(tt)   
    for ii in range(0,N-1): 
        b = math.sqrt((T_vector[0,ii] - P_vector[0,ii])**2 + (T_vector[1, ii] - P_vector[1,ii])**2) # ||T' - P'||
        x_prime = ((k_value*a_a)/b)*(T_vector[0, ii] - P_vector[0, ii])
        y_prime = ((k_value*a_a)/b)*(T_vector[1, ii] - P_vector[1, ii])
        x_n0 = P_vector[0, ii]
        y_n0 = P_vector[1, ii]
        ## Forward Euler Equation
        x_n1 = x_n0 + h*x_prime
        y_n1 = y_n0 + h*y_prime
       ## Putting n+1 value into designated matrix
        P_vector[0, ii+1] = x_n1
        P_vector[1, ii+1] = y_n1
        T_n1 = T(t[ii+1])
        T_vector[0, ii+1] = T_n1[0]
        T_vector[1, ii+1] = T_n1[1]
Distance_vector = np.zeros((2, N))
def DistanceBetween():
    t, h = np.linspace(time_range[0], time_range[1], N, retstep=True)
    for ii in range(0, N):
        dist = math.sqrt((P_vector[0, ii] - T_vector[0, ii])**2 + (P_vector[1, ii] - T_vector[1, ii])**2)
        Distance_vector[0, ii] = t[ii]
        Distance_vector[1, ii] = dist       
T_distance = []
P_distance = []
k_crit = [k_value]
PP_vector = np.zeros((2, N));T_vector = np.zeros((2, N))
def ComputeCriticalK():
    for ii in range(0,N-1):
        T_now = math.sqrt((T_vector[0, ii] - T_vector[0, ii+1])**2 + (T_vector[1, ii] - T_vector[1, ii+1])**2)
        T_distance.append(T_now)
        P_now = math.sqrt((P_vector[0, ii] - P_vector[0, ii+1])**2 + (P_vector[1, ii] - P_vector[1, ii+1])**2)
        P_distance.append(P_now)
    # Putting initial values into P matrix
    PP_vector[0, 0] = P_0[0];PP_vector[1, 0] = P_0[1]
    t, h = np.linspace(time_range[0], time_range[1], N, retstep=True)
    # Creating new variable for initial time value
    tt = t[0]
    time_Vector = T(tt)
    T_vector[0, 0] = time_Vector[0]
    T_vector[1, 0] = time_Vector[1]
    a_a = a(tt)
    if math.sqrt((P_vector[0, -1] - T_vector[0, -1])**2 + (P_vector[1, -1] - T_vector[1, -1])**2) < 0.001:
        new_k = k_value*(math.sqrt((P_vector[0, -1] - T_vector[0, -1])**2 + (P_vector[1, -1] - T_vector[1, -1])**2)+sum(P_distance))/sum(P_distance)
        print(new_k)
    else:
        for ii in range(0,N-1): 
            b = math.sqrt((T_vector[0,ii] - PP_vector[0,ii])**2 + (T_vector[1, ii] - PP_vector[1,ii])**2)
            x_prime = (((k_value)*a_a)/b)*(T_vector[0, ii] - PP_vector[0, ii])
            y_prime = (((k_value)*a_a)/b)*(T_vector[1, ii] - PP_vector[1, ii])
            x_n0 = PP_vector[0, ii]
            y_n0 = PP_vector[1, ii]
            ## Forward Euler Equation
            x_n1 = x_n0 + h*x_prime
            y_n1 = y_n0 + h*y_prime
            ## Putting n+1 value into designated matrix
            PP_vector[0, ii+1] = x_n1
            PP_vector[1, ii+1] = y_n1
            T_n1 = T(t[ii+1])
            T_vector[0, ii+1] = T_n1[0]
            T_vector[1, ii+1] = T_n1[1] 
        counter = 0
        while math.sqrt((PP_vector[0, -1] - T_vector[0, -1])**2 + (PP_vector[1, -1] - T_vector[1, -1])**2) > 0.0001:
            PP_distance = []
            for ii in range(0,N-1):
                P_now = math.sqrt((PP_vector[0, ii] - PP_vector[0, ii+1])**2 + (PP_vector[1, ii] - PP_vector[1, ii+1])**2)
                PP_distance.append(P_now)
            new_k = k_crit[counter]*(math.sqrt((PP_vector[0, -1] - T_vector[0, -1])**2 + (PP_vector[1, -1] - T_vector[1, -1])**2)+sum(P_distance))/sum(P_distance)
            k_crit.append(new_k)
            counter += 1
            print('while---', new_k)
            for ii in range(0,N-1):
                b = math.sqrt((T_vector[0, ii] - PP_vector[0, ii])**2 + (T_vector[1, ii] - PP_vector[1, ii])**2)
                x_prime = (((new_k)*a_a)/b)*(T_vector[0, ii] - PP_vector[0, ii])
                y_prime = (((new_k)*a_a)/b)*(T_vector[1, ii] - PP_vector[1, ii])
                x_n0 = PP_vector[0, ii]
                y_n0 = PP_vector[1, ii]
                ## Forward Euler Equation
                x_n1 = x_n0 + h*x_prime
                y_n1 = y_n0 + h*y_prime
                ## Putting n+1 value into designated matrix
                PP_vector[0, ii+1] = x_n1
                PP_vector[1, ii+1] = y_n1

def Plot():
    fig = plt.figure()
    plt1 = fig.add_subplot(221)
    plt2 = fig.add_subplot(222)
    plt1.set_title('x = Position, y = Time')
    plt2.set_title('x = Time, y = Distance')
    plt1.plot(T_vector[0,], T_vector[1,], 'b', P_vector[0,], P_vector[1,], 'g')
    plt2.plot(Distance_vector[0,], Distance_vector[1,], 'r')
    plt1.legend(['Exact', 'Forward Euler'])
    plt2.legend(['Distance'])
def Plot2():
    fig = plt.figure()
    plt1 = fig.add_subplot(221)
    plt1.set_title('x = Position, y = Time')
    plt1.plot(T_vector[0,], T_vector[1,], 'b', PP_vector[0,], PP_vector[1,], 'g')
    plt1.legend(['Exact', 'Forward Euler'])
PursuitCurve()
DistanceBetween()
#ComputeCriticalK()
Plot()
Plot2()
'''
plt.plot(tFE, yFE, 'b', tHM, yHM, 'r', tAB, yAB, 'y', tFE, np.exp(-tFE), 'g')
plt.legend(['Forward Euler', 'Heun', 'Adams-Bashforth 2', 'Exact'])
plt.xlabel('Time (t)')
plt.ylabel('Solution (x)')
'''
