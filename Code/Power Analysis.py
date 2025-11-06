"""
EE 430 Power Analytical Methods of Power Systems - Fall 2025
Term Project - Newton-Raphson Algorithm
Joshua Consenz - 11/3/25

Creates a five bus system with six transmission lines, and implements the Newton-Raphson algorithm to solve the system
from an initial state to a steady state.
"""

import math
import cmath
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Functions import *
from bus import Bus
from t_line import T_line

# Base values and tolerance of the system
baseMVA = 100
V_Tolerance = 0.05

def to_polar(y : complex):
    """
    Takes a complex number, y, and converts it into a tuple with magnitude and phase in degrees
    :param y: input complex number in rectangular format
    :return: a tuple of the complex number as magnitude and phase in degrees
    """
    r, theta = cmath.polar(y)
    r = round(r,3)
    theta = round(math.degrees(theta),3)
    return r, theta

def build_ybus_rect(buses, t_lines):
    """
    Function for creating a Ybus matrix from a system if buses and transmission lines
    :param buses: 1xn array of every bus in the system using the Bus class
    :param t_lines: 1xk array of the transmission lines of the system connecting buses using the T_line class
    :return: an nxn Ybus matrix using the admittance property of each transmission line
    """
    y_bus = np.zeros((len(buses),len(buses)), complex) # create empty Ybus matrix of size nxn
    row, col = y_bus.shape # gets row and column length of Ybus matrix
    for i in range(row):
        for j in range(col):
            # Iterative loop for checking each transmission line in the t_lines array
            for k in range(len(t_lines)):
                # Checks and assigns the off-diagonal elements based on what buses the transmission line connects
                if ((t_lines[k].start == buses[i] and t_lines[k].end == buses[j]) or
                        (t_lines[k].start == buses[j] and t_lines[k].end == buses[i])):
                    g_temp = round(t_lines[k].Gsh,3)
                    b_temp = round(t_lines[k].Bsh,3)
                    y_bus[i,j] = -1*complex(g_temp,b_temp)
                # Checks and assigns diagonal elements if the transmission line touches the relevant bus
                elif (t_lines[k].start == buses[i] or t_lines[k].end == buses[j]) and i ==j:
                    g_temp = round(t_lines[k].Gsh, 3)
                    b_temp = round(t_lines[k].Bsh, 3)
                    y_bus[i, j] += complex(g_temp,b_temp)

    return y_bus

def build_ybus_polar(y_bus_rect):
    """
    Converts a Ybus in rectangular complex format to a Ybus with complex numbers in polar form
    :param y_bus_rect: an nxn matrix of a Ybus with complex values in rectangular format
    :return: an nxn matrix of the input Ybus matrix with the complex numbers changed to polar form as a tuple
    """
    y_bus_polar = np.zeros(y_bus_rect.shape,tuple)
    row, col = y_bus_polar.shape
    for i in range(row):
        for j in range(col):
            y_bus_polar[i,j] = to_polar(y_bus_rect[i,j])
    return y_bus_polar

# Create the buses using the Bus data type from the given data
Alan = Bus("SL", 0.98, 0, 0, 0, 0, 0, 0)
Betty = Bus("PV", 1.00, 210, 50, 0, 0, 0, 1)
Clyde = Bus("PQ", 1.00, 0, 0, 110, 85, 150, 2)
Doug = Bus("PQ", 1.00, 0, 0, 100, 95, 50, 3)
Eve = Bus("PQ", 1.00, 0, 0, 150, 120, 0, 4)

# Creates the Transmission lines using the T_line data type from given data
AB = T_line(Alan, Betty, 0.009, 0.041, 0.000, 0.000, 125)
BE = T_line(Betty, Eve, 0.006, 0.037, 0.000, 0.000, 250)
AD = T_line(Alan, Doug, 0.007, 0.055, 0.000, 0.000, 200)
DE = T_line(Doug, Eve, 0.006, 0.045, 0.000, 0.000, 125)
DC = T_line(Doug, Clyde, 0.011, 0.061, 0.000, 0.000, 80)

CE = T_line(Clyde, Eve, 0.010, 0.051, 0.000, 0.000, 75)
Y_AA = AB.admittance + AD.admittance
Y_BB = AB.admittance + BE.admittance
Y_CC = DC.admittance + CE.admittance
Y_DD = AD.admittance + DE.admittance + DC.admittance

Y_EE = BE.admittance + DE.admittance + CE.admittance

y_bus_1 = np.array(
        [[  Y_AA,           -1*AB.admittance,           0,          -1*AD.admittance,   0               ],
         [  -1*AB.admittance,   Y_BB,                   0,              0,              -1*BE.admittance],
         [  0,                  0,                      Y_CC,       -1*DC.admittance,   -1*CE.admittance],
         [  -1*AD.admittance,   0,              -1*DC.admittance,      Y_DD,            -1*DE.admittance],
         [  0,              -1*BE.admittance,   -1*CE.admittance,   -1*DE.admittance,   Y_EE             ]])

y_bus_2 = build_ybus_rect([Alan, Betty, Clyde, Doug, Eve], [AB, BE, AD, DE, DC, CE])
print(y_bus_2)

test = build_ybus_polar(y_bus_2)
print(test)

