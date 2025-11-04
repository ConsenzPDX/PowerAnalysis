"""
EE 430 Power Analytical Methods of Power Systems - Fall 2025
Term Project - Newton-Raphson Algorithm
Joshua Consenz - 11/3/25

Creates a five bus system with six transmission lines, and implements the Newton-Raphson algorithm to solve the system
from an initial state to a steady state.
"""


import string
import numpy as np
import matplotlib.pyplot as plt

from Functions import *
from bus import Bus
from t_line import T_line

# Base values and tolerance of the system
baseMVA = 100
V_Tolerance = 0.05

# Create the buses using the Bus data type from the given data
Alan = Bus("SL", 0.98, 0, 0, 0, 0, 0)
Betty = Bus("PV", 1.00, 210, 50, 0, 0, 0)
Clyde = Bus("PQ", 1.00, 0, 0, 110, 85, 150)
Doug = Bus("PQ", 1.00, 0, 0, 100, 95, 50)
Eve = Bus("PQ", 1.00, 0, 0, 150, 120, 0)

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

y_bus = np.array(
        [[  Y_AA,           -1*AB.admittance,           0,          -1*AD.admittance,   0               ],
         [  -1*AB.admittance,   Y_BB,                   0,              0,              -1*BE.admittance],
         [  0,                  0,                      Y_CC,       -1*DC.admittance,   -1*CE.admittance],
         [  -1*AD.admittance,   0,              -1*DC.admittance,      Y_DD,            -1*DE.admittance],
         [  0,              -1*BE.admittance,   -1*CE.admittance,   -1*DE.admittance,   Y_EE             ]])


print(diag(y_bus))
print(off_diag(y_bus))
print(col(y_bus))
print(row(y_bus))


