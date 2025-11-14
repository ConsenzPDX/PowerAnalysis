"""
EE 430 Power Analytical Methods of Power Systems - Fall 2025
Term Project - Newton-Raphson Algorithm
Joshua Consenz - 11/14/25

System conditions and inputs to the system to test the algorithms I made
"""
from numpy.ma.core import zeros_like

from newton_raphson import *

"""
==============================
Five Bus PQ System
==============================
"""

# Base values and tolerance of the system
name = "FiveBus_PQ"
baseMVA = 100
V_Tolerance = 0.05

# Create the buses using the Bus data type from the given data
Alan = Bus("Alan", "SL", 0.98, 0, 0, 0, 0, 0)
Betty = Bus("Betty", "PV", 1.00, 210, 50, 0, 0, 10)
Clyde = Bus("Clyde", "PQ", 1.00, 0, 0, 110, 85, 150)
Doug = Bus("Doug", "PQ", 1.00, 0, 0, 100, 95, 50)
Eve = Bus("Eve", "PQ", 1.00, 0, 0, 150, 120, 0)

# Creates the Transmission lines using the T_line data type from given data
AB = T_line(Alan, Betty, 0.009, 0.041, 0.000, 0.000, 125)
BE = T_line(Betty, Eve, 0.006, 0.037, 0.000, 0.000, 250)
AD = T_line(Alan, Doug, 0.007, 0.055, 0.000, 0.000, 200)
DE = T_line(Doug, Eve, 0.006, 0.045, 0.000, 0.000, 125)
DC = T_line(Doug, Clyde, 0.011, 0.061, 0.000, 0.000, 80)
CE = T_line(Clyde, Eve, 0.010, 0.051, 0.000, 0.000, 75)

# Collect buses and transmission lines into arrays to pass to looping function
busArray = np.array([Alan, Betty, Clyde, Doug, Eve])
tLineArray = np.array([AB, BE, AD, DE, DC, CE])

FiveBus_PQ = Newton_Raphson(busArray, tLineArray, baseMVA, V_Tolerance, 100, 0.001, name)

prnt = np.zeros_like(FiveBus_PQ)
for i in range(len(FiveBus_PQ)):
    prnt[i] = round(FiveBus_PQ[i], 3)
print("Final Unknown matrix:", prnt)

ybusRect = build_ybus_rect(busArray, tLineArray)
# print(build_ybus_polar(ybusRect))

"""
=========================================
Homework 3 Problem 1 Test Case
=========================================
"""

# Uno = Bus("Uno", "SL", 1.00, 0, 0, 0, 0, 0)
# Dos = Bus("Dos", "PQ", 1.00, 0, 0, 0.9, 0.5, 0)
# Tres = Bus("Tres", "PV", 1.01, 1.3, 0, 0, 0, 1.0)
#
# UD = T_line(Uno, Dos, 0, 0.1, 0, 0, 1)
# UT = T_line(Uno, Tres, 0, 0.25, 0, 0, 1)
# DT = T_line(Dos, Tres, 0, 0.2, 0, 0, 1)
# #
# Tres_1 = Bus("Tres_1", "PV", 1.01, 1.3, 0, 0, 0, 1.0)
# Dos_1 = Bus("Dos_1", "PQ", 1.00, 0, 0, 0.9, 0.5, 0)
# Uno_1 = Bus("Uno_1", "SL", 1.00, 0, 0, 0, 0, 0)
#
# UD_1 = T_line(Uno_1, Dos_1, 0, 0.1, 0, 0, 1)
# UT_1 = T_line(Uno_1, Tres_1, 0, 0.25, 0, 0, 1)
# DT_1 = T_line(Dos_1, Tres_1, 0, 0.2, 0, 0, 1)
# #
# # # HURRAY!!!!! This works and mirrors my homework problem
# print("\n")
# print(Newton_Raphson(np.array([Uno, Dos, Tres]), np.array([UD, UT, DT]), 1, 0.05, 100, "HW3"))
# print("\n")
# print(Newton_Raphson(np.array([Tres_1, Uno_1, Dos_1]), np.array([UD_1, UT_1, DT_1]), 1, 0.01))

"""
======================
Class Example
======================
"""

