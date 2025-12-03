"""
EE 430 Power Analytical Methods of Power Systems - Fall 2025
Term Project - Newton-Raphson Algorithm
Joshua Consenz - 11/19/25

System conditions and inputs to the system to test the algorithms I made
"""

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
criterion = 1e-15
FiveBus_PQ = Newton_Raphson(busArray, tLineArray, baseMVA, V_Tolerance, 100, criterion, name)

prnt = np.zeros_like(FiveBus_PQ)
for i in range(len(FiveBus_PQ)):
    prnt[i] = round(FiveBus_PQ[i], 3)
print("Final Unknown matrix:", FiveBus_PQ)