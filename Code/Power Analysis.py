"""
EE 430 Power Analytical Methods of Power Systems - Fall 2025
Term Project - Newton-Raphson Algorithm
Joshua Consenz - 11/6/25

Creates a five bus system with six transmission lines, and implements the Newton-Raphson algorithm to solve the system
from an initial state to a steady state.
"""

import math
import cmath
import string
from typing import Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Functions import *
from bus import Bus
from t_line import T_line

# Base values and tolerance of the system
baseMVA = 100
V_Tolerance = 0.05

def to_polar(y : complex) -> tuple:
    """
    Takes a complex number, y, and converts it into a tuple with magnitude and phase in degrees
    :param y: input complex number in rectangular format
    :return: a tuple of the complex number as magnitude and phase in degrees
    """
    r, theta = cmath.polar(y)
    r = round(r,3)
    theta = round(math.degrees(theta),3)
    return r, theta

def build_ybus_rect(buses: np.ndarray, t_lines: np.ndarray) -> np.ndarray:
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

def build_ybus_polar(y_bus_rect: np.ndarray) -> np.ndarray:
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
     # element [0] in the tuple is magnitude, and [1] is angle
    return y_bus_polar

def get_index(buses: np.ndarray, i: int) -> int:
    """
    Function to find bus i in buses array regardless of
    :param buses: array of system buses
    :param i: bus index we are looking for
    :return: bus i's position in the buses array
    """
    index = 0
    for j in range(len(buses)):
        if buses[j].index == i:
            index = j
            break # Exit the loop when we've found the right value
    return index

"""
====================================
Jacobian Matrix Element Functions
====================================
"""
def J1_NE(bus1: Bus, bus2: Bus, ybus: np.ndarray) -> float:
    """
    Jacobian Value J1 when i != j
    NE = Not Equals
    J1 = P_i/d_j = V_i * Y_ij * V_j * sin(d_i - d_j - theta_ij)
    :param bus1: Bus i we are calculating the partial derivative of real power for
    :param bus2: Bus j whose angle we are using with respect to, for the partial derivative
    :param ybus: Ybus matrix of the system, which contains the necessary admittance values
    :return: the partial derivative of P_i with respect to d_j for Jacobian J1
    """
    v_i = bus1.volts
    v_j = bus2.volts
    y_ij = ybus[bus1.index,bus2.index][0]
    d_i = bus1.angle
    d_j = bus2.angle
    theta_ij = ybus[bus1.index,bus2.index][1]
    j1 = v_i * y_ij * v_j * math.sin(d_i - d_j - theta_ij)
    return round(j1,3)

def J2_NE(bus1: Bus, bus2: Bus, ybus: np.ndarray) -> float:
    """
    Jacobian Value J2 when i != j
    NE = Not Equals
    J2 = P_i/V_j = V_i * Y_ij * cos(d_i - d_j - theta_ij)
    :param bus1: Bus i we are calculating the partial derivative of real power for
    :param bus2: Bus j whose voltage we are using with respect to, for the partial derivative
    :param ybus: Ybus matrix of the system, to access specific admittance values
    :return: the partial derivative of P_i with respect to V_j for Jacobian J2
    """
    v_i = bus1.volts
    y_ij = ybus[bus1.index, bus2.index][0]
    d_i = bus1.angle
    d_j = bus2.angle
    theta_ij = ybus[bus1.index, bus2.index][1]
    j2 = v_i * y_ij * math.cos(d_i - d_j - theta_ij)
    return round(j2,3)

def J3_NE(bus1: Bus, bus2: Bus, ybus: np.ndarray) -> float:
    """
    Jacobian Value J3 when i != j
    NE = Not Equals
    J3 = Q_i/d_j = -1 * V_i * Y_ij * V_j * cos(d_i - d_j - theta_ij)
    :param bus1: Bus i we are calculating the partial derivative of reactive power for
    :param bus2: Bus j whose angle we are using with respect to, for the partial derivative
    :param ybus: Ybus matrix of the system, to access specific admittance values
    :return: the partial derivative of Q_i with respect to d_j for Jacobian J3
    """
    v_i = bus1.volts
    v_j = bus2.volts
    y_ij = ybus[bus1.index, bus2.index][0]
    d_i = bus1.angle
    d_j = bus2.angle
    theta_ij = ybus[bus1.index, bus2.index][1]
    j3 = -1 * v_i * y_ij * v_j * math.cos(d_i - d_j - theta_ij)
    return round(j3, 3)

def J4_NE(bus1: Bus, bus2: Bus, ybus: np.ndarray) -> float:
    """
    Jacobian Value J4 when i != j
    NE = Not Equals
    J4 = Q_i/V_j = V_i * Y_ij * sin(d_i - d_j - theta_ij)
    :param bus1: Bus i we are calculating the partial derivative of reactive power for
    :param bus2: Bus j whose voltage we are using with respect to, for the partial derivative
    :param ybus: Ybus matrix of the system, to access specific admittance values
    :return: the partial derivative of Q_i with respect to V_j for Jacobian J4
    """
    v_i = bus1.volts
    y_ij = ybus[bus1.index, bus2.index][0]
    d_i = bus1.angle
    d_j = bus2.angle
    theta_ij = ybus[bus1.index, bus2.index][1]
    j4 = v_i * y_ij * math.sin(d_i - d_j - theta_ij)
    return round(j4,3)

def J1_E(buses: np.ndarray, ybus: np.ndarray, i: int) -> float:
    """
    Jacobian Value J1 when i == j
    E = equals
    J1 = P_i/d_i = -1* V_i * sum(j = 0 -> len(buses), j != i, Y_ij * v_j * sin(d_i - d_j - theta_ij)
    :param buses: array of buses in the system
    :param ybus: Ybus matrix of the system
    :param i: index of the values we want to calculate the partial derivative of
    :return: The poartial derivative of P_i with respect to d_i for Jacobian J1
    """
    bus_i = buses[get_index(buses,i)]
    v_i = bus_i.volts
    d_i = bus_i.angle
    j1 = 0
    for bus_j in buses:
        if bus_j.index != i:
            v_j = bus_j.volts
            y_ij = ybus[bus_i.index, bus_j.index][0]
            d_j = bus_j.angle
            theta_ij = ybus[bus_i.index, bus_j.index][1]
            j1 += y_ij * v_j * math.sin(d_i - d_j - theta_ij)
    j1 *= -1 * v_i
    return round(j1,3)

def J2_E(buses: np.ndarray, ybus: np.ndarray, i: int) -> float:
    """
    Jacobian Value J2 when i == j
    E = equals
    J1 = P_i/V_i = V_i * Y_ii * cos(theta_ii) + sum(j = 0 -> len(buses), Y_ij * v_j * cos(d_i - d_j - theta_ij)
    :param buses: array of buses in the system
    :param ybus: Ybus matrix of the system
    :param i: index of the values we want to calculate the partial derivative of
    :return: The poartial derivative of P_i with respect to V_i for Jacobian J2
    """
    bus_i = buses[get_index(buses,i)]
    v_i = bus_i.volts
    d_i = bus_i.angle
    y_ii = ybus[bus_i.index, bus_i.index][0]
    theta_ii = ybus[bus_i.index, bus_i.index][1]
    j2 = v_i * y_ii * math.cos(theta_ii)
    for bus_j in buses:
        v_j = bus_j.volts
        y_ij = ybus[bus_i.index, bus_j.index][0]
        d_j = bus_j.angle
        theta_ij = ybus[bus_i.index, bus_j.index][1]
        j2 += y_ij * v_j * math.cos(d_i - d_j - theta_ij)
    return round(j2,3)

def J3_E(buses: np.ndarray, ybus: np.ndarray, i: int) -> float:
    """
    Jacobian Value J3 when i == j
    E = equals
    J3 = Q_i/d_i = V_i * sum(j = 0 -> len(buses), j != i, Y_ij * v_j * cos(d_i - d_j - theta_ij)
    :param buses: array of buses in the system
    :param ybus: Ybus matrix of the system
    :param i: index of the values we want to calculate the partial derivative of
    :return: The poartial derivative of Q_i with respect to d_i for Jacobian J1
    """
    bus_i = buses[get_index(buses,i)]
    v_i = bus_i.volts
    d_i = bus_i.angle
    j3 = 0
    for bus_j in buses:
        if bus_j.index != i:
            v_j = bus_j.volts
            y_ij = ybus[bus_i.index, bus_j.index][0]
            d_j = bus_j.angle
            theta_ij = ybus[bus_i.index, bus_j.index][1]
            j3 += y_ij * v_j * math.cos(d_i - d_j - theta_ij)
    j3 *= v_i
    return round(j3, 3)

def J4_E(buses: np.ndarray, ybus: np.ndarray, i: int) -> float:
    """
    Jacobian Value J4 when i == j
    E = equals
    J4 = Q_i/V_i = -1 * V_i * Y_ii * sin(theta_ii) + sum(j = 0 -> len(buses), Y_ij * v_j * sin(d_i - d_j - theta_ij)
    :param buses: array of buses in the system
    :param ybus: Ybus matrix of the system
    :param i: index of the values we want to calculate the partial derivative of
    :return: The poartial derivative of Q_i with respect to V_i for Jacobian J2
    """
    bus_i = buses[get_index(buses,i)]
    v_i = bus_i.volts
    d_i = bus_i.angle
    y_ii = ybus[bus_i.index, bus_i.index][0]
    theta_ii = ybus[bus_i.index, bus_i.index][1]
    j4 = -1 * v_i * y_ii * math.cos(theta_ii)
    for bus_j in buses:
        v_j = bus_j.volts
        y_ij = ybus[bus_i.index, bus_j.index][0]
        d_j = bus_j.angle
        theta_ij = ybus[bus_i.index, bus_j.index][1]
        j4 += y_ij * v_j * math.sin(d_i - d_j - theta_ij)
    return round(j4, 3)

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

busArray = np.array([Alan, Betty, Clyde, Doug, Eve])
tLineArray = np.array([AB, BE, AD, DE, DC, CE])

unoredered = np.array([Eve, Doug, Clyde, Betty, Alan])

yBusRect = build_ybus_rect(busArray, tLineArray)
yBusPolar = build_ybus_polar(yBusRect)

# Impose Flat Start
Alan.__voltAngle__(0)
Betty.__voltAngle__(0)
Clyde.__voltAngle__(0)
Doug.__voltAngle__(0)
Eve.__voltAngle__(0)

print(J1_NE(Clyde, Doug, yBusPolar))
print(J2_NE(Clyde, Doug, yBusPolar))
print(J3_NE(Clyde, Doug, yBusPolar))
print(J4_NE(Clyde, Doug, yBusPolar))
print(J1_E(busArray, yBusPolar, 1))
print(J1_E(unoredered, yBusPolar, 1))
print(J2_E(busArray, yBusPolar, 2))
print(J2_E(unoredered, yBusPolar, 2))
print(J3_E(busArray, yBusPolar, 3))
print(J3_E(unoredered, yBusPolar, 3))
print(J4_E(busArray, yBusPolar, 4))
print(J4_E(unoredered, yBusPolar, 4))

