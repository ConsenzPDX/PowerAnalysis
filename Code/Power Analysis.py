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
from time import process_time_ns
from typing import Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Functions import *
from bus import Bus
from t_line import T_line


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

def build_unknown(buses: np.ndarray) -> np.ndarray:
    """
    Builds the unknown matrix in a flat start using tuples.
    unknown[0] = bus index
    unknown[1] = "delta" or "voltage"
    unknown[2] = value of unkonwn matrix
    Flat start conditions are delta = 0, voltage = 1
    :param buses: array of buses in system
    :return: unknown matrix using tuples
    """
    # TODO: recreate function to work without tuples
    delta = np.zeros(buses.shape, tuple)
    voltage = np.zeros(buses.shape, tuple)
    for bus in buses:
        if bus.type == "PV":
            delta[bus.index] = (bus.index, "delta", bus.angle)
        elif bus.type == "PQ":
            delta[bus.index] = (bus.index, "delta", bus.angle)
            voltage[bus.index] = (bus.index, "voltage", bus.volts)
    delta = delta[delta != 0]
    voltage = voltage[voltage != 0]
    unknown = np.concatenate((delta, voltage))
    return unknown

def build_mismatch(buses: np.ndarray) -> np.ndarray:
    """
    Builds the unknown matrix from initial conditions using tuples
    mismatch[0] = bus index
    mismatch[1] = "P" or "Q"
    mismatch[2] = value of mismatch amtrix
    :param buses: array of buses in system
    :return: mismatch matrix using tuples
    """
    # TODO: recreate function to work without tuples
    P = np.zeros(buses.shape, tuple)
    Q = np.zeros(buses.shape, tuple)
    for bus in buses:
        if bus.type == "PV":
            P[bus.index] = (bus.index, "P", bus.netP)
        elif bus.type == "PQ":
            P[bus.index] = (bus.index, "P", bus.netP)
            Q[bus.index] = (bus.index, "Q", bus.netQ)
    P = P[P != 0]
    Q = Q[Q != 0]
    mismatch = np.concatenate((P, Q))
    return mismatch

def calc_p_at_bus(buses: np.ndarray, ybus: np.ndarray, i: int) -> float:
    # TODO: comment this function
    bus_i = buses[get_index(buses,i)]
    v_i = bus_i.volts
    d_i = bus_i.angle
    P = 0
    for bus in buses:
        y_ij = ybus[bus_i.index, bus.index][0]
        v_j = bus.volts
        d_j = bus.angle
        theta_ij = ybus[bus_i.index, bus.index][1]
        P += v_i * y_ij * v_j * math.cos(theta_ij + d_j - d_i)
    return P

def calc_q_at_bus(buses: np.ndarray, ybus: np.ndarray, i: int) -> float:

    # TODO: comment and fill out docstring
    bus_i = buses[get_index(buses, i)]
    v_i = bus_i.volts
    d_i = bus_i.angle
    Q = 0
    for bus in buses:
        y_ij = ybus[bus_i.index, bus.index][0]
        v_j = bus.volts
        d_j = bus.angle
        theta_ij = ybus[bus_i.index, bus.index][1]
        Q += y_ij * v_j * math.sin(theta_ij + d_j - d_i)
    Q *= -1 * v_i
    return Q

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
    :param ybus: Ybus matrix of the system in polar form, which contains the necessary admittance values
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
    :param ybus: Ybus matrix of the system in polar form, to access specific admittance values
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
    :param ybus: Ybus matrix of the system in polar form, to access specific admittance values
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
    :param ybus: Ybus matrix of the system in polar form, to access specific admittance values
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
    :param ybus: Ybus matrix of the system in polar form
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
    :param ybus: Ybus matrix of the system in polar form
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
    :param ybus: Ybus matrix of the system in polar form
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
    :param ybus: Ybus matrix of the system in polar form
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

"""
=============================
Newton-Raphson Algorithm
=============================
"""
def Newton_Raphson(system: np.ndarray, tLines: np.ndarray, baseMVA: float, V_Tolerance: float):
    """

    :param system:
    :param tLines:
    :param baseMVA:
    :param V_Tolerance:
    :return:
    """

    #TODO: The whole tuple structure of my mismatch and unknown matrices is really bad
    # Really I just need to create an index array and a net power array for
    # Also I can't update the bus's net power with the calcualted P and Q values
    # The difference is what goes in, not the calculated
    # As long as I update it again it should be fine, but I don't like using my buses to store temp variables

    #TODO: Probably want to restructure my bus class as well
    # it makes more sense to have the indeces based on how the array is given to the function
    # currently I am fighting with myself by having a strictly defined bus index that can be placed anywhere in the array

    # Creates copy of system to operate on
    buses = np.copy(system)

    # Convert Bus Power Values to per unit
    # for bus in buses:
    #     bus.Pgen = bus.Pgen/baseMVA
    #     bus.Pload = bus.Pload/baseMVA
    #     bus.Qgen = bus.Qgen/baseMVA
    #     bus.Qload = bus.Qload/baseMVA
    #     bus.Qcap = bus.Qcap/baseMVA


    "Step 1: Set up the Y_bus matrix"
    yBusRect = build_ybus_rect(buses, tLines)
    yBusPolar = build_ybus_polar(yBusRect)

    "Step 2: Impose Flat Start"
    for bus in buses:
        bus.__voltAngle__(0)
        bus.__netP__()
        bus.__netQ__()
    unknown = build_unknown(buses) # Create unknown matrix



    #for k in range(10):

    "Step 3: Set up the mismatch matrix"
    mismatch_k = build_mismatch(buses)

    for val in mismatch_k:
        index = get_index(buses, val[0])
        if val[1] == "P":
            buses[index].netP = calc_p_at_bus(buses, yBusPolar, val[0])
        elif val[1] == "Q":
            buses[index].netQ = calc_q_at_bus(buses, yBusPolar, val[0])
    mismatch_k1 = build_mismatch(buses)

    mismatch = np.zeros(len(mismatch_k))
    for j in range(len(mismatch)):
        mismatch[j] = mismatch_k[j][2] - mismatch_k1[j][2]
    print(mismatch_k)
    print(mismatch_k1)
    print(mismatch)



"""
==============================
Initial System Conditions 
==============================
"""

# Base values and tolerance of the system
baseMVA = 100
V_Tolerance = 0.05

# Create the buses using the Bus data type from the given data
Alan = Bus("Alan", "SL", 0.98, 0, 0, 0, 0, 0, 0)
Betty = Bus("Betty", "PV", 1.00, 210, 50, 0, 0, 0, 1)
Clyde = Bus("Clyde", "PQ", 1.00, 0, 0, 110, 85, 150, 2)
Doug = Bus("Doug", "PQ", 1.00, 0, 0, 100, 95, 50, 3)
Eve = Bus("Eve", "PQ", 1.00, 0, 0, 150, 120, 0, 4)

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
unoredered = np.array([Eve, Doug, Clyde, Betty, Alan]) # Desgined to test index != busses position

print(Newton_Raphson(busArray, tLineArray, baseMVA, V_Tolerance))




