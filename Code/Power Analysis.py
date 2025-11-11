"""
EE 430 Power Analytical Methods of Power Systems - Fall 2025
Term Project - Newton-Raphson Algorithm
Joshua Consenz - 11/10/25

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
from copy import copy

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
    theta = round(theta,3)
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
        i_index = get_index(buses, i)
        for j in range(col):
            j_index = get_index(buses, j)
            # Iterative loop for checking each transmission line in the t_lines array
            for k in range(len(t_lines)):
                # Checks and assigns the off-diagonal elements based on what buses the transmission line connects
                if ((t_lines[k].start == buses[i_index] and t_lines[k].end == buses[j_index]) or
                        (t_lines[k].start == buses[j_index] and t_lines[k].end == buses[i_index])):
                    g_temp = round(t_lines[k].Gsh,3)
                    b_temp = round(t_lines[k].Bsh,3)
                    y_bus[i,j] = -1*complex(g_temp,b_temp)
                # Checks and assigns diagonal elements if the transmission line touches the relevant bus
                elif (t_lines[k].start == buses[i_index] or t_lines[k].end == buses[j_index]) and i_index == j_index:
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
    unknown[2] = value of unknown matrix
    Flat start conditions are delta = 0, voltage = 1
    :param buses: array of buses in system
    :return: unknown matrix using tuples
    """
    # TODO: recreate function to work without tuples
    delta = np.zeros(buses.shape, tuple)
    voltage = np.zeros(buses.shape, tuple)
    for bus in buses:
        if bus.type == "PV":
            delta[bus.index] = (bus.index, "d", bus.angle)
        elif bus.type == "PQ":
            delta[bus.index] = (bus.index, "d", bus.angle)
            voltage[bus.index] = (bus.index, "v", bus.volts)
    delta = delta[delta != 0]
    voltage = voltage[voltage != 0]
    unknown = np.concatenate((delta, voltage))
    return unknown

def build_mismatch(buses: np.ndarray) -> np.ndarray:
    """
    Builds the unknown matrix from initial conditions using tuples
    mismatch[0] = bus index
    mismatch[1] = "P" or "Q"
    mismatch[2] = value of mismatch matrix
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
    return round(P, 3)

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
    return round(Q,3)

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
    :return: The partial derivative of P_i with respect to d_i for Jacobian J1
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
    :return: The partial derivative of P_i with respect to V_i for Jacobian J2
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
    :return: The partial derivative of Q_i with respect to d_i for Jacobian J1
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
    :return: The partial derivative of Q_i with respect to V_i for Jacobian J2
    """
    bus_i = buses[get_index(buses,i)]
    v_i = bus_i.volts
    d_i = bus_i.angle
    y_ii = ybus[bus_i.index, bus_i.index][0]
    theta_ii = ybus[bus_i.index, bus_i.index][1]
    j4 = -1 * v_i * y_ii * math.sin(theta_ii)
    for bus_j in buses:
        v_j = bus_j.volts
        y_ij = ybus[bus_i.index, bus_j.index][0]
        d_j = bus_j.angle
        theta_ij = ybus[bus_i.index, bus_j.index][1]
        j4 += y_ij * v_j * math.sin(d_i - d_j - theta_ij)
    return round(j4, 3)

def get_jacobian_i_and_j(Jn: np.ndarray, indices: np.ndarray) -> np.ndarray:
    """
    Creates a matrix of the same shape as Jn, and returns a matrix i_and_j with the 1D array of indices converted
    into a 2D array with the bus indices for the partial derivative of Jn
    :param Jn:
    :param indices:
    :return:
    """
    i_and_j = np.zeros_like(Jn, tuple)
    for i in range(Jn.shape[0]):
        for j in range(Jn.shape[1]):
            i_and_j[i, j] = (indices[i], indices[j])
    return i_and_j

def build_J1(buses: np.ndarray, ybus: np.ndarray, n: int) -> np.ndarray:
    J1 = np.zeros((n-1,n-1))
    indices = np.zeros(len(buses))

    # Calculate indices the partial derivative needs
    for i in range(len(buses)):
        if buses[i].type != "SL":
            indices[i] = buses[i].index
        else:
            indices[i] = 99
    indices = indices[indices != 99]

    i_and_j = get_jacobian_i_and_j(J1, indices)

    # Assign values to J1 submatrix
    for i in range(J1.shape[0]):
        for j in range(J1.shape[1]):
            if i == j:
                J1[i,j] = J1_E(buses, ybus, i_and_j[i,j][0])
            elif i != j:
                bus_i = buses[get_index(buses, i_and_j[i,j][0])]
                bus_j = buses[get_index(buses, i_and_j[i,j][1])]
                J1[i,j] = J1_NE(bus_i, bus_j, ybus)

    return J1

def build_J2(buses: np.ndarray, ybus: np.ndarray, n: int, m: int) -> np.ndarray:
    """

    :param buses:
    :param ybus:
    :param n:
    :param m:
    :return:
    """
    J2 = np.zeros((n-1,n-1-m))
    indices = np.zeros(len(buses))
    for i in range(len(buses)):
        if buses[i].type != "SL":
            indices[i] = get_index(buses,i)
        else:
            indices[i] = 99
    indices = indices[indices != 99]

    i_and_j = get_jacobian_i_and_j(J2, indices)

    for i in range(J2.shape[0]):
        for j in range(J2.shape[1]):
            if i == j:
                J2[i,j] = J2_E(buses, ybus, i_and_j[i,j][0])
            elif i != j:
                bus_i = buses[get_index(buses, i_and_j[i,j][0])]
                bus_j = buses[get_index(buses, i_and_j[i,j][1])]
                J2[i,j] = J2_NE(bus_i, bus_j, ybus)

    return J2

def build_J3(buses: np.ndarray, ybus: np.ndarray, n: int, m: int) -> np.ndarray:
    """

    :param buses:
    :param ybus:
    :param n:
    :param m:
    :return:
    """
    J3 = np.zeros((n-1-m,n-1))

    indices = np.zeros(len(buses))
    for i in range(len(buses)):
        if buses[i].type != "SL":
            indices[i] = get_index(buses, i)
        else:
            indices[i] = 99
    indices = indices[indices != 99]

    i_and_j = get_jacobian_i_and_j(J3, indices)

    for i in range(J3.shape[0]):
        for j in range(J3.shape[1]):
            if i == j:
                J3[i, j] = J3_E(buses, ybus, i_and_j[i, j][0])
            elif i != j:
                bus_i = buses[get_index(buses, i_and_j[i, j][0])]
                bus_j = buses[get_index(buses, i_and_j[i, j][1])]
                J3[i, j] = J3_NE(bus_i, bus_j, ybus)

    return J3

def build_J4(buses: np.ndarray, ybus: np.ndarray, n: int, m: int) -> np.ndarray:
    """

    :param buses:
    :param ybus:
    :param n:
    :param m:
    :return:
    """
    J4 = np.zeros((n-1-m,n-1-m))

    indices = np.zeros(len(buses))
    for i in range(len(buses)):
        if buses[i].type == "PQ":
            indices[i] = get_index(buses, i)
        else:
            indices[i] = 99
    indices = indices[indices != 99]

    i_and_j = get_jacobian_i_and_j(J4, indices)

    for i in range(J4.shape[0]):
        for j in range(J4.shape[1]):
            if i == j:
                J4[i, j] = J4_E(buses, ybus, i_and_j[i, j][0])
            elif i != j:
                bus_i = buses[get_index(buses, i_and_j[i, j][0])]
                bus_j = buses[get_index(buses, i_and_j[i, j][1])]
                J4[i, j] = J4_NE(bus_i, bus_j, ybus)

    return J4

def create_jacobian(buses: np.ndarray, ybus: np.ndarray) -> np.ndarray:
    n = len(buses)  # number of buses
    m = 0  # number of PV buses
    for bus in buses:
        if bus.type == "PV":
            m += 1

    J1 = build_J1(buses, ybus, n)
    J2 = build_J2(buses, ybus, n, m)
    J3 = build_J3(buses, ybus, n, m)
    J4 = build_J4(buses, ybus, n, m)

    top_half = np.hstack([J1, J2])
    bottom_half = np.hstack([J3, J4])
    Jacobian = np.vstack([top_half, bottom_half])
    return Jacobian

"""
=============================
Newton-Raphson Algorithm
=============================
"""
def Newton_Raphson(buses: np.ndarray, tLines: np.ndarray, baseMVA: float, V_Tolerance: float) -> np.ndarray:
    """

    :param buses:
    :param tLines:
    :param baseMVA:
    :param V_Tolerance:
    :return:
    """

    #TODO: The whole tuple structure of my mismatch and unknown matrices is really bad
    # Really I just need to create an index array and a net power array for
    # Also I can't update the bus's net power with the calculated P and Q values
    # The difference is what goes in, not the calculated
    # As long as I update it again it should be fine, but I don't like using my buses to store temp variables

    #TODO: Probably want to restructure my bus class as well
    # it makes more sense to have the indices based on how the array is given to the function
    # currently I am fighting with myself by having a strictly defined bus index that can be placed anywhere in the array

    # Convert Bus Power Values to per unit
    for bus in buses:
        bus.Pgen = bus.Pgen/baseMVA
        bus.Pload = bus.Pload/baseMVA
        bus.Qgen = bus.Qgen/baseMVA
        bus.Qload = bus.Qload/baseMVA
        bus.Qcap = bus.Qcap/baseMVA


    "Step 1: Set up the Y_bus matrix"
    yBusRect = build_ybus_rect(buses, tLines)
    yBusPolar = build_ybus_polar(yBusRect)

    "Step 2: Impose Flat Start"
    for bus in buses:
        bus.__voltAngle__(0)
        bus.__netP__()
        bus.__netQ__()
    unknown_k = build_unknown(buses) # Create unknown matrix

    # Loop to determine convergence. Stops after n iterations in case convergence isn't reached
    #for k in range(10):

    "Step 3: Set up the mismatch matrix"
    mismatch_specified = build_mismatch(buses)

    # Creates a dummy copy of the buses to store values and use in calculations
    buses_copy = np.zeros_like(buses)
    for j in range(len(buses)):
        buses_copy[j] = buses[j]

    # Create calculated values for the mismatch matrix
    for val in mismatch_specified:
        index = get_index(buses_copy, val[0])
        if val[1] == "P":
            buses_copy[index].netP = calc_p_at_bus(buses_copy, yBusPolar, val[0])
        elif val[1] == "Q":
            buses_copy[index].netQ = calc_q_at_bus(buses_copy, yBusPolar, val[0])
    mismatch_calculated = build_mismatch(buses_copy)

    # Update buses with the new P and Q values using the calculated and specified matrices
    # and create the mismatch matrix
    mismatch = np.zeros(len(mismatch_specified))
    for j in range(len(mismatch_specified)):
        mismatch[j] = mismatch_specified[j][2] - mismatch_calculated[j][2]
        index = get_index(buses, mismatch_specified[j][0])
        if mismatch_specified[j][1] == "P":
            buses[index].netP = mismatch[j]
        elif mismatch_specified[j][1] == "Q":
            buses[index].netQ = mismatch[j]

    "Step 4: Create and fill in the Jacobian"
    Jacobian = create_jacobian(buses, yBusPolar)
    # TODO: verify this is correct with hand calculations for original input system

    "Update Unknown Matrix"
    J_inverse = np.linalg.inv(Jacobian)
    vals = np.zeros_like(unknown_k)
    for j in range(len(unknown_k)):
        vals[j] = unknown_k[j][2]
    unknown_k1 = np.linalg.matmul(J_inverse, mismatch) + vals

    # Round all values of unknown k+1 matrix
    for i in range(len(unknown_k1)):
        unknown_k1[i] = round(unknown_k1[i],3)

    # Assign voltage and angle values to buses
    for i in range(len(unknown_k)):
        index = get_index(buses, unknown_k[i][0])
        if unknown_k[i][1] == "d":
            buses[index].angle = unknown_k1[i]
        elif unknown_k[i][1] == "v":
            buses[index].volts = unknown_k1[i]


    return unknown_k1


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
unOrdered = np.zeros_like(busArray) # Designed to test index != buses position
un_tLineArray = np.zeros_like(tLineArray)
for i in range(len(busArray)):
    unOrdered[i] = busArray[i]
for i in range(len(tLineArray)):
    un_tLineArray[i] = tLineArray[i]

unOrdered = unOrdered[::-1]
un_tLineArray = un_tLineArray[::-1]

# Test system from HW 3
Uno = Bus("Uno", "SL", 1.00, 0, 0, 0, 0, 0, 0)
Dos = Bus("Dos", "PQ", 1.00, 0, 0, 0.9, 0.5, 0, 1)
Tres = Bus("Tres", "PV", 1.01, 1.3, 0, 0, 0, 0, 2)

UD = T_line(Uno, Dos, 0, 0.1, 0, 0, 1)
UT = T_line(Uno, Tres, 0, 0.25, 0, 0, 1)
DT = T_line(Dos, Tres, 0, 0.2, 0, 0, 1)

Tres_1 = Bus("Tres_1", "PV", 1.01, 1.3, 0, 0, 0, 0, 2)
Dos_1 = Bus("Dos_1", "PQ", 1.00, 0, 0, 0.9, 0.5, 0, 1)
Uno_1 = Bus("Uno_1", "SL", 1.00, 0, 0, 0, 0, 0, 0)

UD_1 = T_line(Uno_1, Dos_1, 0, 0.1, 0, 0, 1)
UT_1 = T_line(Uno_1, Tres_1, 0, 0.25, 0, 0, 1)
DT_1 = T_line(Dos_1, Tres_1, 0, 0.2, 0, 0, 1)

# HURRAY!!!!! This works and mirrors my homework problem
print(Newton_Raphson(np.array([Uno, Dos, Tres]), np.array([UD, UT, DT]), 1, 0.01))

print(Newton_Raphson(np.array([Tres_1, Uno_1, Dos_1]), np.array([UD_1, UT_1, DT_1]), 1, 0.01))

print(Newton_Raphson(busArray, tLineArray, baseMVA, V_Tolerance))

# These results differ, so my index catch isn't fully foolproof
print(Newton_Raphson(unOrdered, un_tLineArray, baseMVA, V_Tolerance))




