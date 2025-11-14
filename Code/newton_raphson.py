"""
EE 430 Power Analytical Methods of Power Systems - Fall 2025
Term Project - Newton-Raphson Algorithm
Joshua Consenz - 11/14/25

Creates a five bus system with six transmission lines, and implements the Newton-Raphson algorithm to solve the system
from an initial state to a steady state.
"""

import math
import cmath
import string
import time
from typing import Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy

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
    theta = theta
    return r, theta

# Function is useless, kept for posterity
# def get_index(buses: np.ndarray, i: int) -> int:
#     """
#     Function to find bus i in buses array regardless of
#     :param buses: array of system buses
#     :param i: bus index we are looking for
#     :return: bus i's position in the buses array
#     """
#     index = 0
#     for j in range(len(buses)):
#         if buses[j].index == i:
#             index = j
#             break # Exit the loop when we've found the right value
#     return index

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
                elif (t_lines[k].start == buses[i] or t_lines[k].end == buses[j]) and i == j:
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
    # TODO: check ybus is implemented correctly
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
    """
    Calculates the Power at a specific bus, i, based on its connections with other buses
    P_i = |V_i| * sum(j = 1 -> n, |Y_ij| * |V_j| * cos(theta_ij - d_j - d_i)
    :param buses: system information of the buses
    :param ybus: Ybus matrix
    :param i: bus we are interested in
    :return: Power at bus i
    """
    bus_i = buses[i]
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
    """
    Calculates Q at bus, i, based on the system and connections of the bus
    Q_i = -1 * |V_i| * sum(j = 1 -> n, |Y_ij| * |V_j| * sin(theta_ij - d_j - d_i)
    :param buses: system information of the buses
    :param ybus: Ybus matrix
    :param i: bus we are interested in
    :return: Q at bus, i
    """
    bus_i = buses[i]
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
    bus_i = buses[i]
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
    bus_i = buses[i]
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
    bus_i = buses[i]
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
    bus_i = buses[i]
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
    :param Jn: The empty array of the submatrix for some n
    :param indices: The indices of the buses that are needed for the Jacobian element
    :return: 2D array where each index of the partial derivative is mapped to its position in the jacobian submatrix
    """
    i_and_j = np.zeros_like(Jn, tuple)
    for i in range(Jn.shape[0]):
        for j in range(Jn.shape[1]):
            i_and_j[i, j] = (int(indices[i]), int(indices[j]))
    return i_and_j

def build_J1(buses: np.ndarray, ybus: np.ndarray, n: int) -> np.ndarray:
    """
    Builds the Jacobian submatrix J1
    :param buses: System information of buses
    :param ybus: Ybus matrix of the system in polar form
    :param n: number of buses in the system
    :return: J1 submatrix
    """
    J1 = np.zeros((n-1,n-1))
    indices = np.zeros(len(buses), int)

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
                bus_i = buses[i_and_j[i,j][0]]
                bus_j = buses[i_and_j[i,j][1]]
                J1[i,j] = J1_NE(bus_i, bus_j, ybus)

    return J1

def build_J2(buses: np.ndarray, ybus: np.ndarray, n: int, m: int) -> np.ndarray:
    """
    Builds the Jacobian submatrix J2
    :param buses: System information
    :param ybus: Ybus Matrix
    :param n: number of buses
    :param m: number of PV buses
    :return: J2 submatrix
    """
    J2 = np.zeros((n-1,n-1-m))
    indices = np.zeros(len(buses))
    for i in range(len(buses)):
        if buses[i].type != "SL":
            indices[i] = i
        else:
            indices[i] = 99
    indices = indices[indices != 99]

    i_and_j = get_jacobian_i_and_j(J2, indices)

    for i in range(J2.shape[0]):
        for j in range(J2.shape[1]):
            if i == j:
                J2[i,j] = J2_E(buses, ybus, i_and_j[i,j][0])
            elif i != j:
                bus_i = buses[i_and_j[i,j][0]]
                bus_j = buses[i_and_j[i,j][1]]
                J2[i,j] = J2_NE(bus_i, bus_j, ybus)

    return J2

def build_J3(buses: np.ndarray, ybus: np.ndarray, n: int, m: int) -> np.ndarray:
    """
    Builds the Jacobian submatrix J3
    :param buses: System information
    :param ybus: Ybus Matrix
    :param n: number of buses
    :param m: number of PV buses
    :return: submatrix J3
    """
    J3 = np.zeros((n-1-m,n-1))

    indices = np.zeros(len(buses))
    for i in range(len(buses)):
        if buses[i].type != "SL":
            indices[i] = i
        else:
            indices[i] = 99
    indices = indices[indices != 99]

    i_and_j = get_jacobian_i_and_j(J3, indices)

    for i in range(J3.shape[0]):
        for j in range(J3.shape[1]):
            if i == j:
                J3[i, j] = J3_E(buses, ybus, i_and_j[i, j][0])
            elif i != j:
                bus_i = buses[i_and_j[i, j][0]]
                bus_j = buses[i_and_j[i, j][1]]
                J3[i, j] = J3_NE(bus_i, bus_j, ybus)

    return J3

def build_J4(buses: np.ndarray, ybus: np.ndarray, n: int, m: int) -> np.ndarray:
    """
    Builds the Jacobian submatrix J4
    :param buses: System information
    :param ybus: Ybus Matrix
    :param n: number of buses
    :param m: number of PV buses
    :return: J4 submatrix
    """
    J4 = np.zeros((n-1-m,n-1-m))

    indices = np.zeros(len(buses))
    for i in range(len(buses)):
        if buses[i].type == "PQ":
            indices[i] = i
        else:
            indices[i] = 99
    indices = indices[indices != 99]

    i_and_j = get_jacobian_i_and_j(J4, indices)

    for i in range(J4.shape[0]):
        for j in range(J4.shape[1]):
            if i == j:
                J4[i, j] = J4_E(buses, ybus, i_and_j[i, j][0])
            elif i != j:
                bus_i = buses[i_and_j[i, j][0]]
                bus_j = buses[i_and_j[i, j][1]]
                J4[i, j] = J4_NE(bus_i, bus_j, ybus)

    return J4

def create_jacobian(buses: np.ndarray, ybus: np.ndarray) -> np.ndarray:
    """
    Creates the Jacobian matrix by using submatrix subroutines, and merging the submatrices together
    :param buses: Array of buses in the system
    :param ybus: Ybus matrix
    :return: Complete Jacobian matrix
    """
    n = len(buses)  # number of buses
    m = 0  # number of PV buses
    for bus in buses:
        if bus.type == "PV":
            m += 1
    # TODO: check
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
def Newton_Raphson(buses: np.ndarray, tLines: np.ndarray, base_mva: float, vTolerance: float, iterations = 10, name = "System") -> np.ndarray:
    """
    Newton-Raphson Algorithm designed to operate on a system of buses and transmission lines.
    The buses and Transmission lines are designed to use the Bus and T_line classes to build them.
    The Algorithm iterates a number of times equal to the iteration parameter, and checks for system convergence based
    on the Tolerance parameter.
    :param buses: Numpy Array of the buses that describe the system being analyzed
    :param tLines: Numpy Array of the transmission lines that connect the busses in the system
    :param base_mva: Base real power value to convert all P and Q values into per unit values
    :param vTolerance: Tolerance to check is an iteration has converged to a solution
    :param iterations: The maximum number of iterations allowed the program, if unspecified iterations = 10
    :param name: Name of the system
    :return: The unknown matrix after the system has converged, or the maximum number of iterations has been reached
    """

    # Set start time, to calculate how long the program needs to run for
    start_time = time.time()

    # Convergence criterion
    criterion = 0.01

    # Print Start of Algorithm message
    print(f"Beginning Newton-Raphson Power Flow Solution for {name}")

    # Set the index of each bus with the order it was fed into the system
    for i in range(len(buses)):
        buses[i].__setIndex__(int(i))

    # Convert Bus Power Values to per unit
    for bus in buses:
        bus.Pgen = bus.Pgen / base_mva
        bus.Pload = bus.Pload / base_mva
        bus.Qgen = bus.Qgen / base_mva
        bus.Qload = bus.Qload / base_mva
        bus.Qcap = bus.Qcap / base_mva


    "Step 1: Set up the Y_bus matrix"
    yBusRect = build_ybus_rect(buses, tLines) # Build the ybus using the complex numbers in rectangular form
    yBusPolar = build_ybus_polar(yBusRect) # Converts the rectangular complex numbers to a tuple representing polar form

    "Step 2: Impose Flat Start"
    for bus in buses:
        # If the bus is a PQ bus, set the voltage to flat start value, v = 1
        if bus.type == "PQ":
                bus.__setVoltage__(1)
        bus.__voltAngle__(0) # Set the voltage angle of the bus to 0 radians
        bus.__netP__() # Calculates the net real Power based on the generated and demanded power of the bus
        bus.__netQ__() # Calculates the net reactive power of the bus using the generated and demanded values

    # Boolean that is checked every for loop to see if system has converged
    converged = False

    # Loop to determine convergence. Stops after iterations in case convergence isn't reached
    for k in range(iterations):

        unknown_k = build_unknown(buses)  # Create unknown matrix from current conditions

        "Step 3: Set up the mismatch matrix"
        # Create the specified half of the mismatch matrix
        mismatch_specified = build_mismatch(buses)

        # Creates a dummy copy of the buses to store values and use in calculations
        buses_copy = np.zeros_like(buses)
        for j in range(len(buses)):
            buses_copy[j] = copy.deepcopy(buses[j])

        # Create calculated values for the mismatch matrix
        for val in mismatch_specified:
            index = val[0]
            if val[1] == "P":
                buses_copy[index].netP = calc_p_at_bus(buses_copy, yBusPolar, val[0])
            elif val[1] == "Q":
                buses_copy[index].netQ = calc_q_at_bus(buses_copy, yBusPolar, val[0])
        # Build the calculated portion of the mismatch matrix
        mismatch_calculated = build_mismatch(buses_copy)

        # Update buses with the new P and Q values using the calculated and specified matrices
        # and create the complete mismatch matrix
        mismatch = np.zeros(len(mismatch_specified))
        for j in range(len(mismatch_specified)):
            # Subtract specified - calculated for each bus
            mismatch[j] = mismatch_specified[j][2] - mismatch_calculated[j][2]
            index = mismatch_specified[j][0]
            # TODO: Talk to Midrar in office hours about if I update the P and Q of the buses during each iteration
            # Now update the buses with the new net power values
            # if mismatch_specified[j][1] == "P":
            #     buses[index].netP = mismatch[j]
            # elif mismatch_specified[j][1] == "Q":
            #     buses[index].netQ = mismatch[j]

        "Step 4: Create and fill in the Jacobian"
        Jacobian = create_jacobian(buses, yBusPolar)
        # TODO: verify this is correct with hand calculations for original input system

        "Update Unknown Matrix"
        J_inverse = np.linalg.inv(Jacobian) # Invert the Jacobian for the calculation
        # Create a matrix of only values for the mismatch matrix
        vals = np.zeros_like(unknown_k)
        for j in range(len(unknown_k)):
            vals[j] = unknown_k[j][2]
        # Calculate the new unknown matrix values
        unknown_k1 = np.linalg.matmul(J_inverse, mismatch) + vals

        # Round all values of unknown k+1 matrix
        for i in range(len(unknown_k1)):
            unknown_k1[i] = round(unknown_k1[i],3)

        # Assign voltage and angle values to buses
        for i in range(len(unknown_k)):
            index = unknown_k[i][0]
            if unknown_k[i][1] == "d":
                buses[index].angle = unknown_k1[i]
            elif unknown_k[i][1] == "v":
                buses[index].volts = unknown_k1[i]

        # Check for convergence
        if all(abs(power) < criterion for power in mismatch):
            converged = True
        else:
            converged = False

        # Check Q at the PV bus
        for bus in buses:
            Q_k1 = calc_q_at_bus(buses, yBusPolar, bus.index) / base_mva
            if bus.type == "PV" and abs(Q_k1) > abs(bus.Qcap) :
                print(f"PV bus, {bus.name}, changed to a PQ bus after iteration {k}")
                bus.type = "PQ"

        print(mismatch)
        # Message to the user if the system converged
        if converged:
            end_time = time.time()
            length = end_time - start_time
            print(f"System converged in {k} Iterations over {round(length, 3)} seconds")
            break

    # Message to the user in case the system didn't converge
    if not converged:
        end_time = time.time()
        length = end_time - start_time
        print(f"System did not converge after {k} Iterations over {round(length,3)} seconds")
    return unknown_k1