"""
EE 430 Power Analytical Methods of Power Systems - Fall 2025
Term Project - Newton-Raphson Algorithm
Joshua Consenz - 11/3/25

Defines the functions DAIG, OFF-DIAG, COL, and ROW for use in the main Newton-Raphson algorithm
"""
import numpy as np
import matplotlib.pyplot as plt

def diag(a):
    """
    Finds the diagonal elements of an input matrix, and returns a vector of the elements
    :param a: nxn matrix
    :return: Vector of diagonal elements
    """
    row, col = a.shape
    array = np.zeros(col, complex)

    for i in range(row):
        for j in range(col):
            #r = a[i, j].real
            #im = a[i, j].imag
            if i == j:
                array[i] = a[i,j]
    return array

def off_diag(a):
    """
    Finds the off-diagonal elements of an input matrix, and returns a vector of the elements
    :param a: nxn matrix
    :return: vector of off-diagonal elements
    """
    row, col = a.shape
    size = row*col-row
    array = np.zeros(size,complex)

    k = 0
    for i in range(row):
        for j in range(col):
            if i !=j:
                array[k] = a[i, j]
                k += 1

    return array

def col(a):
    """
    Finds each off-diagonal element in the array and returns a vector of the column indices
    :param a: an nxn matrix
    :return: vector of column indices
    """
    row, col = a.shape
    size = row*col - row
    array = np.zeros(size)
    k = 0
    for i in range(row):
        for j in range(col):
            if i != j and a[i,j] != 0:
                array[k] = j
                k += 1
    # Returns the array with all zeros removed from the array
    return array[array != 0]

def row(a):
    """
    Finds the off-diagonal, non-zero elements in each row of an nxn matrix, and returns a vector of the number of
    each of these elements in each row
    :param a: an nxn matrix
    :return: vector of the number of each non-zero, off-diagonal elements in each row
    """
    row, col = a.shape
    array = np.zeros(row)

    for i in range(row):
        for j in range(col):
            if (i != j) and (a[i, j] != 0):
                array[i] += 1
    return array