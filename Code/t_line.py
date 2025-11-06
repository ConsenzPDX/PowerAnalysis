"""
EE 430 Power Analytical Methods of Power Systems - Fall 2025
Term Project - Newton-Raphson Algorithm
Joshua Consenz - 11/6/25

Class for creating buses in the main program Power Analysis.py
Currently is only comprised of a constructor
"""
import string
import numpy as np
from bus import Bus

class T_line:
    def __init__(self, start:Bus, end:Bus, Rse: float, Xse:float, Gsh: float, Bsh: float, Rating:int):
        """
        Data type describing a transmission line with a start(from) and end(t0) bus, as well as impedance, admittance,
        and power rating
        :param start: starting(from) bus of the transmission line
        :param end: ending(t) bus of the transmission line
        :param Rse: Resistance
        :param Xse: Reactance
        :param Gsh: Conductance
        :param Bsh: Susceptance
        :param Rating: Power (MVA) rating of the transmission line
        """

        # Assigns the values from the constructor to the new T_line's values
        self.start = start
        self.end = end
        self.Rse = Rse
        self.Xse = Xse
        self.Gsh = Gsh
        self.Bsh = Bsh
        self.Rating = Rating

        self.impedance = complex(self.Rse, self.Xse)
        self.admittance = 1/self.impedance
        self.Gsh = self.admittance.real
        self.Bsh = self.admittance.imag