"""
EE 430 Power Analytical Methods of Power Systems - Fall 2025
Term Project - Newton-Raphson Algorithm
Joshua Consenz - 11/3/25

Class for creating transmission lines in the main program Power Analysis.py
Currently is only comprised of a constructor
"""
import string
import numpy as np

class Bus:
    def __init__(self, type: string, volts: float, Pgen: float, Qgen: float, Pload:float, Qload:float, Qcap:float, index:int):
        """
        Constructor for the Bus class. Creates a bus with a type, voltage, power flow values, and reactive power cap
        :param type: SL(slack), PV(gen), or PQ(load) are the used values. This is mainly for me to use, not the program
        :param volts: Per unit voltage of the bus
        :param Pgen: Real power generated at the bus
        :param Qgen: Reactive power generated at the bus
        :param Pload: Real power consumed by the bus
        :param Qload: Reactive power consumed by the bus
        :param Qcap: Cap of the reactive power the bus can handle
        """

        self.type = type
        self.volts = volts
        self.Pgen = Pgen
        self.Qgen = Qgen
        self.Pload = Pload
        self.Qload = Qload
        self.Qcap = Qcap
        self.index = index
