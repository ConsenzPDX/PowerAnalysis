"""
EE 430 Power Analytical Methods of Power Systems - Fall 2025
Term Project - Newton-Raphson Algorithm
Joshua Consenz - 11/10/25

Class for creating transmission lines in the main program newton_raphson.py
Currently is only comprised of a constructor
"""
import string
import numpy as np

class Bus:
    def __init__(self, name: str, type: str, volts: float, Pgen: float, Qgen: float, Pload:float, Qload:float, Qcap:float):
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

        self.name = name
        self.type = type
        self.volts = volts
        self.Pgen = Pgen
        self.Qgen = Qgen
        self.Pload = Pload
        self.Qload = Qload
        self.Qcap = Qcap

    def __setIndex__(self, index):
        self.index = index

    def __setVoltage__(self, voltage):
        self.volts = voltage

    def __voltAngle__(self, angle):
        self.angle = angle

    def __netP__(self):
        self.netP = self.Pgen - self.Pload

    def __netQ__(self):
        self.netQ = self.Qgen - self.Qload