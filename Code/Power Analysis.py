import numpy as np
import matplotlib.pyplot as plt
from Functions import diag
from Functions import off_diag
from Functions import col
from Functions import row

class bus():
    def __init__(self, type, v, d, P, Q):
        """
        Constructor for the bus class
        :param type: string input for the bus type: slack, gen, load
        :param v: voltage at the bus
        :param d: voltage angle at the bus
        :param P: net real power of the bus
        :param Q: net reactive power of the bus
        """
        self.type = type
        self.v = v
        self.d = d
        self.P = P
        self.Q = Q