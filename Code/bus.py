import string
import numpy as np

class Bus:
    def __init__(self, type: string, volts: float, Pgen: float, Qgen: float, Pload:float, Qload:float, Qcap:float):
        """

        :param type:
        :param volts:
        :param Pgen:
        :param Qgen:
        :param Pload:
        :param Qload:
        :param Qcap:
        """

        self.type = type
        self.volts = volts
        self.Pgen = Pgen
        self.Qgen = Qgen
        self.Pload = Pload
        self.Qload = Qload
        self.Qcap = Qcap