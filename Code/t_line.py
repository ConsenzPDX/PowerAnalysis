import string
import numpy as np
from bus import Bus

class T_line:
    def __init__(self, start:Bus, end:Bus, Rse: float, Xse:float, Gsh: float, Bsh: float, Rating:int):
        """
        Data type describing a transmission line with a start(from) and end(t0) bus, as well as impedance, admittance,
        and power rating
        :param start:
        :param end:
        :param Rse:
        :param Xse:
        :param Gsh:
        :param Bsh:
        :param Rating: Power (MVA) rating of the transmission line
        """
        self.start = start
        self.end = end
        self.Rse = Rse
        self.Xse = Xse
        self.Gsh = Gsh
        self.Bsh = Bsh
        self.Rating = Rating
