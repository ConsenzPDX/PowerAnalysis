import string
import numpy as np
import matplotlib.pyplot as plt
from Functions import diag
from Functions import off_diag
from Functions import col
from Functions import row
from bus import Bus
from t_line import T_line

Alan = Bus("SL", 0.98, 0, 0, 0, 0, 0)
Betty = Bus("PV", 1.00, 210, 50, 0, 0, 0)
Clyde = Bus("PQ", 1.00, 0, 0, 110, 85, 150)
Doug = Bus("PQ", 1.00, 0, 0, 100, 95, 50)
Eve = Bus("PQ", 1.00, 0, 0, 150, 120, 0)

AB = T_line(Alan, Betty, 0.009, 0.041, 0.000, 0.000, 125)
BE = T_line(Betty, Eve, 0.006, 0.037, 0.000, 0.000, 250)
AD = T_line(Alan, Doug, 0.007, 0.055, 0.000, 0.000, 200)
DE = T_line(Doug, Eve, 0.006, 0.045, 0.000, 0.000, 125)
DC = T_line(Doug, Clyde, 0.011, 0.061, 0.000, 0.000, 80)
CE = T_line(Clyde, Eve, 0.010, 0.051, 0.000, 0.000, 75)