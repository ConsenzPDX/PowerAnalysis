import numpy as np
import matplotlib.pyplot as plt
from Functions import diag
from Functions import off_diag
from Functions import col
from Functions import row

A = np.array([[0,2,3],[0,5,6],[0,0,9]])
print(A.ndim)
print(A.shape)
B = diag(A)
print(B)
print(off_diag(A))
print(col(A))
print(row(A))