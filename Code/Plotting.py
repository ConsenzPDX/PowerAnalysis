import numpy as np
import matplotlib.pyplot as plt

criterion = np.array([1.00E-02, 1.00E-03, 1.00E-04, 1.00E-05, 1.00E-06, 1.00E-07, 1.00E-08, 1.00E-09, 1.00E-10, 1.00E-11, 1.00E-12, 1.00E-13, 1.00E-14])
iterations = np.array([6, 9, 12, 14, 17, 19, 22, 24, 28, 29, 33, 34, 38])

plt.figure(figsize=(8,5))
plt.plot(criterion, iterations, marker='o')
plt.title('Newton-Raphson Iterations vs Convergence Criterion')
plt.xlabel('Criterion')
plt.ylabel('Number of Iterations')
plt.xscale('log')
plt.ylim(0, 40)
plt.grid(True)
plt.show()