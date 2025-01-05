import numpy as np
import matplotlib.pyplot as plt

# Suppose you stored iteration index and residual in each line, or just residual
# For instance, if you only store the residual in each line:
residuals = np.loadtxt("RESVEC.dat")
iterations = np.arange(len(residuals))

plt.figure()
plt.semilogy(iterations, residuals, label='Richardson Residual')
plt.xlabel('Iteration')
plt.ylabel('Residual (log scale)')
plt.legend()
plt.title('Convergence History (Richardson)')
plt.grid(True)
plt.show()

# # If you also have ERRVEC.dat:
# errors = np.loadtxt("ERRVEC.dat")
# plt.figure()
# plt.semilogy(iterations, errors, label='Richardson Error')
# plt.xlabel('Iteration')
# plt.ylabel('Error norm (log scale)')
# plt.legend()
# plt.title('Error History (Richardson)')
# plt.grid(True)
# plt.show()
