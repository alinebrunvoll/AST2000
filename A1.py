import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt


def P(a, b, my, sigma):
    f = lambda x: (1 /(np.sqrt(2*np.pi)*sigma)) * np.exp(-0.5 *((x-my)/sigma)**2)
    return integrate.quad(f, a, b)[0]


P(-1, 1, 0, 1)
P(-2, 2, 0, 1)
P(-3, 3, 0, 1)

# assert FWHM = 2 * sigma * np.sqrt(2*log(2))


#  the Maxwell-Boltzmann Distribution
N = 10**5
sigma = 1; my = 0
f = lambda x: (1 /(np.sqrt(2*np.pi)*sigma)) * np.exp(-0.5 *((x-my)/sigma)**2)

vx = np.linspace(-2.5, 2.5, N)
plt.plot(vx, f(vx))
plt.title('Sannsynlighetstetthet til fart i x-retning (10^4 m/s)')
plt.xlabel('x')
plt.ylabel('Sannsynlighet')
plt.show()

# vx = np.linspace(5, 30)
print(P(5, 30, 0, 1)* N)
