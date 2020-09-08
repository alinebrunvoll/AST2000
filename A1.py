# EGEN KODE
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt


def P(a, b, my, sigma):
    f = lambda x: (1 /(np.sqrt(2*np.pi)*sigma)) * np.exp(-0.5 *((x-my)/sigma)**2)
    return integrate.quad(f, a, b)[0]  # Integrerer sannsynlighetstettheten fra a til b

# Printer ut sannsynligheten for at et tiilfeldig tall er innenfor standardavviket sigma
print(P(-1, 1, 0, 1))
print(P(-2, 2, 0, 1))
print(P(-3, 3, 0, 1))



#  the Maxwell-Boltzmann Distribution
N = 10**5
T = 3000
k = 1.38 * 10**(-23)                                 # Boltzmann-konstanten
m = 2*1.66*10**(-27)

sigma = np.sqrt((k*T) /m); my = 0
f = lambda x: (1 /(np.sqrt(2*np.pi)*sigma)) * np.exp(-0.5 *((x-my)/sigma)**2)

vx = np.linspace(-2.5*10**(4), 2.5*10**(4), N)
plt.plot(vx, f(vx))
plt.title('Sannsynlighetstetthet til hastighet i x-retning (10^4 m/s)')
plt.xlabel('Hastighet')
plt.ylabel('Sannsynlighet')
plt.show()


print(P(5, 30, 0, 1) * N)
# antalltettheteten ^^


v = np.linspace(0, 3*10**4, N)
plt.plot(v, np.sqrt(f(v)**2))
plt.title('Absolutt hastighetsfordeling')
plt.xlabel('Hastighet'); plt.ylabel('Sannsynlighet')
plt.show()
# Dette er ikke i konflikt, da vi kun ser på verdier fra 0 til tre
# da v kun kan være positiv

def P_(m, k, T, a, b):
    p = lambda v: v*(m/(2*np.pi*k*T))**(3/2) * np.exp(-(1/2)*(m*v**2/(k*T))) *4*np.pi*v**2
    return integrate.quad(p, a, b)[0]

k = 1.38 * 10**(-23)
print(P_(1, k, 3000, 0, np.inf))
