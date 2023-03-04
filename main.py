# Calculates the position of the Sun and plots the solar analemma without
# empirical formulas
# Sources: https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Celestial_Mechanics_(Tatum)/
# 09%3A_The_Two_Body_Problem_in_Two_Dimensions/9.05%3A_Position_in_an_Elliptic_Orbit#:~:text=r%3Dl1%2Becosv.
# &text=P2%3D4%CF%802,coordinates%20described%20by%20Equation%209.6.

import math
import numpy as np
import matplotlib.pyplot as plt

# Earth's semi-major axis in m
a = 149597887500
# Earth's semi-minor axis in m
b = 149576999826
# Earth's orbital eccentricity
e = 0.01671
# Universal gravitational constant in N*m*m/kg*kg
G = 6.67*10**-11
# Mass of the Sun in kg
M = 1.989*10**30
# Time at which perihelion occurred
T = 0
# number of steps
d = 1000000
# Time variable
t = np.linspace(0, 365.25*24*60*60, d)
# Coordinates
x = np.linspace(0, 1, d)
y = np.linspace(0, 1, d)

# Kepler's 3rd law for orbital period in s:
P = math.sqrt(4*(math.pi**2)*(a**3)/(G*M))
print(P)

# Mean anomaly in radians
MA = 2*math.pi*(t-T)/P

# Eccentric anomaly E in radians
# The equation for E is transcendental, so we use Newton-Raphson
# First guess:
i = 0
x0 = e*math.sin(MA[0])/(1-e*math.cos(MA[0]))
E_old = MA[0]+x0*(1-0.5*x0**2)

for ma in MA:
    error = 1

    while error > 0.0001:
        E = ma-e*(E_old*math.cos(E_old)-math.sin(E_old))/(1-e*math.cos(E_old))
        error = E-E_old
        E_old = E

    # Plot the ellipse traced by Earth's orbit
    x[i] = a*math.cos(E)
    y[i] = b*math.sin(E)
    i = i+1

plt.plot(x, y)
#plt.scatter(t, MA)
plt.show()
