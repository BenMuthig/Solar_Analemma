# Calculates the position of the Sun and plots the solar analemma without
# empirical formulas
# Sources: https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Celestial_Mechanics_(Tatum)/
# 09%3A_The_Two_Body_Problem_in_Two_Dimensions/9.05%3A_Position_in_an_Elliptic_Orbit

# To do:
# 1. Clean up
# 2. Add code to individual functions
# 3. Bug fix: Sun altitude >90°/latitude <23.5°
# 4. Add functionality: longitude
# 5. Increase accuracy? (eg. stellar aberration, refraction, etc...)

import math
import numpy as np
import matplotlib.pyplot as plt

# Earth's semi-major axis in m
a = 149597887500
# Earth's orbital eccentricity
e = 0.0167086
# Earth's semi-minor axis in m
# b = 149576999826
b = a*math.sqrt(1-e**2)
# print(b)
# Universal gravitational constant in N*m*m/kg*kg
G = 6.67*10**-11
# Mass of the Sun in kg
M = 1.989*10**30
# Time at which perihelion occurred
# Most recent: 2023/1/4 1600h UTC
T0 = 0
# Time in Greenwich sidereal time: 2256h/82560s
GT0 = 82560
# number of steps
d = 366
# Time variable
# Vernal equinox 2023/3/20 2124h UTC
# Difference to perihelion: 6499441s
t = np.linspace(0+20*60*60, 365*24*60*60+20*60*60, d)
# t = np.linspace(6499440, 6499440, d)
# Regular seconds to sidereal seconds
t2s = 24*60*60/86164.0905
# Coordinates
x = np.linspace(0, 1, d)
y = np.linspace(0, 1, d)
# Eccentric anomaly v
v = np.linspace(0, 1, d)
# Right ascention
RA = np.linspace(0, 1, d)
# Declination
dec = np.linspace(0, 1, d)
# Azimuth angle
az = np.linspace(0, 1, d)
# Altitude/elevation angle
alt = np.linspace(0, 1, d)
# Angle in ecliptic of vernal equinox relative to perihelion: 1.32612508rad
ve = 1.32612508
# Earth axial tilt 2023/1/1
at = 23.4362871883*math.pi/180
# Observer latitude
lat = 65*math.pi/180

# Kepler's 3rd law for orbital period in s:
P = math.sqrt(4*(math.pi**2)*(a**3)/(G*M))
# print(P/(24*60*60))

# Function to calculate ecliptic coordinates given a time in seconds since perihelion
def ecliptic_coords(t):
    # Mean anomaly in radians
    MA = 2*math.pi*(t-T0)/P

    # Eccentric anomaly E in radians
    # The equation for E is transcendental, so we use Newton-Raphson
    # First guess:
    i = 0
    x0 = e*math.sin(MA[0])/(1-e*math.cos(MA[0]))
    E_old = MA[0]+x0*(1-0.5*x0**2)

    for ma in MA:
        error = 1

        while error > 0.0001:
            E = (ma-e*(E_old*math.cos(E_old)-math.sin(E_old)))/(1-e*math.cos(E_old))
            error = E-E_old
            E_old = E

        # Plot the ellipse traced by Earth's orbit
        x[i] = a*math.cos(E)
        y[i] = b*math.sin(E)

        # Angle zeroed at vernal equinox, and corrected for Earth vs Sun
        v[i] = (math.atan2((math.sqrt(1-e**2)*math.sin(E)), (math.cos(E)-e)) - ve - math.pi)
        RA[i] = math.atan2(math.cos(at)*math.sin(v[i]), math.cos(v[i])) + math.pi
        dec[i] = -math.asin(math.sin(at)*math.sin(v[i]))
        #print(i)
        #h = RA[i]*180/math.pi/15
        #min = (abs(h)-math.floor(abs(h)))*60
        #sec = (min-math.floor(min))*60
        #print(str(math.floor(abs(h))*(h/abs(h))) + "h " + str(math.floor(min)) + "min " + str(sec) + "s")
        #deg = dec[i]*180/math.pi
        #mins = (abs(deg)-math.floor(abs(deg)))*60
        #secs = (mins-math.floor(mins))*60
        #print(str(math.floor(abs(deg))*(deg/abs(deg))) + "° " + str(math.floor(mins)) + "' " + str(secs) + "''")

        # Local sidereal time (currently Greenwich) converted to radians
        LST = (((GT0 + t[i]*t2s) % (24*60*60)) / (60*60)) * 15 * math.pi/180
        H = LST - RA[i]
        print(i)
        #h = H*180/math.pi/15
        #min = (abs(h)-math.floor(abs(h)))*60
        #sec = (min-math.floor(min))*60
        #print(str(math.floor(abs(h))*(h/abs(h))) + "h " + str(math.floor(min)) + "min " + str(sec) + "s")
        alt[i] = math.asin(math.sin(dec[i])*math.sin(lat)+math.cos(dec[i])*math.cos(lat)*math.cos(H))
        az[i] = math.asin(-math.sin(H)*math.cos(dec[i])/math.cos(alt[i]))
        print(alt[i]*180/math.pi)
        i = i+1


ecliptic_coords(t)
#altt = math.asin(math.sin(22.5*math.pi/180)*math.sin(45*math.pi/180)+
#                 math.cos(22.5*math.pi/180)*math.cos(45*math.pi/180)*math.cos(60.9*math.pi/180))
#print(altt*180/math.pi)

#plt.scatter(x, y)
#plt.scatter(t, v)
#plt.scatter(t, dec)
#plt.scatter(RA, dec)
plt.plot(az, alt)
plt.axis('equal')
plt.show()
