# Calculates the position of the Sun and plots the solar analemma without
# empirical formulas
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
# Universal gravitational constant in N*m*m/kg*kg
G = 6.67*10**-11
# Mass of the Sun in kg
M = 1.989*10**30
# Time at which perihelion occurred
# Most recent: 2023/1/4 1600h UTC
# Time in Greenwich sidereal time: 2256h/82560s
ST0 = 82560
# number of steps
d = 366
# Time variable
# Vernal equinox 2023/3/20 2124h UTC
# Difference to perihelion: 6499441s
T = np.linspace(0+20*60*60, 365*24*60*60+20*60*60, d)
# Regular seconds to sidereal seconds
t2s = 24*60*60/86164.0905
# Coordinates
x = np.linspace(0, 1, d)
y = np.linspace(0, 1, d)
# Eccentric anomaly v
v = np.linspace(0, 1, d)
# Right ascension
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


def rad2deg(rad):
    return rad*180/math.pi


# Uses abs in floor to ensure negative numbers are handled correctly
def deg2hms(deg):
    h = math.floor(abs(deg))*(deg/abs(deg))
    min = (abs(deg)-math.floor(abs(deg)))*60
    sec = (min-math.floor(min))*60
    min = math.floor(min)
    return h, min, sec


# Calculates ecliptic coordinates given a time in seconds since perihelion
def ecliptic_coords(t):
    # Mean anomaly in radians
    ma = 2*math.pi*t/P

    # Eccentric anomaly E in radians
    # The equation for E is transcendental, so we use Newton-Raphson
    # First guess:
    x0 = e*math.sin(ma)/(1-e*math.cos(ma))
    E_old = ma+x0*(1-0.5*x0**2)
    error = 1
    while error > 0.0001:
        E = (ma - e * (E_old * math.cos(E_old) - math.sin(E_old))) / (1 - e * math.cos(E_old))
        error = E - E_old
        E_old = E

    # Rectangular coordinates of ellipse traced by Earth's orbit
    x = a * math.cos(E)
    y = b * math.sin(E)

    # True anomaly zeroed at vernal equinox, and corrected for Earth vs Sun
    v = (math.atan2((math.sqrt(1 - e ** 2) * math.sin(E)), (math.cos(E) - e)) - ve - math.pi)

    return x, y, v


# Calculates Sun's position in equatorial coordinates using true anomaly
def equatorial_coords(v):
    RA = math.atan2(math.cos(at) * math.sin(v), math.cos(v)) + math.pi
    dec = -math.asin(math.sin(at) * math.sin(v))
    return RA, dec


# Calculates local sidereal time (currently fixed at prime meridian)
# Arguments:
# ST0: sidereal time 0
# t: time since ST0
# lon: longitude
def lst_calc(ST0, t):
    return (((ST0 + t * t2s) % (24 * 60 * 60)) / (60 * 60)) * 15 * math.pi / 180


# Calculates horizontal coordinates using equatorial coordinates, observers latitude, and local sidereal time
def horizontal_coords(RA, dec, lat, LST):
    H = LST - RA
    alt = math.asin(math.sin(dec) * math.sin(lat) + math.cos(dec) * math.cos(lat) * math.cos(H))
    az = math.asin(-math.sin(H) * math.cos(dec) / math.cos(alt))
    return az, alt


i = 0
for t in T:
    x, y, v = ecliptic_coords(t)
    RA, dec = equatorial_coords(v)
    LST = lst_calc(ST0, t)
    az[i], alt[i] = horizontal_coords(RA, dec, lat, LST)
    i = i+1


#plt.scatter(x, y)
#plt.scatter(t, v)
#plt.scatter(t, dec)
#plt.scatter(RA, dec)
plt.plot(az, alt)
plt.axis('equal')
plt.show()
