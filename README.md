# Solar_Analemma
Calculates the solar analemma based on latitude.
Iterates through a time array and calculates ecliptic, equatorial, and horizontal coordinates.
Earth's position relative to the Sun is calculated using Kepler's laws rather than using an empirical formula to describe the Sun's position.

To do:
- Currently all calculations are based on the prime meridian. Functionality to specify longitude could be added.
- Various additional calculations could be done to increase accuracy slightly (factoring refraction, stellar aberration, etc...)

Sources for calculation:
- Earth's position relative to Sun in ecliptic: https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Celestial_Mechanics_(Tatum)/09%3A_The_Two_Body_Problem_in_Two_Dimensions/9.05%3A_Position_in_an_Elliptic_Orbit
- Converting from ecliptic to equatorial coordinates: https://en.wikipedia.org/wiki/Position_of_the_Sun
- Converting from equatorial to horizontal coordinates: https://sceweb.sce.uhcl.edu/helm/WEB-Positional%20Astronomy/Tutorial/Conversion/Conversion.html

Sources for hardcoding initial values and testing/debugging:
- Sidereal time: https://aa.usno.navy.mil/data/siderealtime
- Testing RA, declination, and altitude angle: https://theskylive.com/sun-info