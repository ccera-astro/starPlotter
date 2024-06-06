# Usage coordinateConverter.py  H  Az Al      
#    or coordinateConverter.py  E  RA Dec
#    or coordinateConverter.py  G  Lon Lat
# where H, E, G is the input coordinate system (Horizontal, Equatorial, Galactic)
#   and the numerical inputs are in degrees
# obs is a pyephem observer object

from math import degrees, radians, sin, cos, asin, atan2
import ephem
import time
from datetime import datetime, timedelta
import sys

def H2E(obs, dAz, dAl) :
    (RA, dec) = obs.radec_of(str(dAz),str(dAl))
    RA_deg = degrees(RA)
    if RA_deg > 180. : RA_deg -= 360.
    return RA_deg, degrees(float(dec))

def E2G( RA, dec) :
    RA_hours = RA/15.
    eCoords = ephem.Equatorial(str(RA_hours),str(dec))
    gCoords = ephem.Galactic(eCoords)
    return degrees(gCoords.lon), degrees(gCoords.lat)

def E2H(obs,  RA, dec) :
    RA_hours = RA/15.
    eCoords = ephem.Equatorial(str(RA_hours),str(dec))
    fakeStar = ephem.FixedBody()
    fakeStar._ra = eCoords.ra
    fakeStar._dec = eCoords.dec
    fakeStar.compute(obs)
    return degrees(float(fakeStar.az)), degrees(float(fakeStar.alt))

def G2E(gLon, gLat) :
    gCoords = ephem.Galactic(str(gLon),str(gLat))
    eCoords = ephem.Equatorial(gCoords)
    RA = degrees(float(eCoords.ra))
    if RA > 180. : RA -= 360.
    return RA , degrees(eCoords.dec)

def getObserver() :
    princeton=ephem.Observer()
    princeton.lat='40.344892'
    princeton.lon='-74.651692'
    return princeton

def coordinateConverterFunction(coordinateSystem, phi, theta, offsetHours) :
    obs = getObserver()
    #obs.date = datetime.now() + timedelta(hours=4)
    obs.date = datetime.utcnow() + timedelta(hours=offsetHours)
    if coordinateSystem == 'H' :
        (az, alt) = (phi, theta)
        (RA, dec) = H2E(obs, az, alt)
        (gLon, gLat) = E2G(RA, dec)
    elif coordinateSystem == 'E' :
        (RA, dec) = (phi, theta) 
        (az, alt) = E2H(obs, RA, dec)
        (gLon, gLat) = E2G(RA, dec)
    elif coordinateSystem == 'G' :
        (gLon, gLat) = (phi, theta) 
        (RA, dec) = G2E(gLon, gLat) 
        (az, alt) = E2H(obs, RA, dec)
    else :
        print("ERROR in coordinateConverter.py invalid coordinate system = " + coordinateSystem)
        return 

    print("Coordinate Summary for " + str(datetime.now()+timedelta(hours=offsetHours)) + " Eastern Time") 
    print("Horizontal:  Az ={0:8.2f} Alt ={1:8.2f}".format(az,alt))
    print("Equatorial:  RA ={0:8.2f} dec ={1:8.2f}".format(RA,dec))
    print("  Galactic: lon ={0:8.2f} lat ={1:8.2f}".format(gLon,gLat))
    return az, alt, RA, dec, gLon, gLat

if __name__ == '__main__':
    # parse command-line arguments
    arglen = len(sys.argv)
    if (arglen < 4) :
        print("ERROR in coordinateConverter.py invalid number of arguments = " + arglen + " (should be 4).")
    else :
        offsetHours = 0
        if arglen >= 5 : offsetHours = float(sys.argv[4])
        coordinateConverterFunction(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]),offsetHours)
                                 
