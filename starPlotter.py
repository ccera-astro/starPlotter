from math import degrees, radians, pi, sqrt, atan, cos, sin, acos
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
import ephem
import coordinateConverter 
import time 
import os
import glob 

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--animation",type=float,default=1.0,help="Animation speed-up factor.")
    parser.add_argument("-p","--printLevel",type=int,help="Print level.")
    parser.add_argument("-t","--threshold",type=float,default=1.,help="S1400 threshold.")
    parser.add_argument("-o","--observatory",default='Carp',help="Observatory name: Jadwin, Belmar, Carp")
    parser.add_argument("--name",help="Name of object to intensify.")
    return parser.parse_args()

def getObserver(observatory) :
    obs=ephem.Observer()
    if observatory.lower() == 'belmar' :
        obs.lat='40.184693'
        obs.lon='-74.056116'
    elif observatory.lower() == 'jadwin' :
        obs.lat='40.34488'
        obs.lon='-74.65155'
    elif observatory.lower() == 'carp' :
        obs.lat='45.35025'
        obs.lon='-76.05609'
    else :
        print("ERROR in getObserver.  No observatory called {0:s}.  Exiting".format(observatory))
        exit() 
    return obs

def getSun(obs) :
    sun = ephem.Sun()
    sun.compute(obs)
    az, al = degrees(sun.az), degrees(sun.alt)
    az = -az + 90.
    return az, al

def getMoon(obs) :
    moon = ephem.Moon()
    moon.compute(obs)
    az, al = degrees(moon.az), degrees(moon.alt)
    az = -az + 90.
    return az, al

def getPolaris(obs) :
    polaris = ephem.readdb("Polaris,f|M|F7,2:31:48.704,89:15:50.72,2.02,2000")
    polaris.compute(obs)
    az, al = degrees(polaris.az), degrees(polaris.alt)
    az = -az + 90.
    return az, al

def getCasA(obs) :
    CasA = ephem.readdb("CasA,f|J|F7,23:23:25.8,58:48:00,2.02,2000")
    CasA.compute(obs)
    az, al = degrees(CasA.az), degrees(CasA.alt)
    #print("In getCasA: az = {0:.2f} al = {1:.2f}".format(az,al))
    az = -az + 90.
    return az, al

def getM31(obs) :
    M31 = ephem.readdb("M31,f|G,0:42:44,+41:16:8,4.16,2000,11433|3700|35")
    M31.compute(obs)
    az, al = degrees(M31.az), degrees(M31.alt)
    az = -az + 90.
    return az, al

def getM33(obs) :
    M33 = ephem.readdb("M33,f|G,1:33:50,+30:39:37,5.72,2000")
    M33.compute(obs)
    az, al = degrees(M33.az), degrees(M33.alt)
    az = -az + 90.
    return az, al

def getGalaxy(obs) :
    phi, theta = [], []
    for gLong in np.linspace(-180.,180,180) :
        (RA, dec)  = coordinateConverter.G2E(gLong, 0.)
        (az, al)  = coordinateConverter.E2H(obs, RA, dec)
        phi.append(radians(-az + 90.))
        theta.append(min(90.-al,90.))
    return phi, theta

class dummyTelescope() :
    def __init__(self) :
        self.az, self.alt = 10., 70.
        return 
    def position(self) :
        self.az += 1.0
        self.az = self.az % 360.
        self.alt += 1.0
        self.alt = self.alt % 90. 
        return self.az, self.alt
    
def getTelescope(obs,dummy) :
    import socket
    name = socket.gethostname().lower()
    carp_site = ('receiver' in name) or ('motion' in name)  
    if carp_site :
        import xmlrpc.client as xml 
        rpc = xml.ServerProxy("http://172.22.121.35:9090")
        values = rpc.query_both_axes()
        alt, az = values[0], values[1]
        return 90. - az, alt
    else :
        az, alt = dummy.position()
        return 90. - az, alt 
    
def telescopeInMotion(args,obs,dummy,lastTime,lastUVW) :
    #print("In telescopeInMotion(): lastTime={0:f} lastUVW={1:s}".format(lastTime,str(lastUVW)))
    try:
        phi, tht = getTelescope(obs,dummy) 
        #print "tht={0:f} radians".format(tht)
        cs, sn = cos(radians(tht)), sin(radians(tht)) 
        phi = radians(phi)
        UVW = [sn*cos(phi), sn*sin(phi), cs]
        angle = degrees(acos(lastUVW[0]*UVW[0]  + lastUVW[1]*UVW[1] + lastUVW[2]*UVW[2]))
        if angle > 1. :    # dish is moving
            return True, time.time(), UVW
        else :
            return False, time.time(), UVW
    except:
        print("TelescopeInMotion failed.")

    return False, lastTime, lastUVW

    
def getPulsarLists(obs,args) :
    #;PSRJ;RAJ;DECJ;P0;S1400;
    #1;J0332+5434;03:32:59.3;+54:34:43.5;0.714520;203.00;
    pulsar, pulsarP0, pulsarS1400 = {}, {}, {}
    for line in open('pulsarList.txt','r').readlines()[1:] :
        if len(line) < 20 : continue 
        psr = ephem.FixedBody()
        #psr._epoch="1950/1/1 00:00:00"
        vals = line.split(';')
        name, psr._ra, psr._dec, P0, S1400 = vals[1], ephem.hours(vals[2]), ephem.degrees(vals[3]), float(vals[4]), float(vals[5])
        (az, al)  = coordinateConverter.E2H(obs, degrees(psr._ra), degrees(psr._dec))
        if (al > 0.) :
            pulsar[name] = psr
            pulsarP0[name] = P0
            pulsarS1400[name] = S1400

    return pulsar, pulsarP0, pulsarS1400

def printPulsarList(pulsar, pulsarP0, pulsarS1400) :
    count = 0
    for psr in pulsar.keys() :
        if pulsarS1400[psr] < args.threshold : continue
        count += 1
        print("Psr: {0:11s} RA={1:6.2f} dec={2:6.2f} P0={3:9.6f} S1400={4:6.2f}".format(
            psr,degrees(pulsar[psr]._ra),degrees(pulsar[psr]._dec), pulsarP0[psr], pulsarS1400[psr])) 
        
    print("{0:d} pulsars in total".format(count))


# begin execution here 
args = getArgs()
obs = getObserver(args.observatory)
dummy = dummyTelescope() 
obs.date = ephem.now() 

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},figsize=(11, 7), facecolor='#0D5973')

ax.set_position([0.0, 0.1, 0.7, 0.8])
ax.set_rlabel_position(-25.)  # Move radial labels away from plotted line
ax.grid(True,color='white')
ax.set_facecolor('black')
ax.set_xticks(ax.get_xticks())
ax.xaxis.set_ticklabels(['E','','N','','W','','S',' '],color='white') 
ax.set_yticks([10.,20.,30.,40.,50.,60.,70.,80.,90.])
ax.yaxis.set_ticklabels(['','','60','','','30','','','0'],color='white')

# add the telescope direction
az, al = getTelescope(obs,dummy) 
telescopeMarker1, = ax.plot([radians(az)], [90.-al], '+', color = 'red', markersize=20, label='Telescope' )
telescopeMarker2, = ax.plot([radians(az)], [90.-al], 'o', markerfacecolor = 'none', markeredgecolor = 'red', markersize=8 )

az, al = getSun(obs)
#sunMarker1, = ax.plot(radians(az), 90.-al, 'y*', markersize=20.)
sunMarker2, = ax.plot(radians(az), 90.-al, 'o', color="#ffff00", markersize=18.0, label="Sun")

az, al = getMoon(obs) 
moonMarker, = ax.plot(radians(az), 90.-al, 'wo', markersize=18., label="Moon")

az, al = getPolaris(obs) 
polarisMarker, = ax.plot(radians(az), 90.-al, 'w*', markersize=10., label="Polaris") 

az, al = getCasA(obs) 
casAMarker, = ax.plot(radians(az), 90.-al, 'r*', markersize=10., label="Cas A")

az, al = getM31(obs) 
M31Marker, = ax.plot(radians(az), 90.-al, 'wd', markersize=10., label="M31") 

az, al = getM33(obs) 
M33Marker, = ax.plot(radians(az), 90.-al, 'gd', markersize=10., label="M33")  

#draw the galaxy
phi, theta = getGalaxy(obs) 
galaxyLine, = ax.plot(phi,theta, 'g-') 

#galactic longitude markers
galacticLongitudeMarkers, galacticLongitudeText = {}, {}
for gLong in np.linspace(-160.,180.,18) :
    (RA, dec)  = coordinateConverter.G2E(gLong, 0.)
    (az, al)  = coordinateConverter.E2H(obs, RA, dec)
    if (al > 0. ):
        gm, = ax.plot(radians(-az+90.), 90.-al, 'g+', markersize=10.)
        galacticLongitudeMarkers[gLong] = gm 
        txt = "  {0:.0f}".format(gLong)
        gt = ax.text(radians(-az+90.), 90.-al, txt, color='g')
        galacticLongitudeText[gLong] = gt

# dummy point to include zero 
ax.plot(90.,0.,'w.')
   
# create a box onto which the pulsar list can be printed 
ax2 = fig.add_subplot(2,9,16,facecolor='#0D5973')
ax2.get_xaxis().set_visible(False) 
ax2.get_yaxis().set_visible(False)
ax2.axis('off')

pulsar, pulsarP0, pulsarS1400 = getPulsarLists(obs,args)

# make a list of the top five pulsars
count, N = 0, 5
txt = ['A','B','C','D','E']
sortedPulsarList = sorted(pulsarS1400, key=pulsarS1400.get, reverse=True)
topPulsarList = sortedPulsarList[0:N]

#add the pulsars
nObject = len(pulsar)
pulsarMarkers, pulsarLabels = {}, {} 
for p in pulsar.keys() :
    psr = pulsar[p]
    S1400 = pulsarS1400[p]
    if S1400 < args.threshold : continue
    (az, al)  = coordinateConverter.E2H(obs, degrees(psr._ra), degrees(psr._dec))
    az = -az + 90.
    pm, = ax.plot(radians(az), 90.-al, 'o', color = '#90F5F5', markersize=max(sqrt(pulsarS1400[p]),1.) )
    pulsarMarkers[p] = pm  
    if p in topPulsarList :
        ii = topPulsarList.index(p) 
        (az, al)  = coordinateConverter.E2H(obs, degrees(psr._ra), degrees(psr._dec))
        az = -az + 90.
        pl = ax.text(radians(az), 45., txt[ii], color = '#90F5F5', size=16)
        pulsarLabels[p] = pl 

# time text
yPos, dy = 0.9, 0.12 
timeText = ax2.text(.1, yPos , args.observatory + ':  ' + str(obs.date) + ' UTC', color='white', fontsize=14)
yPos -= 0.05
for i in range(N) :
    yPos -= dy
    psr = sortedPulsarList[i]
    st1 = "{0:s}: {1:s}  {2:6.1f} mJy  {3:5.1f} ms".format(txt[i],psr,pulsarS1400[psr],1000.*pulsarP0[psr])
    ax2.text(-0.25, yPos, st1, color='#90F5F5', fontsize=14)
    
angle = radians(40.)
ax.legend(loc="upper left",fontsize=16, labelcolor='white', facecolor='#0D5973', bbox_to_anchor=(.6 + np.cos(angle)/1.5, .5 + np.sin(angle)/1.5))

plt.show(block=False)

# animate the display 
lastTime, lastUVW = 0., [0., 0., -1.]
iter = 0
os.system("echo GO > starGO")
while len(glob.glob('starGO')) > 0 :
    firstIteration = True
    for i in range(12) :
        if args.animation > 1.01 : time.sleep(1.)
        #else : time.sleep(10.0)
            
        inMotion, lastTime, lastUVW = telescopeInMotion(args,obs,dummy,lastTime,lastUVW)
        #inMotion = False 
        if firstIteration or inMotion :
            firstIteration = False 
            az, al = getTelescope(obs,dummy)
            # drop a bread crumb to keep track of where the telescope has been
            telescopeMarker, = ax.plot([radians(az)], [90.-al], 'ro', markersize=4 )
            ax.set_rmax(90.0)
            break

    obs = getObserver(args.observatory)
    if args.animation > 1. :
        obs.date = ephem.now() + args.animation*1.5*iter/86160.
        iter += 1

    timeText.set_text(args.observatory + ':  ' + str(obs.date) + ' UTC')
    
    az, al = getTelescope(obs,dummy)    
    telescopeMarker1.set_data([radians(az)],[90.-al])
    telescopeMarker2.set_data([radians(az)],[90.-al])

    az, al = getSun(obs)
    #sunMarker1.set_data([radians(az)],[90.-al])
    sunMarker2.set_data([radians(az)],[90.-al]) 

    az, al = getMoon(obs) 
    moonMarker.set_data([radians(az)],[90.-al])

    az, al = getPolaris(obs) 
    polarisMarker.set_data([radians(az)],[90.-al])

    az, al = getCasA(obs) 
    casAMarker.set_data([radians(az)],[90.-al])

    az, al = getM31(obs) 
    M31Marker.set_data([radians(az)], [90.-al]) 

    phi, theta = getGalaxy(obs) 
    galaxyLine.set_data(phi,theta) 

    #galactic longitude markers
    for gLong in galacticLongitudeMarkers.keys() :
        (RA, dec)  = coordinateConverter.G2E(gLong, 0.)
        (az, al)  = coordinateConverter.E2H(obs, RA, dec)
        if al < 0. : al = -1000.
        galacticLongitudeMarkers[gLong].set_data(radians(-az+90.), 90.-al)
        galacticLongitudeText[gLong].set_position([radians(-az+90.), 90.-al]) 

    #pulsars
    for p in pulsarMarkers.keys() :
        psr = pulsar[p]
        (az, al)  = coordinateConverter.E2H(obs, degrees(psr._ra), degrees(psr._dec))
        az = -az + 90.
        pulsarMarkers[p].set_data(radians(az), 90.-al)
        if p in topPulsarList :
            if al < 0. : al = -1000.
            pulsarLabels[p].set_position([radians(az), 95.-al])
        
    plt.pause(5.0)












