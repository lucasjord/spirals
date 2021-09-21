#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import struct, os, math, random, matplotlib, datetime, argparse
from matplotlib import rc, font_manager
from numpy import pi, arcsin, cos, sin, percentile, exp

from mpl_toolkits.basemap import Basemap

from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u

import matplotlib.animation as animation


matplotlib.rcParams.update({'font.size': 14})
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)

##############################################################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose",
                        help="increase output verbosity",
                        action="count",default=0)
    parser.add_argument('tec_map',
                        help='jplg TEC map e.g. jplgDDD0.YYi',
                        type=str)
    parser.add_argument('-m','--maser',
                        help='Maser number or name e.g. G232.62 or s1',default=None)
    args = parser.parse_args()
    '''
    Load in specified TEC map. Only can load in one at a time atm.
    '''
    t, long, lat, tec, rms = read_jpltec(args.tec_map)
    '''
    Get maser of interest. Defaults to G232.62 aka s1 at the moment.
    '''
    source  = SkyCoord(ra='07:32:09.79',dec='-16:58:12.4',unit=(u.hourangle,u.deg)) #source='G232.62+0.99',
    gst     = 15*t.sidereal_time('apparent').value
    
    # deletes old animations if they exist.
    try:
        animate1.event_source.stop()
    except NameError:
        ''
    '''
    Makes basemap of Australia-NZ area. Currently hard-coded.
    '''
    x, y    = np.meshgrid(long, lat)
    fig, ax = plt.subplots(1,figsize=(14,7))
    ax = Basemap(llcrnrlon=70.,llcrnrlat=-61.,urcrnrlon=220.,urcrnrlat=15.,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='l',projection='merc',ax=ax)
    ax.drawcoastlines()
    ax.drawstates()
    ax.drawparallels(np.arange(-90,90,10),labels=[1,1,0,1]);
    ax.drawmeridians(np.arange(-180,180,15),labels=[1,1,0,1]);
    '''
    Put antenna locations on map. Currently hard-coded.
    '''
    a_lat = [-42.804,-31.868,-29.047,-14.375,-36.43]
    a_lon = [147.440,133.809,115.350,132.150,174.66]
    col = ['g','orange','b','r','y']
    for i in range(len(a_lat)):
        ax.scatter(a_lon[i],a_lat[i],latlon=True,color=col[i],zorder=99)
    '''
    Put on axis objects- 0:TEC colormap, 1:source position, 2:date/time, 
    3:antenna-source LoS (W.I.P) and colorbar (fixed -2 to 40)
    '''
    m = 0
    plot0 = [ax.pcolor(x,y,tec[:,:,m], latlon=True,cmap='magma',vmin=-2,vmax=40)]
    plot1 = [ax.scatter([source.ra.deg-gst[m]],[source.dec.deg],latlon=True,color='w',s=50,zorder=99)]
    plot2 = [fig.text(0.175,0.15,s=f'{t.iso[m]}',color='w',size=12)]
    plot3 = []
    #for i in range(len(a_lat)):
    #    plot3.append(ax.drawgreatcircle(a_lon[i],a_lat[i],
    #    source.ra.deg-gst[m],source.dec.deg,color='w',ls='-.',lw=0.25))
    c = plot0[0]
    ax.colorbar(c,pad='8%',label=R'$I_e$ (TECU)')
    '''
    Create animation. 
    '''
    animate1 = animation.FuncAnimation(fig, update, range(tec.shape[2]), 
                                   fargs=(t,gst,source,x,y,tec,plot0,plot1,plot2,plot3,fig,ax),interval=500)
    plt.show()
##############################################################################################

def get_file(path):
    #opens and external file and makes it into a list
    fopen = path
    f=open(fopen, 'r+')
    g=list(f)
    g=map(lambda s: s.strip(), g)
    return np.array(list(g))

def splitt(old_list):
    #splits the list entries into sublists
    new_list=[]
    for i in old_list:
        new_list+=[i.split()]
    return new_list

def update(m,t,gst,source,x,y,tec,plot0,plot1,plot2,plot3,fig,ax):
    plot0[0].remove()
    plot0[0] = ax.pcolor(x,y,tec[:,:,m], latlon=True,cmap='magma',vmin=-2,vmax=40)
    plot1[0].remove()
    plot1[0] = ax.scatter([source.ra.deg-gst[m]],[source.dec.deg],latlon=True,color='w',s=50,zorder=99)
    plot2[0].remove()
    plot2[0] = fig.text(0.175,0.15,s=f'{t.iso[m]}',color='w',size=12)
    #for ln in plot3: ln[0].remove()
    #for i in range(len(a_lat)):
    #    plot3.append(ax.drawgreatcircle(a_lon[i],a_lat[i],
    #    source.ra.deg-gst[m],source.dec.deg,color='w',ls='-.',lw=0.25))

def read_jpltec(tec_loc):
    tecfile  = get_file(tec_loc)
    start    = [s for s in tecfile if 'HEADER' in s][0]
    tecends  = [s for s in tecfile if 'END OF TEC MAP' in s]
    rmsends  = [s for s in tecfile if 'END OF RMS MAP' in s]
    header   = tecfile[:np.where(tecfile==start)[0][0]-1]
    data     = tecfile[ np.where(tecfile==start)[0][0]+1:np.where(tecfile==tecends[-1])[0][0]]
    rdat     = tecfile[ np.where(tecfile==tecends[-1])[0][0]+1:np.where(tecfile==rmsends[-1])[0][0]]

    # read tec map
    tecstart = np.array([np.where(data==s)[0][0] for s in data if 'START OF TEC MAP' in s])
    lat = []
    for m in range(len(tecstart)):
        tmap = data[tecstart[m]+1:tecstart[m]+428]

    # read rms map
    rmsstart = np.array([np.where(rdat==s)[0][0] for s in rdat if 'START OF RMS MAP' in s])
    for m in range(len(rmsstart)):
        rmap = data[rmsstart[m]+1:rmsstart[m]+428]

    t, lat, long = [], np.arange(87.5,-87.6,-2.5), np.arange(-180,180.1,5)
    tec  = np.zeros(shape=(len(lat),len(long),12))
    rms  = np.zeros(shape=(len(lat),len(long),12))
    for m in range(12):
        # format tec map
        tmap        = data[tecstart[m]+1:tecstart[m]+428]
        time_string = tmap[0].split()
        map_t       = datetime.datetime(year=int(time_string[0]),month=int(time_string[1]),day=int(time_string[2]),
                     hour=int(time_string[3]),minute=int(time_string[4]),second=int(time_string[5]))
        latstart = np.array([np.where(tmap==s)[0][0] for s in tmap if 'LAT/LON1/LON2/DLON/H' in s])
        for n in range(len(lat)):
            tec[n,:,m] = np.array([s for d in splitt(tmap[latstart[n]+1:latstart[n]+6]) for s in d],dtype='float')
        ####################################################
        rmap = rdat[rmsstart[m]+1:rmsstart[m]+428]
        rlatstar = np.array([np.where(rmap==s)[0][0] for s in rmap if 'LAT/LON1/LON2/DLON/H' in s])
        for n in range(len(lat)):
            rms[n,:,m] = np.array([s for d in splitt(rmap[rlatstar[n]+1:rlatstar[n]+6]) for s in d],dtype='float')
        t.append(map_t)
    return Time(t,location=('0d', '0d')), long, lat, tec*0.1, rms*0.1

##############################################################################################

if __name__=='__main__':
    main()

##############################################################################################

