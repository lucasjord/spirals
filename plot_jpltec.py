#!/usr/bin/env python3

'''
Written by Lucas J. Hyland 2021/09/21
University of Tasmania, Australia
'''

import matplotlib.pyplot as plt
import numpy as np
import struct, os, math, random, matplotlib, datetime, argparse, difflib
from matplotlib import rc, font_manager
from numpy import pi, arcsin, cos, sin, percentile, exp

from mpl_toolkits.basemap import Basemap

from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u

import matplotlib.animation as animation


matplotlib.rcParams.update({'font.size': 14})
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
#rc('text', usetex=True)

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
                        help='Maser number or name e.g. G232.62 or s1',default='G232.62')
    parser.add_argument('-l','--lower',
                        help='Lowest TEC to show', type=float,default=-2)
    parser.add_argument('-u','--upper',
                        help='Highest TEC to show', type=float,default=40)
    parser.add_argument('-R','--rms',
                        help='Do RMS instead of total', action="count",default=0)
    args = parser.parse_args()

    '''
    Load in specified TEC map. Only can load in one at a time atm.
    '''
    t, longi, lati, tec, rms = read_jpltec(args.tec_map)
    '''
    Get maser of interest. Defaults to G232.62 aka s1 at the moment.
    '''
    maser  = get_maser(args.maser)
    source = SkyCoord(ra=maser.ra,dec=maser.dec,unit=(u.hourangle,u.deg))
    gst    = 15*t.sidereal_time('apparent').value
    
    # deletes old animations if they exist (they shouldn't)
    try:
        animate1.event_source.stop()
    except NameError:
        ''
    '''
    Makes basemap of Australia-NZ area. Currently hard-coded.
    '''
    x, y    = np.meshgrid(longi, lati)
    fig, ax = plt.subplots(1,figsize=(14,7))
    ax = Basemap(llcrnrlon=70.,llcrnrlat=-66.,urcrnrlon=220.,urcrnrlat=5.,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='l',projection='merc',ax=ax)
    ax.drawcoastlines()
    ax.drawstates()
    ax.drawparallels(np.arange(-90,90,10),labels=[1,1,0,1]);
    ax.drawmeridians(np.arange(-180,180,15),labels=[1,1,0,1]);
    '''
    Put antenna locations on map. Currently hard-coded for WASCI
    '''
    ant_lat = [-42.804,-31.868,-29.047,-14.375,-36.43]
    ant_lon = [147.440,133.809,115.350,132.150,174.66]
    col = ['g','orange','b','r','y']
    for i in range(len(ant_lat)):
        ax.scatter(ant_lon[i],ant_lat[i],latlon=True,color=col[i],zorder=99)
    '''
    Put on axis objects- 0:TEC colormap, 1:source position, 2:date/time, 
    3:antenna-source LoS (W.I.P) and colorbar (fixed -2 to 40)
    '''
    m = 0
    if args.rms>0:
        plot0 = [ax.pcolor(x,y,rms[:,:,m], latlon=True,cmap='magma',vmin=args.lower,vmax=args.upper)]
    else:
        plot0 = [ax.pcolor(x,y,tec[:,:,m], latlon=True,cmap='magma',vmin=args.lower,vmax=args.upper)]

    plot1 = [ax.scatter([source.ra.deg-gst[m]],[source.dec.deg],latlon=True,color='w',s=50,zorder=99)]
    plot2 = [fig.text(0.175,0.15,s=f'{t.iso[m]}',color='w',size=12)]
    plot3 = []
    #for i in range(len(a_lat)):
    #    plot3.append(ax.drawgreatcircle(a_lon[i],a_lat[i],
    #    source.ra.deg-gst[m],source.dec.deg,color='w',ls='-.',lw=0.25))
    c = plot0[0]
    if args.rms>0: ax.colorbar(c,pad='8%',label=R'$I_e$ (TECU)')
    else:          ax.colorbar(c,pad='8%',label=R'$\delta I_e$ (TECU)')
    '''
    Create animation. 
    '''
    if args.rms>0:
        animate1 = animation.FuncAnimation(fig, update, range(tec.shape[2]), 
                                   fargs=(t,gst,source,x,y,rms,plot0,plot1,plot2,plot3,fig,ax,args),interval=500)
    else:
        animate1 = animation.FuncAnimation(fig, update, range(tec.shape[2]), 
                               fargs=(t,gst,source,x,y,tec,plot0,plot1,plot2,plot3,fig,ax,args),interval=500)       
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

class Maser:
    kind = 'maser' 
    def __init__(self, name, ra, dec, vel, flux, alias):
        self.name  = name  
        self.ra    = ra
        self.dec   = dec
        self.vel   = float(vel)
        self.cflux = float(flux)
        self.alias = alias

def get_maser(template):
    masers = np.array([
        ['G232.620+0.996','07:32:09.79','-16:58:12.4', 22.9, 11.5,'s1' ],
        ['G287.371+0.644','10:48:04.44','-58:27:01.0', -1.9, 21.9,'s3' ],
        ['G309.921+0.479','13:50:41.78','-61:35:10.2',-57.9, 57.6,'s4' ],
        ['G323.740-0.263','15:31:45.45','-56:30:50.1',-50.4,346.6,'s5' ],
        ['G327.402+0.445','15:49:19.50','-53:45:13.9',-82.9, 37.7,'s6' ],
        ['G328.254-0.532','15:57:59.75','-53:58:00.4',-36.8, 20.9,'s7' ],
        ['G328.808+0.633','15:55:48.45','-52:43:06.6',-44.4, 30.6,'s8' ],
        ['G339.622-0.121','16:46:05.99','-45:36:43.3',-33.2, 23.3,'s9' ],
        ['G339.884-1.259','16:52:04.67','-46:08:34.2',-35.6,424.1,'s10'],
        ['G345.505+0.348','17:04:22.91','-40:44:21.7',-14.1, 24.5,'s11'],
        ['G291.274-0.709','11:11:53.35','-61:18:23.7',-30.7, 10.7,'s14'],
        ['G299.772-0.005','12:23:48.97','-62:42:25.3', -6.7, 12.3,'s15'],
        ['G318.948-0.196','15:00:55.40','-58:58:52.1',-36.3, 12.5,'s16'],
        ['G326.475+0.703','15:43:16.64','-54:07:14.6',-38.4, 13.5,'s17'],
        ['G328.237-0.547','15:57:58.28','-53:59:22.7',-44.7, 41.9,'s18'],
        ['G329.029-0.205','16:00:31.80','-53:12:49.6',-36.1, 11.1,'s19'],
        ['G332.295+2.280','16:05:41.72','-49:11:30.3',-23.7, 10.5,'s20'],
        ['G337.920-0.456','16:41:06.05','-47:07:02.5',-38.6, 12.7,'s21'],
        ['G345.010+1.792','16:56:47.58','-40:14:25.8',-17.0, 14.2,'s22'],
        ['G348.550-0.979','17:19:20.41','-39:03:51.6',-10.4, 10.5,'s23'],
        ['G352.630-1.067','17:31:13.91','-35:44:08.7', -3.3, 17.6,'s24']])
    try:
        match = difflib.get_close_matches(template, masers[:,0])[0]
        index = masers[:,0]==match
        return Maser(*masers[index,:][0])
    except IndexError:
        try: 
            match = difflib.get_close_matches(template, masers[:,-1])[0]
            index = masers[:,-1]==match
            return Maser(*masers[index,:][0])
        except IndexError:
            print('Cannot identify maser, defaulting to G232.62')
            return Maser(*masers[0,:])

def update(m,t,gst,source,x,y,mapp,plot0,plot1,plot2,plot3,fig,ax,args):
    plot0[0].remove()
    plot0[0] = ax.pcolor(x,y,mapp[:,:,m], latlon=True,cmap='magma',vmin=args.lower,vmax=args.upper)
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

