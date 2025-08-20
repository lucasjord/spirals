#!/usr/bin/env python3

import numpy as np

'''
Convert an az/el to parallactic angle (and latitude)
'''

def calc_sind(az_r,el_r,lat_r):
    sind = np.sin(el_r)*np.sin(lat_r) + np.cos(el_r)*np.cos(lat_r)*np.cos(az_r)
    return sind

def calc_HA(az_r,el_r,lat_r):
    nom_HA = np.sin(np.pi-az_r)
    dom_HA = np.cos(np.pi-az_r)*np.sin(lat_r) + np.tan(el_r)*np.cos(lat_r)
    HA = np.arctan2(nom_HA,dom_HA)
    return HA

def calc_PA(az_r,el_r,lat_r):
    sind = calc_sind(az_r,el_r,lat_r)
    HA   = calc_HA(az_r,el_r,lat_r)
    nom_PA = np.sin(HA)
    dom_PA = np.tan(lat_r)*np.cos(np.arcsin(sind)) - sind*np.cos(HA)
    PA = np.arctan2(nom_PA,dom_PA)
    return PA

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("azimuth",type=float)
    parser.add_argument("elevation",type=float)
    parser.add_argument("-t","--telescope",type=str,default=None)
    parser.add_argument("-l","--latitude",type=float,default=-999.)
    args = parser.parse_args()
    if args.telescope==None and args.latitude==-999.9:
        sys.exit('Need to specify telescope or latitude')
    tel2lat = {'cd':-31.87,'hb':-42.80,'ke':-14.375,'yg':-29.04}

    if args.telescope!=None:
        try: latitude = tel2lat[args.telescope]/57.2
        except IndexError:
            sys.exit('Unknown telescope {}'.format(args.telescope))
    elif args.latitude!=-999.9:
        latitude = args.latitude/57.2
    az_r = args.azimuth/57.2
    el_r = args.elevation/57.2
    #######
    PA = calc_PA(az_r,el_r,latitude)
    print(PA*57.2)

if __name__=='__main__':
    main()