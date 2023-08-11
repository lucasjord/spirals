#!/usr/bin/env ParselTongue

####################################################################################################

'''
    The only purpose of this script is to find, image, jmfit multiview calibrators then return 
    average offset. 

    Lucas J. Hyland
'''

####################################################################################################

from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
from Wizardry.AIPSData import AIPSUVData as WAIPSUVData, AIPSImage as WAIPSImage
import numpy as np
import os, pdb, sys, argparse

####################################################################################################

def main():
    print('\n\n')

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbosity",
                        help="increase output verbosity",
                        action="count",default=0)
    parser.add_argument("experiment",
                        help="experiment code e.g. X420A. Case sensitive.",
                        type=str)
    parser.add_argument("aipsid",
                        help="AIPS USERID e.g. 666. Between 1-46656.",
                        type=int)
    parser.add_argument("-g","--gain",
                        help="CLEAN loop gain, 0-1. Default 0.1",
                        type=float,default=0.1)
    parser.add_argument("-n","--niter",
                        help="Number of CLEAN iterations. Default 200.",
                        type=int,default=200)
    parser.add_argument("-c","--cell",
                        help="Pixel size in image grid, units arcsec.",
                        type=float,default=0.0001)
    parser.add_argument("-i","--imsize",
                        help="Image size in pixels. 256, 512, 1024, 2048 etc. Default 1024.",
                        type=int,default=1024)
    parser.add_argument("-s","--seq",
                        help="Sequence of created images. Default 99. Will overwrite <seq> images if they exist.",
                        type=int,default=99)
    parser.add_argument("-z","--nozap",
                        help="Zap beam files? Enable to keep IBM001 entries.",
                        action="count",default=0)
    parser.add_argument("-d","--delete",
                        help="Enable to delete ICL001 entries after fitting.",
                        action="count",default=0)
    parser.add_argument("-D","--disk",
                        help="Disk where splits are located",
                        type=int,default=1)

    args = parser.parse_args()

    AIPS.userno = args.aipsid

    '''
        First check whether experiment exists in AIPS at given AIPS USER#
    '''
    success = 0
    for ftype in ["_G","_C","_L"]:
        for expclass in ["UVDATA"]:
            for disk in range(len(AIPS.disks)-1):
                for seq in range(2):
                    indata = AIPSUVData(args.experiment.upper()+ftype,expclass,disk+1,seq+1)
                    if indata.exists():
                        if args.verbosity>0:
                            print("Found file {}.{}.{}.{}".format(args.experiment.upper()+ftype,
                                   expclass, disk+1, seq+1))
                        success = success + 1
                        break
                    else:
                        if args.verbosity>1: print("Did not find file {}.{}.{}.{}".format(
                                                    args.experiment.upper()+ftype,
                                                    expclass, disk+1, seq+1))
                        continue
    if success<=0:
        sys.exit("No catalogues with INNAME {} in AIPSid {}!".format(args.experiment.upper(),
                                                                     args.aipsid))
    else:
        if args.verbosity>0:
            print("Found files in {}, continuing...".format(args.aipsid))
        pass


    '''
        Check whether data has been calibrated sufficiently
    '''
    if AIPSUVData(args.experiment.upper()+"_L","UVDATA",1,2).exists():
        if args.verbosity>0: print("Found file {}.UVDATA.1.2".format(args.experiment.upper()+"_L"))
        indata = AIPSUVData(args.experiment.upper()+"_L","UVDATA",1,2)

    elif AIPSUVData(args.experiment.upper()+"_L","UVDATA",2,2).exists():
        if args.verbosity>0: print("Found file {}.UVDATA.2.2".format(args.experiment.upper()+"_L"))
        indata = AIPSUVData(args.experiment.upper()+"_L","UVDATA",2,2)

    else: sys.exit("Cannot find CVELed data {}.UVDATA.*.2".format(args.experiment.upper()+"_L"))

    try:
        check_sncl(indata,4,8)
    except RuntimeError:
        sys.exit("Data is not calibrated sufficiently for this step CL>=8, SN>=4 required")


    '''
        Now we need to identify matching source names. Will assume J=calibrator and G=target
    '''
    calibrators = [s for s in indata.sources if "J" in s]
    if args.verbosity>0: print("Found {} calibrators.".format(len(calibrators)))
    if args.verbosity>1: print("Calibrators: {}".format(calibrators))
    targets     = [s for s in indata.sources if "G" in s]
    if len(targets)>1: sys.exit("More than one possible target found, targets: {}".format(targets))
    target     = targets[-1]
    if args.verbosity>0: print("Found target {}".format(target))


    '''
        Now look for SPLIT catalogue entries. They should be in disk 1
    '''
    cdata = []
    for i in range(len(calibrators)):
        if AIPSUVData(calibrators[i],"SPLIT",args.disk,1).exists():
            if args.verbosity>0: print("Found file {}.SPLIT.{}.1".format(calibrators[i],args.disk))
            cdata.append(AIPSUVData(calibrators[i],"SPLIT",args.disk,1))
        else:
            #pdb.set_trace()
            print("Could not find catalogue {}.SPLIT.{}.1".format(calibrators[i],args.disk))


    '''
        For each valid split file, we are going to image it.
    '''
    orig_stdout = sys.stdout
    null = open(os.devnull, 'w')
    if args.verbosity>0:
    	for cal_split in cdata:
    		print cal_split

#    print args
#    pdb.set_trace()
    for cal_split in cdata:
        if AIPSImage(cal_split.name,"ILC001",args.disk,args.seq).exists():
            print("Deleting old image {}.ILC001.{}.{}".format(cal_split.name,args.disk,args.seq))
            AIPSImage(cal_split.name,"ILC001",args.disk,args.seq).zap()
        #if args.verbosity==0: sys.stdout = null
        _image(cal_split, args.niter, args.gain, args.cell, args.imsize, args.seq)
        #if args.verbosity==0: sys.stdout = orig_stdout
        if not args.nozap>0: _zapbeam(cal_split.name,args.disk,args.seq)


    '''
        Now to look for peaks in the emission, store values and RMS.
        We will be using Wizardry instead of JMFIT
    '''

    X, Y  = np.meshgrid(np.arange(1,args.imsize+.5,1),np.arange(1,args.imsize+.5,1))
    x     = np.zeros(shape=(len(cdata),))
    y     = np.zeros(shape=(len(cdata),))
    snr   = np.zeros(shape=(len(cdata),))
    i=-1
    print ''
    for cal_split in cdata:
        i=i+1
        if AIPSImage(cal_split.name,"ICL001",args.disk,args.seq).exists(): 
            wim = WAIPSImage(cal_split.name,"ICL001",args.disk,args.seq)
        else: 
            continue
        wim.squeeze()
        rms = 3*wim.pixels.std()
        if rms>wim.header.datamax: continue
        ew  = -(X-wim.header.crpix[0])*args.cell
        ns  =  (Y-wim.header.crpix[1])*args.cell
        ind = np.where(wim.pixels==wim.header.datamax)
        snr[i] = wim.header.datamax/rms
        x[i]   = ew[ind]
        y[i]   = ns[ind]
        print("{0:s} {1:+7.4f} {2:+7.4f} {3:4.1f}".format(wim.name,
        													ew[ind][0],
        													ns[ind][0],
        													snr[i]))
        if args.delete>0: wim.zap()

    if args.verbosity>0:
        print ''
        print -x
        print y
        print snr

    x_off1 = np.nansum(x*snr**1)/np.nansum(snr**1)
    y_off1 = np.nansum(y*snr**1)/np.nansum(snr**1)

    x_off2 = np.nansum(x*snr**2)/np.nansum(snr**2)
    y_off2 = np.nansum(y*snr**2)/np.nansum(snr**2)

    print('\n#############################################')
    print(" Target shift is ra={0:+7.4f}, dec={1:+7.4f}".format(-x_off2, -y_off2))
    print("    Add shift as given to maser pos_shift")
    print('#############################################\n')

###################################################
### AIPS TASKS DEFINITIONS
###################################################

def _image(indata,niter=200,gain=0.1,cell=1.0e-4,imsize=1024,seq=99):
    imagr            = AIPSTask("imagr")
    imagr.default
    imagr.indata     = indata
    imagr.docalib    = -1
    imagr.outseq     = seq
    imagr.outname    = indata.name
    imagr.outdisk    = indata.disk
    imagr.cell[1:]   = [cell,cell]
    imagr.imsize[1:] = [imsize,imsize]
    imagr.gain       = gain
    imagr.niter      = niter
    imagr.dotv       = -1
    #imagr.inp()
    imagr()

def _zapbeam(inname,disk,seq=99):
    beam=AIPSImage(inname,"IBM001",disk,seq)
    if beam.exists(): beam.zap()


def check_sncl(indata,sn,cl):
    if (indata.table_highver("AIPS CL")>=cl and
        indata.table_highver("AIPS SN")>=sn):
            print("Tables are okay")
    
    if indata.table_highver("AIPS CL")<cl:
        print('CL bad '+indata.table_highver("AIPS CL"),cl)
        raise RuntimeError("Not enough CL tables")

    if indata.table_highver("AIPS SN")<sn:
        print('SN bad '+indata.table_highver("AIPS SN"),sn)
        raise RuntimeError("Not enough SN tables")


###################################################
###################################################

if __name__=="__main__":
    main()

###################################################
###################################################
