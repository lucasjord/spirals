#!/home/observer/anaconda2/bin/ParselTongue

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
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose",
                        help="increase output verbosity",
                        action="count",default=0)
    parser.add_argument("experiment",
                        help="experiment code e.g. X420A. Case sensitive.",
                        type=str)
    parser.add_argument("aipsid",
                        help="AIPS USERID e.g. 666. Between 1-46656",
                        type=int)
    parser.add_argument("-g","--gain",
                        help="CLEAN loop gain, 0-1",
                        type=float,default=0.3)
    parser.add_argument("-n","--niter",
                        help="Number of CLEAN iterations",
                        type=int,default=200)
    parser.add_argument("-c","--cell",
                        help="Pixel size in image grid, units arcsec",
                        type=float,default=0.0001)
    parser.add_argument("-i","--imsize",
                        help="Image size in pixels",
                        type=int,default=1024)
    parser.add_argument("-q","--seq",
                        help="Sequence of created images. Default 99. Will overwrite <seq> images if they exist.",
                        type=int,default=99)
    parser.add_argument("-z","--nozap",
                        help="Zap beam files? Enable to NOT delte the beams.",
                        type=count,default=0)

    args = parser.parse_args()

    AIPS.userno = args.aipsid


    '''
        First check whether experiment exists in AIPS at given AIPS USER#
    '''
    success = 0
    for ftype in ["_G","_C","_L"]
        for expclass in ["UVDATA"]:
            for disk in range(3):
                for seq in range(3):
                    indata = AIPSUVData(args.experiment+ftype,expclass,disk,seq)
                    if indata.exists(): 
                        #if arg.verbosity>0: print(f"Found file {args.experiment+ftype}.{expclass}.{disk}.{seq}")
                        success = success + 1
                        break
                    else: 
                        #if arg.verbosity>1: print(f"Did not find file {args.experiment+ftype}.{expclass}.{disk}.{seq}")
                        continue
    if success<=0: sys.exit("No catalogues with INNAME {} in AIPSid {}!".format(args.experiment,args.aipsid))
    else: 
        if arg.verbosity>0: print("Found files in {}, continuing...".format(args.aipsid))
        pass
    

    '''
        Check whether data has been calibrated sufficiently
    '''
    if AIPSUVData(args.experiment+"_L","UVDATA",1,2).exists():
        if arg.verbosity>0: print("Found file {}."UVDATA".1.2".format(args.experiment+"_L"))
        indata = AIPSUVData(args.experiment+"_L","UVDATA",1,2)

    elif AIPSUVData(args.experiment+"_L","UVDATA",2,2).exists():
        if arg.verbosity>0: print("Found file {}."UVDATA".2.2".format(args.experiment+"_L"))
        indata = AIPSUVData(args.experiment+"_L","UVDATA",2,2)

    else: sys.exit("Cannot find CVELed data {}."UVDATA".*.2".format(args.experiment+"_L"))

    try:
        check_sncl(indata,8,4)
    except RuntimeError:
        sys.exit("Data is not calibrated sufficiently for this step CL>=8, SN>=4 required")


    '''
        Now we need to identify matching source names. Will assume J=calibrator and G=target
    '''
    calibrators = [s for s in indata.sources if "J" in s]
    if arg.verbosity>0: print("Found {} calibrators.".format(len(calibrators)))
    if arg.verbosity>1: print("Calibrators: {}".format(calibrators))
    targets     = [s for s in indata.sources if "G" in s]
    if len(targets)>1: sys.exit("More than one possible target found, targets: {}".format(targets))
    target     = targets[-1]
    if arg.verbosity>0: print("Found target {}".format(target))
    

    '''
        Now look for SPLIT catalogue entries. They should be in disk 1
    '''
    cdata = []   
    for i in range(len(calibrators)):
        if AIPSUVData(calibrators[i],"SPLIT",1,1).exists():
            #if arg.verbosity>0: print("Found file {}.SPLIT.1.1".format(calibrators[i]))
            cdata.append(AIPSUVData(calibrators[i],"SPLIT",1,1))
        else:
            print("Could not find catalogue {}.SPLIT.1.1".format(calibrators[i]))


    '''
        For each valid split file, we are going to image it.
    '''
    for cal_split in cdata:
        if AIPSImage(cal_split.name,"ILC001",1,args.seq).exists():
            if arg.verbosity>0: print("Deleting old image {}.ILC001.1.{}".format(cal_split.name,args.seq))
            AIPSImage(cal_split.name,"ILC001",1,args.seq).zap()
        _image(cal_split,args.gain,args.niter,args.cell,args.imsize,args.seq)
        _zapbeam(cal_split.name,args.seq)
    

    '''
        Now to look for peaks in the emission, store values and RMS.
        We will be using Wizardry instead of JMFIT
    '''

    X, Y  = np.meshgrid(np.arange(1,args.imsize+.5,1),np.arange(1,args.imsize+.5,1))
    x   = np.zeros(shape=(len(cdata),))
    y   = np.zeros(shape=(len(cdata),))
    snr = np.zeros(shape=(len(cdata),))
    i=-1
    for cal_split in cdata:
        i=i+1
        if AIPSImage(cal_split.name,"ICL001",1,args.seq).exists(): wim = WAIPSImage(cal_split.name,"ICL001",1,args.seq)
        else: continue
        wim.squeeze()
        rms = 3*wim.pixels.std()
        if rms>wim.header.datamax: continue
        ew  = (X-wim.header.crpix[0])*-1e-4
        ns  = (Y-wim.header.crpix[1])*+1e-4     
        ind = np.where(wim.pixels==wim.header.datamax)
        snr[i] = wim.header.datamax/rms
        x[i]   = ew[ind]
        y[i]   = ns[ind]
    
    x_off1 = (x*snr**1).sum()/((snr**1).sum())
    y_off1 = (y*snr**1).sum()/((snr**1).sum())

    x_off2 = (x*snr**2).sum()/((snr**2).sum())
    y_off2 = (y*snr**2).sum()/((snr**2).sum())

    print("Shift is ra={0:+6.4f}, dec={1:+6.4f}".format(-x_off2, -y_off2))
    print("Add shift as given to maser pos_shift")


###################################################
### AIPS TASKS DEFINITIONS
###################################################

def _image(indata,niter=200,gain=0.3,cell=1.0e-4,imsize=1024,seq=99):
    imagr            = AIPSTask("imagr")
    imagr.default
    imagr.indata     = indata
    imagr.docalib    = -1
    imagr.outseq     = seq
    imagr.outname    = indata.name
    imagr.cell[1:]   = [cell,cell]
    imagr.imsize[1:] = [imsize,imsize]
    imagr.gain       = gain
    imagr.niter      = niter
    imagr.dotv       = -1
    imagr()

def _zapbeam(inname,seq=99):
    beam=AIPSImage(inname,"IBM001",1,seq)
    if beam.exists():
        beam.zap()

def check_sncl(indata,sn,cl):
    if (indata.table_highver("AIPS CL")>=cl and
        indata.table_highver("AIPS SN")>=sn):
            print("Tables are okay")
    
    if indata.table_highver("AIPS CL")<cl:
        print(indata.table_highver("AIPS CL"),cl)
        raise RuntimeError("Not enough CL tables")

    if indata.table_highver("AIPS SN")<sn:
        print(indata.table_highver("AIPS SN"),sn)
        raise RuntimeError("Not enough SN tables")

def _jmfit(inname=",inseq=1,fitout="):
    jmfit           =  AIPSTask("jmfit")
    jmfit.default
    image           =  AIPSImage(inname,"ICL001",EXPDISK,inseq)
    jmfit.indata    =  image
    jmfit.fitout    =  fitout
    jmfit.doprint   = -4
    jmfit.blc[1:]  = [216.0, 216.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    jmfit.trc[1:]  = [296.0, 296.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    jmfit.niter     =  500
    if image.exists(): jmfit()

def check_sncl(indata,sn,cl):
    if (indata.table_highver("AIPS CL")==cl and
        indata.table_highver("AIPS SN")==sn):
            print("Tables are fine")
    
    if indata.table_highver("AIPS CL")<cl:
        print(indata.table_highver("AIPS CL"),cl)
        raise RuntimeError("Not enough CL tables")

    if indata.table_highver("AIPS CL")>cl:
        while indata.table_highver("AIPS CL")>cl:
            indata.zap_table("AIPS CL", 0)

    if indata.table_highver("AIPS SN")<sn:
        print(indata.table_highver("AIPS SN"),sn)
        raise RuntimeError("Not enough SN tables")

    if indata.table_highver("AIPS SN")>sn:
        while indata.table_highver("AIPS SN")>sn:
            indata.zap_table("AIPS SN", 0)


###################################################
###################################################

def get_file(path):
    #opens and external file and makes it into a list
    fopen = path
    f=open(fopen, "r+")
    g=list(map(lambda s: s.strip(), list(f)))
    return np.array(g)

###################################################
###################################################

if __name__=="__main__":
    main()

###################################################
###################################################