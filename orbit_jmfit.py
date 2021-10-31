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
                        help="CLEAN loop gain, 0-1. Default 0.3.",
                        type=float,default=0.3)
    parser.add_argument("-n","--niter",
                        help="Number of CLEAN iterations. Default 200.",
                        type=int,default=200)
    parser.add_argument("-c","--cell",
                        help="Pixel size in image grid, units arcsec.",
                        type=float,default=0.0001)
    parser.add_argument("-i","--imsize",
                        help="Image size in pixels. 256, 512, 1024, 2048 etc. Default 1024.",
                        type=int,default=1024)
    parser.add_argument("-q","--seq",
                        help="Sequence of created images. Default 99. Will overwrite <seq> images if they exist.",
                        type=int,default=99)
    parser.add_argument("-z","--nozap",
                        help="Zap beam files? Enable to keep IBM001 entries.",
                        action="count",default=0)
    parser.add_argument("-d","--delete",
                        help="Enable to delete ICL001 entries after fitting.",
                        action="count",default=0)

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
                        if args.verbosity>0: print("Found file {}.{}.{}.{}".format(args.experiment.upper()+ftype,expclass,disk+1,seq+1))
                        success = success + 1
                        break
                    else: 
                        if args.verbosity>1: print("Did not find file {}.{}.{}.{}".format(args.experiment.upper()+ftype,expclass,disk+1,seq+1))
                        continue
    if success<=0: sys.exit("No catalogues with INNAME {} in AIPSid {}!".format(args.experiment.upper(),args.aipsid))
    else: 
        if args.verbosity>0: print("Found files in {}, continuing...".format(args.aipsid))
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
        if AIPSUVData(calibrators[i],"SPLIT",1,1).exists():
            #if args.verbosity>0: print("Found file {}.SPLIT.1.1".format(calibrators[i]))
            cdata.append(AIPSUVData(calibrators[i],"SPLIT",1,1))
        else:
            print("Could not find catalogue {}.SPLIT.1.1".format(calibrators[i]))


    '''
        For each valid split file, we are going to image it.
    '''
    orig_stdout = sys.stdout
    null = open(os.devnull, 'w')
    for cal_split in cdata:
        if AIPSImage(cal_split.name,"ILC001",1,args.seq).exists():
            if args.verbosity>0: 
                print("Deleting old image {}.ILC001.1.{}".format(cal_split.name,args.seq))
            AIPSImage(cal_split.name,"ILC001",1,args.seq).zap()
        sys.stdout = null
        _image(cal_split,args.niter,args.gain,args.cell,args.imsize,args.seq)
        sys.stdout = orig_stdout
        if not args.nozap>0: _zapbeam(cal_split.name,args.seq)
        

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
        if AIPSImage(cal_split.name,"ICL001",1,args.seq).exists(): 
            wim = WAIPSImage(cal_split.name,"ICL001",1,args.seq)
        else: 
            continue
        wim.squeeze()
        rms = 3*wim.pixels.std()
        if rms>wim.header.datamax: continue
        ew  = (X-wim.header.crpix[0])*-1e-4
        ns  = (Y-wim.header.crpix[1])*+1e-4     
        ind = np.where(wim.pixels==wim.header.datamax)
        snr[i] = wim.header.datamax/rms
        x[i]   = ew[ind]
        y[i]   = ns[ind]
        print("{} {} {}".format(wim.name,ew[ind][0],ns[ind][0]))


        if args.delete>0: wim.zap()
    
    x_off1 = (x*snr**1).sum()/((snr**1).sum())
    y_off1 = (y*snr**1).sum()/((snr**1).sum())

    x_off2 = (x*snr**2).sum()/((snr**2).sum())
    y_off2 = (y*snr**2).sum()/((snr**2).sum())

    print('\n#############################################')
    print(" Target shift is ra={0:+6.4f}, dec={1:+6.4f}".format(-x_off2, -y_off2))
    print("    Add shift as given to maser pos_shift")
    print('#############################################\n')

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