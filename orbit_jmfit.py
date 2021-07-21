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
import numpy as np
import os, pdb, sys, argparse

####################################################################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose",
                        help="increase output verbosity",
                        action="count",default=0)
    parser.add_argument('experiment',
                        help='experiment code e.g. X420A. Case sensitive.',
                        type=str)
    parser.add_argument('aipsid',
                        help='AIPS USERID e.g. 666. Between 1-46656',
                        type=int)
    args = parser.parse_args()

    AIPS.userno = args.aipsid

    '''
        First check whether experiment exists in AIPS at given AIPS USER#
    '''
    success = 0
    for ftype in ['_G','_C','_L']
        for expclass in ['UVDATA']:
            for disk in range(3):
                for seq in range(3):
                    indata = AIPSUVData(args.experiment+ftype,expclass,disk,seq)
                    if indata.exists(): 
                        #if arg.verbosity>0: print(f'Found file {args.experiment+ftype}.{expclass}.{disk}.{seq}')
                        success = success + 1
                        break
                    else: 
                        #if arg.verbosity>1: print(f'Did not find file {args.experiment+ftype}.{expclass}.{disk}.{seq}')
                        continue
    if success<=0: sys.exit('No catalogues with INNAME {} in AIPSid {}!'.format(args.experiment,args.aipsid))
    else: 
        if arg.verbosity>0: print('Found files in {}, continuing...'.format(args.aipsid))
        pass
    
    '''
        Check whether data has been calibrated sufficiently
    '''
    if AIPSUVData(args.experiment+'_L','UVDATA',1,2).exists():
        if arg.verbosity>0: print('Found file {}.'UVDATA'.1.2'.format(args.experiment+'_L'))
        indata = AIPSUVData(args.experiment+'_L','UVDATA',1,2)

    elif AIPSUVData(args.experiment+'_L','UVDATA',2,2).exists():
        if arg.verbosity>0: print('Found file {}.'UVDATA'.2.2'.format(args.experiment+'_L'))
        indata = AIPSUVData(args.experiment+'_L','UVDATA',2,2)

    else: sys.exit("Cannot find CVELed data {}.'UVDATA'.*.2".format(args.experiment+'_L'))

    try:
        check_sncl(indata,8,4)
    except RuntimeError:
        sys.exit(f'Data is not calibrated sufficiently for this step CL>=8, SN>=4 required')

    '''
        Now we need to identify matching source names. Will assume J=calibrator and G=target
    '''
    calibrators = [s for s in indata.sources if "J" in s]
    if arg.verbosity>0: print('Found {} calibrators.'.format(len(calibrators)))
    if arg.verbosity>1: print('Calibrators: {}'.format(calibrators))
    targets     = [s for s in indata.sources if "G" in s]
    if len(targets)>1: sys.exit('More than one target found, targets: {}'.format(targets))
    target     = targets[-1]
    if arg.verbosity>0: print('Found target {}'.format(target))
    '''
        Now look for SPLIT catalogue entries. They should be in disk 1
    '''
    cdata = []   
    for i in range(len(calibrators)):
        cdata.append(AIPSUVData(calibrators[i],'SPLIT',1,1))

    for i in orbits_of_interest:
        sources = [middle_sources[i]]+orbit_sources[i]
        if do_del==1:
            for k in range(len(sources)):
                _zapimage(sources[k],1)
                _zapimage(sources[k],2)
                _zapimage(sources[k],3)
                _zapimage(sources[k],4)
                _zapimage(sources[k],5)
                _zapimage(sources[k],6)
                _zapimage(sources[k],7)
                _zapimage(sources[k],8)
                _zapimage(sources[k],9)
                _zapimage(sources[k],10)

    if do_image1==1:
        for i in orbits_of_interest:
            sources = [middle_sources[i]]+orbit_sources[i]
            for j in range(len(orbit_sources[i])):
                _image(middle_sources[i],orbit_sources[i][j],2)
                _zapbeam(orbit_sources[i][j])
            _image(middle_sources[i],middle_sources[i],2)
            _zapbeam(middle_sources[i])

    if do_jmfit==1:
        for i in orbits_of_interest:        
            sources = [middle_sources[i]]+orbit_sources[i]
            for inseq in range(3):
                fitout = './phases/'+middle_sources[i]+'_'+str(inseq+1)+'.jmfit'
                if os.path.exists(fitout): 
                    os.remove(fitout)
                for source in sources:
                    _jmfit(inname=source,inseq=inseq+1,fitout=fitout)

###################################################
### AIPS TASKS DEFINITIONS
###################################################

def check_sncl(indata,sn,cl):
    if (indata.table_highver('AIPS CL')>=cl and
        indata.table_highver('AIPS SN')>=sn):
            print('Tables are okay')
    
    if indata.table_highver('AIPS CL')<cl:
        print(indata.table_highver('AIPS CL'),cl)
        raise RuntimeError('Not enough CL tables')

    if indata.table_highver('AIPS SN')<sn:
        print(indata.table_highver('AIPS SN'),sn)
        raise RuntimeError('Not enough SN tables')

def _jmfit(inname='',inseq=1,fitout=''):
    jmfit           =  AIPSTask('jmfit')
    jmfit.default
    image           =  AIPSImage(inname,'ICL001',EXPDISK,inseq)
    jmfit.indata    =  image
    jmfit.fitout    =  fitout
    jmfit.doprint   = -4
    jmfit.blc[1:]  = [216.0, 216.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    jmfit.trc[1:]  = [296.0, 296.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    jmfit.niter     =  500
    if image.exists(): jmfit()

def check_sncl(indata,sn,cl):
    if (indata.table_highver('AIPS CL')==cl and
        indata.table_highver('AIPS SN')==sn):
            print('Tables are fine')
    
    if indata.table_highver('AIPS CL')<cl:
        print(indata.table_highver('AIPS CL'),cl)
        raise RuntimeError('Not enough CL tables')

    if indata.table_highver('AIPS CL')>cl:
        while indata.table_highver('AIPS CL')>cl:
            indata.zap_table('AIPS CL', 0)

    if indata.table_highver('AIPS SN')<sn:
        print(indata.table_highver('AIPS SN'),sn)
        raise RuntimeError('Not enough SN tables')

    if indata.table_highver('AIPS SN')>sn:
        while indata.table_highver('AIPS SN')>sn:
            indata.zap_table('AIPS SN', 0)

def _zapimage(source,inseq):
    img=AIPSImage(source,'ICL001',1,inseq)
    if img.exists():
        img.zap()

def _zapbeam(source):
    beam=AIPSImage(source,'IBM001',1,1)
    if beam.exists():
        beam.zap()

def _image(midsource,source,Gin):
    imagr=AIPSTask('imagr')
    imagr.default
    imagr.inname=midsource
    imagr.source[1]=source
    imagr.inclass='RING'
    imagr.inseq=2
    imagr.indisk=1
    imagr.docalib=1
    imagr.gainu=Gin
    imagr.outseq=0
    #imagr.bmaj=1.6e-3
    #imagr.bmin=1.6e-3
    #imagr.nbox=1
    #imagr.clbox[1][1:]=[108,108,148,148]
    imagr.antenna[1:]= array
    imagr.baseline[1:]=array
    imagr.outname=source
    imagr.cell[1:]=[0.5e-4,0.5e-4]
    imagr.imsize[1:]=[512,512]
    imagr.gain=0.3
    imagr.niter=200
    imagr.dotv=-1
    imagr()

###################################################
###################################################

def get_file(path):
    #opens and external file and makes it into a list
    fopen = path
    f=open(fopen, 'r+')
    g=list(map(lambda s: s.strip(), list(f)))
    return np.array(g)

def get_experiment():
    files = os.listdir('./')
    reffile=0
    for i in range(len(files)):
        if "LST" in files[i]:
            reffile = files[i]
            break
    if reffile==0:
        sys.exit('No AIPS output file to read from')
    readfile = get_file(reffile)
    aipsid   = readfile[0].split()[2]
    experiment = readfile[1].split()[2]
    return int(aipsid), experiment

###################################################
###################################################

if __name__=='__main__':
    main()

###################################################
###################################################