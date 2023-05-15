#!/home/observer/anaconda2/bin/ParselTongue

###############################################################

from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
from Wizardry.AIPSData import AIPSUVData as WAIPSUVData, AIPSImage as WAIPSImage
import numpy as np, os, sys, matplotlib.pyplot as plt, fileinput
from matplotlib import rc, font_manager

np.seterr(divide='ignore', invalid='ignore')

from astropy.coordinates import SkyCoord
import astropy.units as u

from SPiRALS_s001m import *

###############################################################

def main():
    caldata = AIPSUVData(outname[1],outclass[1],outdisk[1],1)

    if not os.path.exists('./multiview'):
       os.makedirs('./multiview')

    caldata2 = AIPSUVData(outname[1],'SPLAT',outdisk[1],1)
    if caldata2.exists(): caldata2.zap()
    runsplat(caldata,sources=target,stokes='I',docal=1)

    runcalib(caldata2,docal=1,snver=1,solmode='P',soltype='L1R',aparm7=1)

    linedata = AIPSUVData(outname[2],outclass[2],outdisk[2],2)
    runcalib(linedata,docal=1,snver=5,solmode='P',soltype='L1R',aparm7=1,chan=channel)

    make_multiviewcontrol(caldata2,linedata,mv_window=mvwin)

    if os.path.exists('./multiview/multiview.TBOUT'): os.remove('./multiview/multiview.TBOUT')
    runtbout(caldata2,'SN',1,'./multiview/multiview.TBOUT')

    if os.path.exists('./multiview/target.TBOUT'): os.remove('./multiview/target.TBOUT')
    runtbout(linedata,'SN',1,'./multiview/target.TBOUT')

    replace('./multiview/multiview.TBOUT',"'INDE'",'-999.9')
    replace('./multiview/target.TBOUT',"'INDE'",'-999.9')

    os.chdir('./multiview/')

    mv_task = os.system(multiview_path+'multiview_4x')
    if mv_task==0:
        print 'MULTV: Appears to have ended successfully'
    else: 
        print 'Multiview failed, please check data'
        raise ValueError

    os.chdir('../')

    runtbin(linedata,'./multiview/target.TBIN')
    runclcal(linedata,6,8,9,'',1,refant)
    runmaimagr(linedata,calsource,200,0.1e-3,1024,channel,1,'',uvwtfn,robust,beam)


###############################################################

def get_file(path):
    #opens and external file and makes it into a list
    fopen = path
    f=open(fopen, 'r+')
    g=list(f)
    g=map(lambda s: s.strip(), g)
    return g

def splitt(old_list):
    #splits the list entries into sublists
    new_list=[]
    for i in old_list:
        new_list+=[i.split()]
    return np.array(new_list)

def runsplat(indata,sources=[''],stokes='HALF',docal=-1):
    spalt = AIPSTask('splat')
    spalt.indata      = indata
    spalt.sources[1:] = sources
    spalt.stokes      = stokes
    spalt.docalib     = docal
    spalt.gainu       = 0
    spalt.aparm[1:]   = [2,0]
    spalt()

def runclcal(indata,snver,gainver,gainuse,interpol,dobtween,refant):
    clcal          = AIPSTask('CLCAL')
    clcal.indata   = indata
    clcal.refant   = refant
    clcal.snver    = snver
    clcal.inver    = 0
    clcal.gainver  = gainver
    clcal.gainuse  = gainuse
    clcal.interpol = interpol
    clcal.dobtween = dobtween
    clcal()

def runmaimagr(indata,source,niter,cz,iz,channel,docal,imna,uvwtfn,robust,beam):
    if imna=='':
        outname = source
    else:
        outname = source[:11-len(imna)]+imna

    imagr = AIPSTask('IMAGR')

    if indata.header['naxis'][3]>1:
        imagr.bif = input('Enter IF: ')
        imagr.eif = imagr.bif
    else:
        imagr.eif = 1
        imagr.bif = 1
    imagr.indata       = indata
    imagr.docal        = docal
    imagr.sourc[1:]    = [source]
    imagr.uvwtfn       = uvwtfn
    imagr.robust       = robust
    imagr.nchav        = 1
    imagr.bchan        = channel
    imagr.echan        = channel
    imagr.cellsize[1:] = [cz,cz]
    imagr.imsize[1:]   = [iz,iz]
    imagr.outna        = outname
    imagr.niter        = niter
    imagr.outdisk      = indata.disk
    imagr.dotv         = -1
    imagr.bmaj         = beam[0]
    imagr.bmin         = beam[1]
    imagr.bpa          = beam[2]
    imagr()
    
def runcalib(indata,sources=[''],gainuse=0,docal=-1,snver=0,solmode='',soltype='',aparm7=0,chan=0):
    calib             = AIPSTask('CALIB')
    calib.indata      = indata
    calib.calsour[1:] = sources
    calib.gainu       = gainuse
    calib.snver       = snver
    calib.solmode     = solmode
    calib.soltype     = soltype
    calib.aparm[3]    = 1
    calib.aparm[7]    = aparm7
    calib()

def runtbout(indata,inext,inver,outtext):
    tbout         = AIPSTask('TBOUT')
    tbout.indata  = indata
    tbout.inext   = inext
    tbout.inver   = inver
    tbout.outtext = outtext
    tbout.docrt   = 999
    tbout()

def runtbin(outdata,intext):
    tbin          = AIPSTask('TBIN')
    tbin.outdata  = outdata
    tbin.intext   = intext
    tbin()
    
def make_multiviewcontrol(indata,in2data,mv_window=30.,outfile='multiview/multiview_control.inp'):
    if os.path.exists(outfile): os.remove(outfile)
    #making su tables
    su_cal = [[s.id__no, s.source, s.raepo, s.decepo] for s in indata.table('SU',1)]
    su_mas = [[s.id__no, s.source, s.raepo, s.decepo] for s in in2data.table('SU',1) if "G" in s.source][0]
    with open(outfile,'w+') as f:
        print >> f, '! {} QSO offsets from target source {}'.format(outname[0][:-2],su_mas[1].strip(' '))
        print >> f,'! '
        print >> f,'! Number of IFs in multiview.TBOUT file'
        print >> f,'   {}'.format(indata.header['naxis'][3])
        print >> f,'! '
        print >> f,'! Time window (min) multiview fitting (eg, 15 => +/-15 min from a target time)'
        print >> f,'  {}'.format(mv_window)
        print >> f,'! Following source numbers from LISTR(scan)...NEEDS TO BE CORRECT'
        print >> f,'! Src#   Name           dX     dY  in degrees'
        for cal in su_cal:
            print >> f,' {:2.0f}     {:<10s}    {:>5.2f}  {:>5.2f}'.format(cal[0],cal[1].strip(' '),
                                                                      (cal[2]-su_mas[2])*np.cos(su_mas[3]*np.pi/180.),
                                                                      cal[3]-su_mas[3])
            
def replace(file, searchExp, replaceExp):
    for line in fileinput.input(file, inplace=1):
        line = line.replace(searchExp, replaceExp)
        sys.stdout.write(line)
        
def mprint(intext,logfile):
    print intext
    f=open(logfile, 'a')
    f.writelines(intext+'\n')
    f.close()

###############################################################

if __name__=='__main__':
    main()

###############################################################