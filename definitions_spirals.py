##############################################################################
#       ParselTongue definition file for BeSSeL-Survey calibration           #
#                                                                            #
#                  written by Andreas Brunthaler                             #
#                                                                            #
# Changes:                                                                   #
#                                                                            #
# 2010/10/12 Fixed another bug in VBGLU                                      #
# 2010/10/26 Set dparm[4]=0 in all fringe after AIPS bug was fixed           #
#            Added checkatmos to look for zero zenith delay offsets          #
# 2010/10/29 Added plotatmos to plot ATMOS.FITS file                         #
#            Minor changes in plot_baseline                                  #
# 2010/12/02 Added check_calsource, runsnflg                                 #
# 2010/12/15 Added main script to file                                       #
# 2010/12/20 Added DELZN, fixed bug in split, and added input model          #
#            for fringe                                                      #
# 2010/12/30 Added zero line in delay/rate plot, added possibility for fixed #
#            range, added printout to logfile, added make_name_atmos         #
# 2011/01/03 Changed position shifts from two steps (DEC, the RA) to one     #
# 2011/01/07 Fixed bug in make_name_atmos                                    #
# 2011/02/09 fit_geoblocks output now written also to file                   #
# 2011/04/13 Changed format of ATMOS_NAME.FITS, implemented SX-geoblocks,    #
#            added Mail input for download options, enable fittp of splited  #
#            files                                                           #
# 2011/04/15 Added check for line data in check_sx                           #
# 2011/05/24 Fixed bug in plotting continuum sources and TECOR               #
# 2011/05/30 Allowed multiple maps for fringe input model (nmaps)            #
# 2011/08/09 Fixed bug in SX geoblocks and added niter for imagr             #
# 2011/11/07 Added check that man_pcal scans are identical, added aipsver,   #
#            added dowload of key files, changed wtthresh to 0.45            #
# 2012/11/09 Added smodel, plot_tables and self-calibration for continuum    #
#            sources                                                         #
# 2013/03/13 Added bandpass calibration. First epoch products, smooth during #
#            split                                                           #
# 2013/12/18 Fixed CVEL behaviour for AIPS versions after 2011, and changed  #
#            how orfit gets it velocity. Also changed CPARM[2] in INDXR to 0.#
# 2014/09/04 Added RDBE check                                                #
# 2015/07/07 Changed PRTAB output for 31DEC2015 version                      #
# 2015/07/22 Added dual_geo = 2                                              #
# 2015/08/15 Fixed bug in dual_geo = 2                                       #
# 2015/12/17 Keeping flags in dual-frequency geoblock data                   #
# 2016/03/23 Fixed plotting bug in RDBE_check. Added choice of TECU model    #
# 2016/04/06 Added weightit=3 and aparm(7)=10 in fringegeo, changed delay    #
#            error for SX geoblocks to 3.                                    #
# 2016/09/07 Changed delay window in manual phasecal                         #
# 2017/11/13 Added LBA apcal, various corrections for LBA data               #
# 2018/08/27 Added YG flag in make station_file.inp                          #
# 2019/06/20 Removed fringecal RR/LL averaging                               #
# 2020/02/24 Added extra vlaglu for 8IF cont data - LJH                      #
# 2020/03/30 Changed runpang and runpang2 to not run on HB/KE linears        #
# 2020/10/03 Removed setjy.optype='VCAL' from runcvel_lba - LJH              #
# 2020/12/14 Added run_snplt_diff to track RR-LL phase drifts  - LJH         #
# 2021/01/12 Reversed runpang & runpang2 changes from 2020/03/30             #
#            Edited bug in mafringe for short source names  - LJH            #
# 2021/03/19 Changed cvel to use input or middle channel of band, now        #
#            seperate channel to fringe fitting - LJH                        #
# 2021/10/03 Added inverse Multiview routines - LJH                          #
##############################################################################

version_date='2021/10/03'

from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
from Wizardry.AIPSData import AIPSUVData as WAIPSUVData
import AIPSTV
import AIPS, os, math, time, sys
from pylab import *

from astropy.coordinates import SkyCoord
import astropy.units as u, fileinput

import pdb #debugger

if 'aipsver' in locals() and globals(): AIPSTask.version = aipsver
else: aipsver = AIPSTask.version

##############################################################################
# Start defs

def check_sx(indata,logfile):
    if indata.header.naxis[3]>1:
        fq=indata.table('AIPS FQ',0)
        fq_span=fq[0]['if_freq'][indata.header.naxis[3]-1]-fq[0]['if_freq'][0]
        frq=(indata.header.crval[2]+0.5*fq_span)/1e9
        if 3>indata.header.crval[2]/1e9>2 and fq_span>1e9:
            return True
        return False
    else: return False

def check_clh(indata,logfile):
    if indata.header.naxis[3]>1:
        fq=indata.table('AIPS FQ',0)
        fq_span=fq[0]['if_freq'][indata.header.naxis[3]-1]-fq[0]['if_freq'][0]
        frq=(indata.header.crval[2]+0.5*fq_span)/1e9
        if 8>indata.header.crval[2]/1e9>4 and fq_span>6e8:
            return True
        return False
    else: return False

def split_sx(indata,split):
    uvcop          = AIPSTask('UVCOP')
    uvcop.indata   = indata
    uvcop.flagver  = indata.table_highver('AIPS FG')
    uvcop.bif      = 1
    uvcop.eif      = 4
    uvcop.outname  = indata.header['observer']+'-G-S'
    uvcop.outclass = 'UVDATA'
    uvcop.outdisk  = uvcop.indisk

    outdata_s = AIPSUVData(uvcop.outname,'UVDATA',1,1)
    if split==1:
        if outdata_s.exists():
            outdata_s.zap()
        uvcop.go()

    uvcop.indata   = indata
    uvcop.flagver  = indata.table_highver('AIPS FG')
    uvcop.bif      = 5
    uvcop.eif      = 8
    uvcop.outname  = indata.header['observer']+'-G-X'
    uvcop.outclass = 'UVDATA'
    uvcop.outdisk  = uvcop.indisk

    outdata_x = AIPSUVData(uvcop.outname,'UVDATA',1,1)
    if split==1:
        if outdata_x.exists():
            outdata_x.zap()
        uvcop.go()

    return outdata_s, outdata_x

def split_clh(indata,split):
    uvcop          = AIPSTask('UVCOP')
    uvcop.indata   = indata
    uvcop.flagver  = indata.table_highver('AIPS FG')
    uvcop.bif      = 1
    uvcop.eif      = 4
    uvcop.outname  = indata.header['observer']+'-G-CL'
    uvcop.outclass = 'UVDATA'
    uvcop.outdisk  = uvcop.indisk

    outdata_l = AIPSUVData(uvcop.outname,'UVDATA',1,1)
    if split==1:
        if outdata_l.exists():
            outdata_l.zap()
        uvcop.go()

    uvcop.indata   = indata
    uvcop.flagver  = indata.table_highver('AIPS FG')
    uvcop.bif      = 5
    uvcop.eif      = 8
    uvcop.outname  = indata.header['observer']+'-G-CH'
    uvcop.outclass = 'UVDATA'
    uvcop.outdisk  = uvcop.indisk

    outdata_h = AIPSUVData(uvcop.outname,'UVDATA',1,1)
    if split==1:
        if outdata_h.exists():
            outdata_h.zap()
        uvcop.go()

    return outdata_l, outdata_h

def mprint(intext,logfile):
    print intext
    f=open(logfile, 'a')
    f.writelines(intext+'\n')
    f.close()

def get_ant(data):
    antennas = {}
    for row in data.table('AN', 0):
        antennas[row.nosta] = row.anname[0:2]
    return antennas

def time_to_hhmmss(time):
    day  = int(time)
    if time>1:
        time=time-int(time)
    hour = int(time*24)
    min  = int(60*(time*24-hour))
    sec  = int(60*(60*(time*24-hour)-min))
    return day,hour,min,sec

def  isinde(number):
    INDE = 3140.892822265625
    return  abs(number-INDE)<1e-12

def get_timerange_tab(indata,table,i):
    sn_table = indata.table(table, 0)
    time1 = sn_table[i].time-0.5*sn_table[i].time_interval
    time2 = sn_table[i].time+0.5*sn_table[i].time_interval
    (day1,hour1,min1,sec1)=time_to_hhmmss(time1)
    (day2,hour2,min2,sec2)=time_to_hhmmss(time2)
    timerange = [day1, hour1, min1, sec1, day2, hour2, min2, sec2]
    return timerange

##############################################################################
# Get the day-of-year integer from the year/month/day
#
def get_day_of_year(year, month, day):
    day_of_year_list = [0,31,59,90,120,151,181,212,243,273,304,334]
    doy = day_of_year_list[month-1] + day
    if(month>2):
        if((year&0x3)==0):
            if((year % 100 != 0) or (year % 400 == 0)):
                doy = doy+1
    return doy

##############################################################################
# Get the day of year from the Year, month, day for the start of observations
#
def get_observation_year_month_day(aips_data):
    date_string = aips_data.header.date_obs
    date_list = date_string.split('-')
    year = int(date_list[0])
    month = int(date_list[1])
    day = int(date_list[2])
    return (year, month, day)

##############################################################################
# Get number of days
#
def get_num_days(indata):
    nx_table = indata.table('AIPS NX', 0)
    n        = len(nx_table)
    num_days = int(nx_table[n-1]['time']+1)
    return num_days

##############################################################################
# Get center UTC
#
def get_utc(indata):
    nx_table = indata.table('AIPS NX', 0)
    n        = len(nx_table)
    ut1      = nx_table[0]['time']
    ut2      = nx_table[n-1]['time']
    utc      = (ut1+ut2)/2.
    return utc

##############################################################################
# Get center frequency in GHz
#
def get_center_freq(indata):
    fq = indata.table('AIPS FQ',0)
    naxis = indata.header['naxis']
    if naxis[3]>1:
        fq_span=fq[0]['if_freq'][indata.header.naxis[3]-1]-fq[0]['if_freq'][0]
        frq=(indata.header.crval[2]+0.5*fq_span)
    else:
        frq=indata.header.crval[2]
    return frq

##############################################################################
# Download TEC maps
#
def get_TEC(year,doy,TECU_model):
    year=str(year)[2:4]
    if doy<10:
        doy='00'+str(doy)
    elif doy<100:
        doy='0'+str(doy)
    else:
        doy=str(doy)
    name=TECU_model+doy+'0.'+year+'i'
    if os.path.exists(name):
        print 'File already there.'
    else:
        #path='ftp://cddis.gsfc.nasa.gov/gps/products/ionex/20'+year+'/'+doy+'/'
        path='ftp://gdc.cddis.eosdis.nasa.gov/gnss/products/ionex/20'+year+'/'+doy+'/'
        #os.popen(r'wget -t 30 -O '+name+'.Z '+path+name+'.Z')
        os.popen(r'curl --insecure -O --ftp-ssl '+path+name+'.Z')
        os.popen(r'uncompress -f '+name+'.Z')

def check_geo(indata):
    nx_table = indata.table('AIPS NX', 0)
    n_block  = 1
    n        = len(nx_table)
    b        = []
    b.append(round(nx_table[0]['time']-0.01,2))
    for i in range(1,n):
        if (nx_table[i]['time']-nx_table[i-1]['time'])>0.02:
            n_block=n_block+1
            b.append(round(nx_table[i]['time']-0.01,2))
    b.append(round(nx_table[n-1]['time']+0.01,2))
    return b

##############################################################################
#
def runTECOR(indata,year,doy,num_days,gainuse,TECU_model):
    year=str(year)[2:4]
    if doy<10:
        doy='00'+str(doy)
    if doy<100:
        doy='0'+str(doy)
    else:
        doy=str(doy)
    name=TECU_model+doy+'0.'+year+'i'
    tecor = AIPSTask('TECOR')
    if os.path.exists(name):
        tecor.infile='PWD:'+name
    tecor.indata=indata
    tecor.nfiles=num_days
    tecor.gainuse = gainuse
    tecor.aparm[1:] = [1,0]
    tecor()

##############################################################################
# Download EOP file
#
def get_eop(eop_path):
    if os.path.exists(eop_path+'usno_finals.erp'):
        age = (time.time() - os.stat(eop_path+'usno_finals.erp')[8])/3600
        mprint('usno_finals.erp exists, not downloaded.',logfile)
    else:
        os.popen(r'curl --insecure -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno_finals.erp')   #wget ftp://cddis.gsfc.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno500_finals.erp http://gemini.gsfc.nasa.gov/solve_save ftp://ftp.lbo.us/pub/staff/wbrisken/EOP
        os.popen(r'curl --insecure -O --ftp-ssl ftp://gdc.cddis.eosdis.nasa.gov/vlbi/gsfc/ancillary/solve_apriori/usno500_finals.erp')
        os.popen(r'mv usno500_finals.erp '+eop_path+'usno_finals2.erp')
        os.popen(r'mv usno_finals.erp '+eop_path+'usno_finals.erp')

##############################################################################
#
def runeops(indata, eop_path):
    eops        = AIPSTask('CLCOR')
    eops.indata = indata
    eops.gainver = 2
    eops.gainuse = 3
    eops.opcode  = 'EOPS'
    eops.infile  = eop_path+'usno_finals.erp'
    eops()

##############################################################################
#
def runuvflg(indata,flagfile,logfile):
    if flagfile!='' and os.path.exists(flagfile):
        uvflg        = AIPSTask('UVFLG')
        uvflg.indata = indata
        uvflg.intext = flagfile
        uvflg.opcode = 'FLAG'
        uvflg.go()
    else:
        mprint('No UVFLG file applied.',logfile)

##############################################################################
#
def runindxr(indata):
    indxr = AIPSTask('indxr')
    indxr.indata = indata
    indxr.cparm[1:] = [0, 0, 1./60.]
    indxr()

##############################################################################
#
def runsnflg(indata, inver, calsource):

    flag_sources = [calsource]
    snflg             = AIPSTask('snflg')
    snflg.indata      = indata
    snflg.flagver     = 0
    snflg.inext       = 'SN'
    snflg.inver       = inver
    snflg.optype      = 'JUMP'
    snflg.dparm[1:]   = [57., 10., 0.]
    snflg()

##############################################################################
#
def run_elvflag(indata,elv_min,logfile):
    uvflg        = AIPSTask('UVFLG')
    uvflg.indata = indata
    uvflg.opcode = 'FLAG'
    uvflg.aparm[1:] = [0,elv_min]
    mprint('#####################################',logfile)
    mprint('Flagging data for Elevations < '+str(elv_min),logfile)
    mprint('#####################################',logfile)
    uvflg.go()

##############################################################################
# Download data from archive
def download(user,passw,file,path,filename):
    if user=='nopass':
        os.popen(r'wget --no-check-certificate -O '+path+
                 file+' ftp://ftp.aoc.nrao.edu/e2earchive/'+
                 file, 'w')
        os.popen(r'mv '+path+file+' '+path+filename)
    else:
        os.popen(r'wget --no-check-certificate --http-user '+user+
                 ' --http-passwd '+passw+' -t45 -O '+path+file+
                 ' https://archive.nrao.edu/secured/'+user+'/'+
                 file, 'w')
        os.popen(r'mv '+path+file+' '+path+filename)

##############################################################################
# Get .key and .sum file from archive
def get_key_file(indata,code):
    mon_dir={1:'jan',2:'feb',3:'mar',4:'apr',5:'may',6:'jun',7:'jul',8:'aug'
             ,9:'sep',10:'oct',11:'nov',12:'dec'}
    path0=('http://www.vlba.nrao.edu/astro/VOBS/astronomy/')
    (year,month,day) = get_observation_year_month_day(indata)
    if code=='':
        code_str         = indata.header.observer.lower()
    else:
        code_str      = code
    path1=path0+mon_dir[month]+str(year)[2:4]+'/'+code_str+'/'+code_str+'.key'
    path2=path0+mon_dir[month]+str(year)[2:4]+'/'+code_str+'/'+code_str+'.sum'
    if os.path.exists(code_str+'.key'):
        pass
    else:
        os.popen(r'wget '+path1)
    if os.path.exists(code_str+'.sum'):
        pass
    else:
        os.popen(r'wget '+path2)

##############################################################################
# Load data and index it
#
def loadindx(filepath,filename,outname,outclass,
             outdisk,nfiles,ncount,doconcat,antname,logfile):
    if os.path.exists(filepath+filename):
        mprint('File exists!',logfile)
    else:
        raise RuntimeError('File does not exists!')

    fitld = AIPSTask('FITLD')
    fitld.datain   = filepath+filename
    fitld.outname  = outname
    fitld.outclass = outclass
    fitld.outseq   = 1
    fitld.outdisk  = int(outdisk)
    fitld.ncount   = ncount
    fitld.nfiles   = nfiles
    fitld.doconcat = doconcat
    fitld.clint    = 1./60.
    fitld.wtthresh = 0.45
    if aipsver!='31DEC09':
        if not antname=='LBA': 
            fitld.antname[1:] = [antname]

    data = AIPSUVData(fitld.outname, fitld.outclass,
                      int(fitld.outdisk), int(fitld.outseq))
    if data.exists():
        data.clrstat()
        data.zap()
        mprint('##############################',logfile)
        mprint('Data already there => deleted!',logfile)
        mprint('##############################',logfile)
    else:
        mprint('#########################',logfile)
        mprint('Data not there => read in',logfile)
        mprint('#########################',logfile)

    fitld.go()

    mprint('################################################',logfile)
    mprint(str(data)+' loaded!',logfile)
    mprint('################################################',logfile)

    if data.exists():
        data.zap_table('AIPS CL',1)
        runindxr(data)
        mprint('#################',logfile)
        mprint('Data new indexed!',logfile)
        mprint('#################',logfile)
    else:
        mprint('No!',logfile)

##############################################################################
# Select one of the inner VLBA antennas
#
def select_refant(indata):
    ant=get_ant(indata)
    if 'LA' in indata.antennas:
        refant='LA'
    elif 'PT' in indata.antennas:
        refant='PT'
    elif 'KP' in indata.antennas:
        refant='KP'
    elif 'FD' in indata.antennas:
        refant='FD'
    elif 'EF' in indata.antennas:
        refant='EF'
    elif 'EB' in indata.antennas:
        refant='EB'
    for i in ant:
        if ant[i]==refant:
            j=i
    return j

##############################################################################
# Select best antenna  based on testfringe
#
def select_refant2(indata, logfile):
    ant=get_ant(indata)
    sn=indata.table('AIPS SN', 2)
    sol=range(max(ant))
    for i in range(len(sol)):
        sol[i]=0
    for i in range(0,len(sn)):
        an=sn[i]['antenna_no']
        if (ant[an]=='LA' or ant[an]=='PT' or
            ant[an]=='KP' or ant[an]=='FD' or
            ant[an]=='EF' or ant[an]=='EB'):
            sol[an-1]+=1
    for i in ant:
        if sol[i-1]==max(sol):
            refant=i
    mprint('#################################',logfile)
    mprint('Reference antenna : '+ant[refant],logfile)
    mprint('#################################',logfile)
    return refant

##############################################################################
# Check table versions
#
def check_sncl(indata,sn,cl,logfile):
    if (indata.table_highver('AIPS CL')==cl and
        indata.table_highver('AIPS SN')==sn):
        mprint('SN and CL tables ok.',logfile)

    if indata.table_highver('AIPS CL')<cl:
        print indata.table_highver('AIPS CL'),cl
        raise RuntimeError('Not enough CL tables')

    if indata.table_highver('AIPS CL')>cl:
        mprint('Deleting old CL tables.',logfile)
        while indata.table_highver('AIPS CL')>cl:
            indata.zap_table('AIPS CL', 0)

    if indata.table_highver('AIPS SN')<sn:
        raise RuntimeError('Not enough SN tables')

    if indata.table_highver('AIPS SN')>sn:
        mprint('Deleting old SN tables.',logfile)
        while indata.table_highver('AIPS SN')>sn:
            indata.zap_table('AIPS SN', 0)

##############################################################################
#
def do_band(indata, bandcal, logfile):
    if bandcal==['']:
        print mprint('No Bandpass calibrator selected.', logfile)
        sys.exit()

    if indata.table_highver('AIPS BP')>0:
        mprint('Deleting old BP tables.',logfile)
        while indata.table_highver('AIPS BP')>0:
            indata.zap_table('AIPS BP', 0)

    bpass              = AIPSTask('BPASS')
    bpass.indata       = indata
    bpass.calsour[1:]  = bandcal
    bpass.docal        = 1
    bpass.bpassprm[4]  =-1
    bpass.bpassprm[5]  = 0
    bpass.bpassprm[10] = 3
    bpass.go()

##############################################################################
# Run testfringe to find best source/antenna
#
def testfringe(indata, refant, flag, logfile):
    fringe             = AIPSTask('FRING')
    fringe.indata      = indata
    fringe.refant      = refant
    fringe.docal       = 1
    fringe.calsour[1:] = ''
    fringe.solint      = 6
    fringe.dparm[8]    = 1
    fringe.dparm[4]    = 0
    fringe.snver       = 2
    a = indata.sources
    bt = check_geo(indata)
    nb = len(bt)-1
    if (flag==1 and nb>2):
        time1 = bt[1]
        time2 = bt[nb-1]
        (day1,hour1,min1,sec1)=time_to_hhmmss(time1)
        (day2,hour2,min2,sec2)=time_to_hhmmss(time2)
        fringe.timerang[1:] = [day1, hour1, min1, sec1, day2, hour2, min2, sec2]
        mprint('#################################',logfile)
        mprint('TIMERANGE: '+str(day1)+' '+str(hour1)+' '+str(min1)+
               ' '+str(sec1)+' '+str(day2)+' '+str(hour2)+' '
               +str(min2)+' '+str(sec2),logfile)
        mprint('#################################',logfile)
    else:
        mprint('#################################',logfile)
        mprint('Whole Timerange',logfile)
        mprint('#################################',logfile)
    fringe()

##############################################################################
# Manual phasecal
#
def geo_man_pcal(indata, refant, logfile):
    fringe          = AIPSTask('FRING')
    fringe.indata   = indata
    fringe.refant   = refant
    fringe.docal    = 1
    fringe.solint   = 6
    fringe.dparm[8] = 1
    fringe.dparm[4] = 0
    fringe.snver    = 0

    (source,timerange)=get_best_scan(indata, logfile, 'geoqual.dat', 1)

    # Delete SN table from test fringe.
    check_sncl(indata, 0, 3, logfile)

    fringe.calsour[1]  = source
    fringe.timerang[1:] = timerange
    fringe()

    sn=indata.table('AIPS SN', 1)
    mprint('###########################################',logfile)
    mprint('Found solutions for '+str(len(sn))+' of '
           +str(len(indata.antennas))+' antennas.',logfile)
    mprint('###########################################',logfile)

##############################################################################
# Manual phasecal
#
def geo_man_pcal_sx(indata_s, indata_x, refant, logfile):
    fringe          = AIPSTask('FRING')
    fringe.indata   = indata_x
    fringe.refant   = refant
    fringe.docal    = 1
    fringe.solint   = 6
    fringe.dparm[8] = 1
    fringe.dparm[4] = 0
    fringe.snver    = 0

    (source,timerange)=get_best_scan(indata_x, logfile, 'geoqual.dat', 1)

    # Delete SN table from test fringe.
    check_sncl(indata_x, 0, 3, logfile)

    fringe.calsour[1]  = source
    fringe.timerang[1:] = timerange
    fringe()

    fringe.indata   = indata_s
    check_sncl(indata_s, 0, 3, logfile)
    fringe()

    sn_s=indata_s.table('AIPS SN', 1)
    sn_x=indata_x.table('AIPS SN', 1)
    mprint('###########################################',logfile)
    mprint('Found solutions for '+str(len(sn_s))+' of '
           +str(len(indata_s.antennas))+' antennas (S-band).',logfile)
    mprint('Found solutions for '+str(len(sn_x))+' of '
           +str(len(indata_x.antennas))+' antennas (X-band).',logfile)
    mprint('###########################################',logfile)

##############################################################################
#
def get_sources(indata):
    su_table = indata.table('AIPS SU', 0)
    max_source = 0

    for i in su_table:
        if i['id__no']>max_source:
            max_source=i['id__no']

    sources=[]
    for i in range(max_source):
        sources.append([])

    for i in su_table:
        sources[i['id__no']-1]=i['source']
    return sources

##############################################################################
# Find best scan for manual phasecal
#
def get_best_scan(indata, logfile, qualfile, do_write):
    sn_table = indata.table('AIPS SN', 0)
    naxis    = indata.header['naxis']
    sources  = get_sources(indata)

    t   = [0]
    qf  = [0]
    snr = []
    tr  = []
    sid = []

    n       = 0
    max_sol = 0

    sid.append(sn_table[0]['source_id'])
    tr.append(get_timerange_tab(indata,'AIPS SN',0))
    snr.append([])

    if naxis[3]>1:
        for j in range(naxis[3]):
            if isinde(sn_table[0]['delay_1'][j])==False:
                t[n]+=1
                qf[n]=qf[n]+1./sn_table[0]['weight_1'][j]**2
                snr[n].append(sn_table[0]['weight_1'][j])

        for i in range(1,len(sn_table)):
            if sn_table[i]['time']==sn_table[i-1]['time']:
                for j in range(naxis[3]):
                    if isinde(sn_table[i]['delay_1'][j])==False:
                        t[n]+=1
                        qf[n]=qf[n]+1./sn_table[i]['weight_1'][j]**2
                        snr[n].append(sn_table[i]['weight_1'][j])
                if t[n]>max_sol:
                    max_sol=t[n]
                    id=n
            else:
                t.append(0)
                qf.append(0)
                snr.append([])
                n+=1
                sid.append(sn_table[i]['source_id'])
                tr.append(get_timerange_tab(indata,'AIPS SN',i))
                for j in range(naxis[3]):
                    if isinde(sn_table[i]['delay_1'][j])==False:
                        t[n]+=1
                        qf[n]=qf[n]+1./sn_table[i]['weight_1'][j]**2
                        snr[n].append(sn_table[i]['weight_1'][j])
                if t[n]>max_sol:
                    max_sol=t[n]
                    id=n

    elif naxis[3]==1:
        if isinde(sn_table[0]['delay_1'])==False:
            t[n]+=1
            qf[n]=qf[n]+1./sn_table[0]['weight_1']**2
            snr[n].append(sn_table[0]['weight_1'])

        for i in range(1,len(sn_table)):
            if sn_table[i]['time']==sn_table[i-1]['time']:
                if isinde(sn_table[i]['delay_1'])==False:
                    t[n]+=1
                    qf[n]=qf[n]+1./sn_table[i]['weight_1']**2
                    snr[n].append(sn_table[i]['weight_1'])
                if t[n]>max_sol:
                    max_sol=t[n]
                    id=n
            else:
                t.append(0)
                qf.append(0)
                snr.append([])
                n+=1
                sid.append(sn_table[i]['source_id'])
                tr.append(get_timerange_tab(indata,'AIPS SN',i))
                if isinde(sn_table[i]['delay_1'])==False:
                    t[n]+=1
                    qf[n]=qf[n]+1./sn_table[i]['weight_1']**2
                    snr[n].append(sn_table[i]['weight_1'])
                if t[n]>max_sol:
                    max_sol=t[n]
                    id=n

    if do_write==1:
        file = './'+qualfile
        f = open(file,'w')
        for i in range(len(t)):
            f.writelines(' '+sources[sid[i]-1]+' Sol: '+str(t[i])+' QF: '+str(round(max(1/(qf[i]-0.00001),0),3))+' '+str(tr[i])+'\n')
        f.close()


    scan=0
    good_scans=[]
    bad_scans=[]
    bad_sources=[]
    for i in range(len(t)):
        if t[i]==max_sol:
            good_scans.append(i)
        elif t[i]<max_sol*0.4:
            bad_scans.append(i)

    scan=good_scans[0]
    source=sources[sid[0]-1]
    timerange=tr[0]

    for i in good_scans:
        if qf[i]<=qf[scan]:
            scan=i
            source=sources[sid[i]-1]
            timerange=tr[i]

    for i in range(len(bad_scans)):
        k=bad_scans[i]
        bad_sources.append(sources[sid[k]-1])
    mprint('#################################',logfile)
    mprint('Bad sources: '+str(bad_sources),logfile)
    mprint('#################################',logfile)

    mprint('#################################',logfile)
    mprint('Manual phase-cal on: '+source,logfile)
    mprint('#################################',logfile)

    max_sol=naxis[3]*naxis[1]*len(indata.antennas)

    mprint('#################################',logfile)
    mprint('TIMERANGE: '+str(timerange),logfile)
    mprint('Max number of solutions: '+str(max_sol),logfile)
    mprint('#################################',logfile)

    return source, timerange

##############################################################################
#
def runlistr(indata):
    listr=AIPSTask('LISTR')
    listr.indata   = indata
    listr.optype   = 'SCAN'
    listr.docrt    = -1
    if os.path.exists(indata.name+'.LST'):
        os.popen('rm '+indata.name+'.LST')
    listr.outprint = 'PWD:'+indata.name.strip()+'.LST'
    listr()

##############################################################################
# Make control_file.inp
#
def make_control_file(indata,refant,max_ant):
    (year, month, day)=get_observation_year_month_day(indata)
    f = open('./geoblock/control_file.inp','w')
    f.writelines('    1.0000                          '+
    '             ! data format (0=old; 1=new)\n'+
    '   30.0000                              ! number interations\n'+
    '    0.5000                              ! loop gain\n'+
    '    0.0000                              ! debug print out level\n')

    if month<10:
        s_month='0'+str(month)
    else:
        s_month=str(month)
    if day<10:
        s_day='0'+str(day)
    else:
        s_day=str(day)

    f.writelines(str(year)+' '+s_month+' '+s_day+
            '                              ! yyyy mm dd\n')
    frq=get_center_freq(indata)/1e9
    f.writelines(str(frq)+'d9                              ! obs freq (Hz)\n')
    f.writelines('    3.0        0.003        1.0         '+
                 '! delay err, rate err, re-weight flag\n')
    f.writelines('    0.0                 '+
                 '               ! Reference UT (hhmmss) for clock\n')
    bt=check_geo(indata)
    print 'Number of baselines is '+str(bt)
    for i in range(len(bt)-1):
        f.writelines('    '+str(bt[i])+'       '+str(bt[i+1])+'              '+
                 '       ! UT time range for block '+str(i+1)+'\n')
    if len(bt)<11:
        for i in range(len(bt)-1,10):
            f.writelines('    0.00       0.00              '+
                      '       ! UT time range for block '+str(i+1)+'\n')

    ant=get_ant(indata)
    nan=len(ant)

    co=1


    for i in range(max_ant):
        if i+1 in ant:
            if (i+1 in ant and i+1!=refant and max_ant<30):
                f.writelines('    0.0        1.0                      ! A'+str(i+1)+' C(cm)\n')
            else:
                f.writelines('    0.0        0.0                      ! A'+str(i+1)+' C(cm)\n')
        else:
            f.writelines('    0.0        0.0                      ! A'+str(i+1)+' C(cm)\n')

    for i in range(max_ant):
        if i+1 in ant:
            if (i+1 in ant and i+1!=refant and co>0 and max_ant<30):
                f.writelines('    0.0        1.0                      ! A'+str(i+1)+' Rcm/h\n')
            else:
                f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Rcm/h\n')
        else:
            f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Rcm/h\n')

    for i in range(max_ant):
        if i+1 in ant:
            if (i+1 in ant and i+1!=refant and co>1 and max_ant<30):
                f.writelines('    0.0        1.0                      ! A'+str(i+1)+' Ac/h2\n')
            else:
                f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Ac/h2\n')
        else:
            f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Ac/h2\n')



    for i in range(len(bt)-1):
        for j in range(max_ant):
            if j+1 in ant:
                if j+1 in ant:
                    f.writelines('    0.0        1.0                      '+
                                 '! A'+str(j+1)+' B'+str(i+1)+' cm\n')
                else:
                    f.writelines('    0.0        0.0                      '+
                                 '! A'+str(j+1)+' B'+str(i+1)+' cm\n')
            else:
                f.writelines('    0.0        0.0                      '+
                             '! A'+str(j+1)+' B'+str(i+1)+' cm\n')

    if len(bt)<11:
        for i in range(len(bt)-1,10):
            for j in range(max_ant):
                f.writelines('    0.0        0.0                      '+
                             '! A'+str(j+1)+' B'+str(i+1)+' cm\n')
    f.close()

##############################################################################
# Make control_file.inp
#
def make_control_file_sx(indata,pr_data,refant,max_ant):
    (year, month, day)=get_observation_year_month_day(indata)
    f = open('./geoblock/control_file_tropos.inp','w')
    f.writelines('    1.0000                          '+
    '             ! data format (0=old; 1=new)\n'+
    '   30.0000                              ! number interations\n'+
    '    0.5000                              ! loop gain\n'+
    '    0.0000                              ! debug print out level\n')

    if month<10:
        s_month='0'+str(month)
    else:
        s_month=str(month)
    if day<10:
        s_day='0'+str(day)
    else:
        s_day=str(day)

    f.writelines(str(year)+' '+s_month+' '+s_day+
            '                              ! yyyy mm dd\n')
    frq=get_center_freq(indata)/1e9
    f.writelines(str(frq)+'d9                              ! obs freq (Hz)\n')
    f.writelines('    3.0        0.003        1.0         '+
                 '! delay err, rate err, re-weight flag\n')
    f.writelines('    0.0                 '+
                 '               ! Reference UT (hhmmss) for clock\n')
    bt=check_geo(indata)
    for i in range(len(bt)-1):
        f.writelines('    '+str(bt[i])+'       '+str(bt[i+1])+'              '+
                 '       ! UT time range for block '+str(i+1)+'\n')
    if len(bt)<11:
        for i in range(len(bt)-1,10):
            f.writelines('    0.00       0.00              '+
                      '       ! UT time range for block '+str(i+1)+'\n')

    ant=get_ant(indata)
    nan=len(ant)

    co=1

    for i in range(max_ant):
        if i+1 in ant:
            if (i+1 in ant and i+1!=refant and max_ant<30):
                f.writelines('    0.0        1.0                      ! A'+str(i+1)+' C(cm)\n')
            else:
                f.writelines('    0.0        0.0                      ! A'+str(i+1)+' C(cm)\n')
        else:
            f.writelines('    0.0        0.0                      ! A'+str(i+1)+' C(cm)\n')

    for i in range(max_ant):
        if i+1 in ant:
            if (i+1 in ant and i+1!=refant and co>0 and max_ant<30):
                f.writelines('    0.0        1.0                      ! A'+str(i+1)+' Rcm/h\n')
            else:
                f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Rcm/h\n')
        else:
            f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Rcm/h\n')

    for i in range(max_ant):
        if i+1 in ant:
            if (i+1 in ant and i+1!=refant and co>1 and max_ant<30):
                f.writelines('    0.0        1.0                      ! A'+str(i+1)+' Ac/h2\n')
            else:
                f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Ac/h2\n')
        else:
            f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Ac/h2\n')



    for i in range(len(bt)-1):
        for j in range(max_ant):
            if j+1 in ant:
                if j+1 in ant:
                    f.writelines('    0.0        1.0                      '+
                                 '! A'+str(j+1)+' B'+str(i+1)+' cm\n')
                else:
                    f.writelines('    0.0        0.0                      '+
                                 '! A'+str(j+1)+' B'+str(i+1)+' cm\n')
            else:
                f.writelines('    0.0        0.0                      '+
                             '! A'+str(j+1)+' B'+str(i+1)+' cm\n')

    if len(bt)<11:
        for i in range(len(bt)-1,10):
            for j in range(max_ant):
                f.writelines('    0.0        0.0                      '+
                             '! A'+str(j+1)+' B'+str(i+1)+' cm\n')
    f.close()


    f = open('./geoblock/control_file_ionos.inp','w')
    f.writelines('    1.0000                          '+
    '             ! data format (0=old; 1=new)\n'+
    '   50.0000                              ! number interations\n'+
    '    0.5000                              ! loop gain\n'+
    '    0.0000                              ! debug print out level\n')

    if month<10:
        s_month='0'+str(month)
    else:
        s_month=str(month)
    if day<10:
        s_day='0'+str(day)
    else:
        s_day=str(day)

    f.writelines(str(year)+' '+s_month+' '+s_day+
            '                              ! yyyy mm dd\n')
    frq=get_center_freq(pr_data)/1e9
    f.writelines(str(frq)+'d9                              ! obs freq (Hz)\n')
    f.writelines('    2.0        0.999        1.0         '+
                 '! delay err, rate err, re-weight flag\n')
    f.writelines('    0.0                 '+
                 '               ! Reference UT (hhmmss) for clock\n')
    bt=check_geo(indata)
    for i in range(len(bt)-1):
        f.writelines('    '+str(bt[i])+'       '+str(bt[i+1])+'              '+
                 '       ! UT time range for block '+str(i+1)+'\n')
    if len(bt)<6:
        for i in range(len(bt)-1,5):
            f.writelines('    0.00       0.00              '+
                      '       ! UT time range for block '+str(i+1)+'\n')

    ant=get_ant(indata)
    nan=len(ant)

    co=1

    for i in range(max_ant):
        if i+1 in ant:
            if (i+1 in ant and i+1!=refant):
                f.writelines('    0.0        0.0                      ! A'+str(i+1)+' C(cm)\n')
            else:
                f.writelines('    0.0        0.0                      ! A'+str(i+1)+' C(cm)\n')
        else:
            f.writelines('    0.0        0.0                      ! A'+str(i+1)+' C(cm)\n')

    for i in range(max_ant):
        if i+1 in ant:
            if (i+1 in ant and i+1!=refant and co>0):
                f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Rcm/h\n')
            else:
                f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Rcm/h\n')
        else:
            f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Rcm/h\n')

    for i in range(max_ant):
        if i+1 in ant:
            if (i+1 in ant and i+1!=refant and co>1):
                f.writelines('    0.0        1.0                      ! A'+str(i+1)+' Ac/h2\n')
            else:
                f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Ac/h2\n')
        else:
            f.writelines('    0.0        0.0                      ! A'+str(i+1)+' Ac/h2\n')



    for i in range(len(bt)-1):
        for j in range(max_ant):
            if j+1 in ant:
                if j+1 in ant:
                    f.writelines('    0.0        1.0                      '+
                                 '! A'+str(j+1)+' B'+str(i+1)+' cm\n')
                else:
                    f.writelines('    0.0        0.0                      '+
                                 '! A'+str(j+1)+' B'+str(i+1)+' cm\n')
            else:
                f.writelines('    0.0        0.0                      '+
                             '! A'+str(j+1)+' B'+str(i+1)+' cm\n')

    if len(bt)<11:
        for i in range(len(bt)-1,10):
            for j in range(max_ant):
                f.writelines('    0.0        0.0                      '+
                             '! A'+str(j+1)+' B'+str(i+1)+' cm\n')
    f.close()

##############################################################################
# Make station_file.inp
#
def make_station_file(indata):
    ant=get_ant(indata)

    f = open('./geoblock/station_file.inp','w')
    for i in ant:
        if i==min(ant):
            f.writelines(str(ant[i])+' '+str(i)+'   ! station_file.inp for '
                         +indata.header.observer+' on '+indata.header.date_obs+'\n')
        else:
            if 'Y' in ant[i]:
                if not 'YG' in ant[i]:
                    pass
                else:
                    f.writelines(str(ant[i])+' '+str(i)+'\n')
            else: f.writelines(str(ant[i])+' '+str(i)+'\n')
    f.close()

##############################################################################
# Make calibrator_file.inp
#
def make_calibrator_file(indata):
    (year, month, day)=get_observation_year_month_day(indata)
    monthlist=['JAN','FEB','MAR','APR','MAY','JUN',
               'JUL','AUG','SEP','OCT','NOV','DEC']

    if day<10:
        strday='0'+str(day)
    else:
        strday=str(day)

    namma=str(year)+monthlist[int(month)-1]+strday

    f = open('./geoblock/calibrator_file.inp','w')
    f.writelines(namma+'_RATE_MDEL.DAT                   '+
                 '                   ! calibrator data file name\n')
    f.writelines(namma+'_SU_TABLE.PRTAB                  '+
                 '                   ! SU table file name\n')
    f.close()

##############################################################################
# Make calibrator_file.inp
#
def make_calibrator_file_sx(indata):
    (year, month, day)=get_observation_year_month_day(indata)
    monthlist=['JAN','FEB','MAR','APR','MAY','JUN',
               'JUL','AUG','SEP','OCT','NOV','DEC']

    if day<10:
        strday='0'+str(day)
    else:
        strday=str(day)

    namma=str(year)+monthlist[int(month)-1]+strday

    f = open('./geoblock/calibrator_file_tropos.inp','w')
    f.writelines(namma+'-tropos.dat                     '+
                 '                   ! calibrator data file name\n')
    f.writelines(namma+'_SU_TABLE.PRTAB                  '+
                 '                   ! SU table file name\n')
    f.close()

    f = open('./geoblock/calibrator_file_ionos.inp','w')
    f.writelines(namma+'-ionos.dat                      '+
                 '                   ! calibrator data file name\n')
    f.writelines(namma+'_SU_TABLE.PRTAB                  '+
                 '                   ! SU table file name\n')
    f.close()

##############################################################################
# Print out SN table
#
def runprtsn(indata):
    prtab           = AIPSTask('PRTAB')
    prtab.indata    = indata
    prtab.inext     = 'SN'
    prtab.invers    = 0
    prtab.docrt     = -1
    prtab.ndig      = 1
    prtab.box[1][1] = 1
    prtab.box[1][2] = 3
    prtab.box[1][3] = 4
    prtab.box[1][4] = 9
    if check_sn_ver(indata)>15:
        prtab.box[2][1] = 15
        prtab.box[2][2] = 17
    else:
        prtab.box[2][1] = 13
        prtab.box[2][2] = 15
    prtab.dohms     = -1
    prtab.doflag    = 0
    (year, month, day)=get_observation_year_month_day(indata)
    monthlist=['JAN','FEB','MAR','APR','MAY','JUN',
               'JUL','AUG','SEP','OCT','NOV','DEC']

    if day<10:
        strday='0'+str(day)
    else:
        strday=str(day)

    namma=str(year)+monthlist[int(month)-1]+strday

    prtab.outprint='PWD:'+namma+'_RATE_MDEL.DAT'
    prtab()

##############################################################################
# Print out SN table
#
def runprtsn_sx(indata_s, indata_x, indata):
    prtab           = AIPSTask('PRTAB')
    prtab.indata    = indata_s
    prtab.inext     = 'SN'
    prtab.invers    = 0
    prtab.docrt     = -1
    prtab.ndig      = 2
    prtab.box[1][1] = 1
    prtab.box[1][2] = 3
    prtab.box[1][3] = 4
    prtab.box[1][4] = 9
    if check_sn_ver(indata_s)>15:
        prtab.box[2][1] = 15
        prtab.box[2][2] = 17
    else:
        prtab.box[2][1] = 13
        prtab.box[2][2] = 15
    prtab.dohms     = -1
    prtab.doflag    = 0
    (year, month, day)=get_observation_year_month_day(indata_s)
    monthlist=['JAN','FEB','MAR','APR','MAY','JUN',
               'JUL','AUG','SEP','OCT','NOV','DEC']

    if day<10:
        strday='0'+str(day)
    else:
        strday=str(day)

    namma=str(year)+monthlist[int(month)-1]+strday

    prtab.outprint='PWD:'+namma+'_low.dat'
    prtab()

    prtab.indata    = indata_x
    prtab.outprint='PWD:'+namma+'_high.dat'
    prtab()

    f_low      = get_center_freq(indata_s)
    f_high     = get_center_freq(indata_x)
    if indata.exists():
        f_tar      = get_center_freq(indata)
    else:
        print 'Target frequency unknonw. Using 6.67 GHz.'
        f_tar = 6.67e9


    f = open('./geoblock/diff_files.inp','w')
    f.writelines(namma+'_low.dat                   '+
                 '                  ! Low frequency data\n')
    f.writelines(namma+'_high.dat                   '+
                 '                 ! High frequency data\n')
    f.writelines(str(f_low/1e9)+'                  '+
                 '                               ! Low frequency\n')
    f.writelines(str(f_high/1e9)+'                  '+
                 '                               ! High frequency\n')
    f.writelines(str(f_tar/1e9)+'                  '+
                 '                             ! Target frequency\n')
    f.close()


##############################################################################
# Print out SU table
#
def runprtsu(indata):
    prtab           = AIPSTask('PRTAB')
    prtab.indata    = indata
    prtab.inext     = 'SU'
    prtab.invers    = 0
    prtab.docrt     = -1
    prtab.box[1][1] = 1
    prtab.box[1][2] = 2
    prtab.box[1][3] = 11
    prtab.box[1][4] = 12
    prtab.box[2][1] = 13
    prtab.box[2][2] = 14
    prtab.box[2][3] = 15
    prtab.dohms     = -1
    prtab.ndig      = 4
    (year, month, day)=get_observation_year_month_day(indata)
    monthlist=['JAN','FEB','MAR','APR','MAY','JUN',
               'JUL','AUG','SEP','OCT','NOV','DEC']

    if day<10:
        strday='0'+str(day)
    else:
        strday=str(day)

    namma=str(year)+monthlist[int(month)-1]+strday

    prtab.outprint='PWD:'+namma+'_SU_TABLE.PRTAB'
    prtab()

##############################################################################
# Run DELZN
#
def run_delzn(indata):
    indata.zap_table('PL', -1)
    if os.path.exists('DELZN.FITS'):
        os.popen('rm DELZN.FITS')

    ant=get_ant(indata)

    delzn           = AIPSTask('DELZN')
    delzn.indata    = indata
    delzn.aparm[1:] = [0,2,2,0,0,0,0,0,1,0]
    delzn.nplots    = len(ant)-1
    delzn.outfile   = 'PWD:DELZN.FITS'
    delzn.dotv      = -1
    delzn()

    if os.path.exists('delzn.ps'):
        os.popen('rm delzn.ps')
    lwpla         = AIPSTask('LWPLA')
    lwpla.indata  = indata
    lwpla.inver   = 0
    lwpla.outfile = 'PWD:delzn.ps'
    lwpla()

    indata.zap_table('PL', -1)

##############################################################################
# Run CLCAL
#
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

def runclcal2(indata,snver,gainver,gainuse,interpol,dobtween,refant, antenna, cals):
    clcal          = AIPSTask('CLCAL')
    clcal.indata   = indata
    clcal.refant   = refant
    clcal.antennas[1:] = antenna
    clcal.calsour[1:]   = cals
    clcal.snver    = snver
    clcal.inver    = 0
    clcal.gainver  = gainver
    clcal.gainuse  = gainuse
    clcal.interpol = interpol
    clcal.dobtween = dobtween
    clcal()

##############################################################################
# Run QUACK
#
def runquack(indata, antennas, time):
    quack              = AIPSTask('QUACK')
    quack.indata       = indata
    quack.antennas[1:] = antennas
    quack.opcode       = 'ENDB'
    quack.aparm[1:]    =  [0,time,0]
    quack()

##############################################################################
# Fringe fit on all geosources
#
def fringegeo(indata, refant):
    if indata.header.telescop=='EVLA':
        print 'Telescope =  EVLA'
        d5 = 3
        d3 = 0
    else:
        d5 = 1
        d3 = 0

    fringe             = AIPSTask('FRING')
    fringe.indata      = indata
    fringe.refant      = refant
    fringe.docal       = 1
    fringe.calsour[1:] = ''
    fringe.solint      = 6
    fringe.weightit    = 3
    fringe.aparm[1:]   = [2, 0, d3, 0, d5, 0, 3]
    fringe.dparm[1:]   = [1, 30, 100, 0]
    fringe.dparm[4]    = 0
    fringe.dparm[8]    = 0
    fringe.snver       = 2
    fringe()

##############################################################################
# Write out qualities of geosources
#
def get_geosource_stats(indata):
    (year, month, day)=get_observation_year_month_day(indata)
    frq      = get_center_freq(indata)
    sn_table = indata.table('AIPS SN', 0)
    naxis    = indata.header['naxis']
    obs      = indata.header['observer']
    sources  = get_sources(indata)

    num_sour=len(indata.sources)
    snr=[]
    for i in range(0,len(sources)+1):
        snr.append([])

    for i in range(1,len(sn_table)):
        for j in range(naxis[3]):
                if isinde(sn_table[i]['weight_1'][j])==False and sn_table[i]['weight_1'][j] > 6:
                    snr[sn_table[i]['source_id']].append(sn_table[i]['weight_1'][j])

    file = './geoblock/'+obs+'_geostat.dat'
    f = open(file,'w')
    f.writelines('!Source   #Det    Max      Min   Avgerage  Median on '+str(year)+'-'+str(month)+'-'+str(day)+' at '+str(round(frq,1))+' GHz \n')
    for i in range(1,len(sources)):
        if len(snr[i])>0 and sources[i]!=[]:
            outprint= sources[i]+'  %4d %6.1f   %6.1f   %6.1f   %6.1f\n' % (len(snr[i]), round(max(snr[i]),1), round(min(snr[i]),1), round(average(snr[i]),1), round(median(snr[i]),1))
            f.writelines(outprint)
        elif sources[i]!=[]:
            outprint= sources[i]+'  %4d %6.1f   %6.1f   %6.1f   %6.1f\n' % (0, 0, 0, 0, 0)
            f.writelines(outprint)
    f.close()

##############################################################################
#
def runsnsmo(indata, inver, outver, refant):
    snsmo           = AIPSTask('SNSMO')
    snsmo.indata    = indata
    snsmo.refant    = refant
    snsmo.inver     = inver
    snsmo.outver    = outver
    snsmo.bparm[1:] = [0, 0, 1, 1, 1]
    snsmo.smotype   = 'VLBI'
    snsmo()

##############################################################################
#
def restore_su(indata, logfile):
    tasav_data=AIPSUVData(indata.name,'TASAV'+str(i),int(indata.disk),1)
    if tasav_data.exists():
        mprint('##########################################', logfile)
        mprint('TASAV file exists, restoring old SU table.', logfile)
        mprint('##########################################', logfile)
        while indata.table_highver('AIPS SU')>0:
            indata.zap_table('AIPS SU', 1)
        runtacop(tasav_data, indata, 'SU', 1, 1, 0)
    else:
        mprint('###############################################',logfile)
        mprint('No TASAV file. Restoring SU table not possible.',logfile)
        mprint('###############################################',logfile)

##############################################################################
#
def restore_fg(indata, logfile):
    tasav_data=AIPSUVData(indata.name,'TASAV'+str(i),int(indata.disk),1)
    if tasav_data.exists():
        mprint('##########################################', logfile)
        mprint('TASAV file exists, restoring old FG table.', logfile)
        mprint('##########################################', logfile)
        while indata.table_highver('AIPS FG')>0:
            indata.zap_table('AIPS FG', 0)
        runtacop(tasav_data, indata, 'FG', 1, 1, 0)
    else:
        mprint('###############################################',logfile)
        mprint('No TASAV file. Restoring SU table not possible.',logfile)
        mprint('###############################################',logfile)

##############################################################################
#
def runtasav(indata, i, logfile):
    tasav         = AIPSTask('TASAV')
    tasav.indata  = indata
    tasav.outna   = indata.name
    tasav.outcla  = 'TASAV'+str(i)
    tasav.outdisk = indata.disk
    tasav_data=AIPSUVData(indata.name,'TASAV'+str(i),int(indata.disk),1)
    if tasav_data.exists():
        mprint('TASAV file exists, do not need save tables', logfile)
    else:
        tasav()

##############################################################################
#
def runtacop(indata, outdata, inext, inver, outver, ncount):
    tacop         = AIPSTask('TACOP')
    tacop.indata  = indata
    tacop.outdata = outdata
    tacop.inext   = inext
    tacop.inver   = inver
    tacop.outver  = outver
    tacop.ncount  = ncount
    tacop()

##############################################################################
#
def get_time():
    t=range(6)
    t[0]=time.localtime()[0]
    t[0]=str(t[0])
    for i in range(1,6):
        t[i]=time.localtime()[i]
        if t[i]<10:
            t[i]='0'+str(t[i])
        else:
            t[i]=str(t[i])
    a=t[3]+':'+t[4]+':'+t[5]+' on '+t[0]+'/'+t[1]+'/'+t[2]
    return a

##############################################################################
# Plot results from geofit
#
def plot_baseline(indata,ant, inter_flag, ptype, doplot_flag, logfile):
    baselines=[]
    str_baselines=[]
    str_baselines2=[]
    data=[]
    start=100
    end=0

    ant_str=get_ant(indata)

    if ant<10:
        str_ant='0'+str(ant)
    else:
        str_ant=str(ant)

    for ant2 in ant_str:
        if ant2<10:
            str_ant2='0'+str(ant2)
        else:
            str_ant2=str(ant2)
        if ant2>ant:
            bl=str_ant+str_ant2
            str_bl=str(ant)+'-'+str(ant2)
            str_bl2=ant_str[ant]+'-'+ant_str[ant2]
        else:
            bl=str_ant2+str_ant
            str_bl=str(ant2)+'-'+str(ant)
            str_bl2=ant_str[ant2]+'-'+ant_str[ant]

        if ptype=='A':
            file='./geoblock/tropos/fit_geodetic_bl'+bl+'.dat'
        elif ptype=='I':
            file='./geoblock/ionos/fit_geodetic_bl'+bl+'.dat'
        if os.path.exists(file):
            f=open(file)
            content=f.read()
            if len(content)>99:
                newdata=loadtxt(file,usecols=(0,1,2,3,4,5,6),comments='!')
                data.append(newdata)
                baselines.append(ant2)
                str_baselines.append(str_bl)
                str_baselines2.append(str_bl2)

    num_baselines=len(baselines)

    f1=figure()

    if ptype=='A':
        toptext=indata.header['observer']+' observed on '+indata.header['date_obs']+' at '+str(round(get_center_freq(indata)/1e9,2))+' GHz'
    elif ptype=='I':
        toptext=indata.header['observer']+' observed on '+indata.header['date_obs']+' for target frequency'
    figtext(0.25,0.97,toptext)


    times2 = []
    block_times=[]
    for tmp in data:
        for tmp2 in tmp:
            times2.append(tmp2[0])
    times2.sort()
    start=min(times2)
    end=max(times2)

    n=1
    block_times.append([min(times2),0])
    for i in range(1,len(times2)):
        if times2[i]-times2[i-1]>1:
            block_times[n-1][1]=times2[i-1]
            block_times.append([times2[i],0])
            n+=1
    block_times[n-1][1]=max(times2)
    n_blocks=len(block_times)

    size=3
    print 'Baseline      rms_delay  rms_rate      data/block'
    print '               [nsec]     [mHz]'

    if len(baselines)>20:
        print 'More than 20 antennas, plotting only 15.'
        b=int(round(len(baselines)/2.))
        num_baselines=b
        num_baselines2=len(baselines)-b
        oneplot = False
    else:
        oneplot = True
        b=len(baselines)

    n=0
    m=0
    for i in baselines[0:b]:
        n+=1
        m+=1
        make_baseline_plot2(num_baselines, data, n, m, start, end, str_baselines, str_baselines2, n_blocks, block_times,f1)

    draw()
    if ptype=='A':
        savefig('delay-rate.ps',bbox_inches='tight')
    elif ptype=='I':
        savefig('delay-rate-ionos.ps',bbox_inches='tight')

    if oneplot==False:
        m=0
        f2=figure()
        toptext=indata.header['observer']+' observed on '+indata.header['date_obs']+' at '+str(round(get_center_freq(indata)/1e9,2))+' GHz'
        figtext(0.25,0.97,toptext)
        for i in baselines[b:]:
            n+=1
            m+=1
            make_baseline_plot2(num_baselines2, data, n, m, start, end, str_baselines, str_baselines2, n_blocks, block_times,f2)

        draw()
        if ptype=='A':
            savefig('delay-rate2.ps')
        elif ptype=='I':
            savefig('delay-rate2-ionos.ps')

    mprint('', logfile)
    mprint('Note: A small number of high delays or rates is usually no problem.',logfile)
    mprint('', logfile)


    if inter_flag==1:
        if int(matplotlib.__version__[0])==0:
            if int(matplotlib.__version__[2:4])<99:
                print 'Close figure to continue.'
                show()
                close()
            else:
                f1.show()
                raw_input('Press enter to close figure and continue. ')
                close()
        else:
            print 'Close figure to continue.'
            show()
            close()

##############################################################################
# Plot results from geofit
#
def make_baseline_plot(num_baselines, data, n, m, start, end, str_baselines, str_baselines2, n_blocks, block_times, f1):
        size=3
        ax=f1.add_subplot(num_baselines,2,2*m-1)
        high_delay=max(max(data[n-1][:,3]),1)
        low_delay=min(min(data[n-1][:,3]),-1)
        if max(data[n-1][:,3]) < 0.2:
            high_delay=max(max(data[n-1][:,3])*1.2,1.0)
            low_delay=min(min(data[n-1][:,3])*1.2,-1.0)
        ax.plot(data[n-1][:,0],data[n-1][:,1],'.b',ms=size) # Delay
        ax.plot(data[n-1][:,0],data[n-1][:,2],'og',ms=size) # Model
        ax.plot(data[n-1][:,0],data[n-1][:,3],'xr',ms=size) # Resid
        trang = end-start
        over = trang/10
        ax.plot([start-over, end+over], [0,0], 'k:')
        ax.set_xlim(start-over,end+over)
        if doplot_flag==1:
            ax.set_ylim(low_delay*1.1,high_delay*1.1)
        else:
            ax.set_ylim(-5,5)
        if m==1:
            title('Multiband delay [nsec]')
        if n<num_baselines:
             ax.xaxis.set_major_locator(NullLocator())
        if n==num_baselines:
             xlabel('UT [hours]')
        ax.text(0.2, 0.65, str_baselines2[n-1], transform=ax.transAxes)

        ax=f1.add_subplot(num_baselines,2,2*m)
        high_rate=max(max(data[n-1][:,6]),10)
        low_rate=min(min(data[n-1][:,6]),-10)
        ax.plot(data[n-1][:,0],data[n-1][:,4],'.b', ms=size)
        ax.plot(data[n-1][:,0],data[n-1][:,5],'.g', ms=size)
        ax.plot(data[n-1][:,0],data[n-1][:,6],'.r', ms=size)
        ax.plot([start-over, end+over], [0,0], 'k:')
        ax.set_xlim(start-over,end+over)
        if doplot_flag==1:
            ax.set_ylim(low_rate*1.1,high_rate*1.1)
        elif doplot_flag==2:
            ax.set_ylim(-5,5)
        else:
            ax.set_ylim(-10,10)
        if m==1:
            title('Rate [mHz]')
        if n<num_baselines:
             ax.xaxis.set_major_locator(NullLocator())
        if n==num_baselines:
             xlabel('UT [hours]')
        ax.text(0.2, 0.65, str_baselines2[n-1], transform=ax.transAxes)

        if max(data[n-1][:,3])>2 or min(data[n-1][:,3])<-2:
            delays=data[n-1][:,3]
            rms=delays.std()
            bad_delays=[]
            for entry in delays:
                if entry>3*rms or entry<-3*rms or entry<-2 or entry>2:
                    bad_delays.append(entry)
            if float(len(bad_delays))/float(len(delays))>0.3:
                warning1='\033[7;33;40m'+str(len(bad_delays))+'/'+str(len(delays))+' high delays\033[0m'
            else:
                warning1='\033[0m'+str(len(bad_delays))+'/'+str(len(delays))+' high delays\033[0m'
        else:
            warning1=''

        if max(data[n-1][:,6])>10 or min(data[n-1][:,6])<-10:
            rates=data[n-1][:,6]
            rms=rates.std()
            bad_rates=[]
            for entry in rates:
                if entry>3*rms or entry<-3*rms or entry<-10 or entry>10:
                    bad_rates.append(entry)
            if float(len(bad_rates))/float(len(rates))>0.3:
                warning2='\033[7;33;40m'+str(len(bad_rates))+'/'+str(len(rates))+' high rates\033[0m'
            else:
                warning2=str(len(bad_rates))+'/'+str(len(rates))+' high rates'
        else:
            warning2=''


        times=data[n-1][:,0]
        j=0

        blocks=[]
        for i in range(n_blocks):
            blocks.append(0)
            for tmp in times:
                if block_times[i][1]>=tmp>=block_times[i][0]:
                    blocks[i]+=1

        if 1 in blocks or 2 in blocks or 3 in blocks:
            warning0='\033[7;33;40m<3 scans\033[0m'
        else:
            warning0=''

        str_blocks=' '
        for k in range(n_blocks):
            if int(blocks[k])<10:
                str_blocks=str_blocks+'   '+str(int(blocks[k]))
            else:
                str_blocks=str_blocks+'  '+str(int(blocks[k]))

        mprint('%5s(%5s)    %4.2f      %4.2f    %18s %8s %17s %10s'% (str_baselines2[n-1],str_baselines[n-1], round(data[n-1][:,2].std(),3),round(data[n-1][:,6].std(),2),str_blocks,warning0,warning1,warning2), logfile)
        draw()

##############################################################################
# Plot results from geofit
#
def make_baseline_plot2(num_baselines, data, n, m, start, end, str_baselines, str_baselines2, n_blocks, block_times, f1):
        size=3
        ax=f1.add_subplot(num_baselines,1,m)
        high_delay=max(max(data[n-1][:,3]),1)
        low_delay=min(min(data[n-1][:,3]),-1)
        if max(data[n-1][:,3]) < 0.2:
            high_delay=max(max(data[n-1][:,3])*1.2,1.0)
            low_delay=min(min(data[n-1][:,3])*1.2,-1.0)
        ax.plot(data[n-1][:,0],data[n-1][:,1],'.b',ms=size) # Delay
        ax.plot(data[n-1][:,0],data[n-1][:,2],'og',ms=size) # Model
        ax.plot(data[n-1][:,0],data[n-1][:,3],'xr',ms=size) # Resid
        trang = end-start
        over = trang/8
        ax.plot([start-over, end+over], [0,0], 'k:')
        ax.set_xlim(start-over,end+over)
        if doplot_flag==1:
            ax.set_ylim(low_delay*1.1,high_delay*1.1)
        elif doplot_flag==2:
            ax.set_ylim(-5,5)
        else:
            ax.set_ylim(-10,10)
        if m==1:
            title('Multiband delay [nsec]')
        if m<num_baselines:
             ax.xaxis.set_major_locator(NullLocator())
        if m==num_baselines:
             xlabel('UT [hours]')
        ax.text(0.01, 0.60, str_baselines2[n-1], transform=ax.transAxes)

        if max(data[n-1][:,3])>2 or min(data[n-1][:,3])<-2:
            delays=data[n-1][:,3]
            rms=delays.std()
            bad_delays=[]
            for entry in delays:
                if entry>3*rms or entry<-3*rms or entry<-2 or entry>2:
                    bad_delays.append(entry)
            if float(len(bad_delays))/float(len(delays))>0.3:
                warning1='\033[7;33;40m'+str(len(bad_delays))+'/'+str(len(delays))+' high delays\033[0m'
            else:
                warning1='\033[0m'+str(len(bad_delays))+'/'+str(len(delays))+' high delays\033[0m'
        else:
            warning1=''

        if max(data[n-1][:,6])>10 or min(data[n-1][:,6])<-10:
            rates=data[n-1][:,6]
            rms=rates.std()
            bad_rates=[]
            for entry in rates:
                if entry>3*rms or entry<-3*rms or entry<-10 or entry>10:
                    bad_rates.append(entry)
            if float(len(bad_rates))/float(len(rates))>0.3:
                warning2='\033[7;33;40m'+str(len(bad_rates))+'/'+str(len(rates))+' high rates\033[0m'
            else:
                warning2=str(len(bad_rates))+'/'+str(len(rates))+' high rates'
        else:
            warning2=''


        times=data[n-1][:,0]
        j=0

        blocks=[]
        for i in range(n_blocks):
            blocks.append(0)
            for tmp in times:
                if block_times[i][1]>=tmp>=block_times[i][0]:
                    blocks[i]+=1

        if 1 in blocks or 2 in blocks or 3 in blocks:
            warning0='\033[7;33;40m<3 scans\033[0m'
        else:
            warning0=''

        str_blocks=' '
        for k in range(n_blocks):
            if int(blocks[k])<10:
                str_blocks=str_blocks+'   '+str(int(blocks[k]))
            else:
                str_blocks=str_blocks+'  '+str(int(blocks[k]))

        mprint('%5s(%5s)    %4.2f      %4.2f    %18s %8s %17s %10s'% (str_baselines2[n-1],str_baselines[n-1], round(data[n-1][:,2].std(),3),round(data[n-1][:,6].std(),2),str_blocks,warning0,warning1,warning2), logfile)
        draw()

##############################################################################
#
def checkatmos(inter_flag,logfile):
    file='ATMOS.FITS'
    data=loadtxt(file,skiprows=1)

    m=0
    for i in data:
        if i[5]!=0:
            m=m+1

    antennas=[]
    for i in range(len(data)):
        tmp=[]
        zero=[]
        nonzero=[]
        if data[i][5]==0:
            antennas.append(data[i][0])
            time_i=data[i][1]+(data[i][2]+data[i][3]/60.+data[i][4]/3600.)/24.
            n=0
            for j in range(len(data)):
                if data[j][0]==data[i][0]:
                    time=data[j][1]+(data[j][2]+data[j][3]/60.+data[j][4]/3600.)/24.
                    tmp.append([data[j][0],time,data[j][5]])
                    if data[j][5]==0 or data[j][5]==999:
                        zero.append(n)
                    else:
                        nonzero.append(n)
                    n=n+1

            if tmp!=[] and nonzero!=[]:
                 for k in zero:
                     if k==0 and time_i==tmp[k][1]:
                         data[i][5]=tmp[nonzero[0]][2]
                     elif k==n-1 and time_i==tmp[k][1]:
                         data[i][5]=tmp[nonzero[len(nonzero)-1]][2]
                     elif  time_i==tmp[k][1]:
                         data[i][5]=999

    n=0
    for i in data:
        if i[5]!=999 and i[5]!=0:
            n=n+1

    if m<len(data):
        mprint('Original ATMOS.FITS file has zero zenith delays.',logfile)
        mprint('Making new ATMOS file.',logfile)
        f=open('NEWATMOS.FITS', 'w')
        f.writelines('   '+str(n)+'\n')
        for i in data:
            if i[5]!=999 and i[5]!=0:
               line1 = ' %2d  %2d %2d %2d %4.1f' % (int(i[0]), int(i[1]), int(i[2]), int(i[3]), i[4])
               line2 = ' %8.3f   %8.3f    %8.5f    %8.5f' % (i[5],i[6],i[7], i[8])
               f.writelines(line1+line2+'\n')
        f.close()

        os.popen('mv ATMOS.FITS ATMOS.FITS.orig')
        os.popen('mv NEWATMOS.FITS ATMOS.FITS')
    else:
        mprint('ATMOS.FITS file ok.',logfile)

    #plotatmos(inter_flag, logfile)

##############################################################################
#
def make_name_atmos(indata):
    ant=get_ant(indata)
    file='ATMOS.FITS'
    data=loadtxt(file,skiprows=1)
    data2=[]
    n=int(max(data[:,0]))
    for j in range(1,n+1):
        for entry in data:
            if entry[0]==j:
                data2.append(entry)

    f=open('ATMOS_NAME.FITS', 'w')
    f.writelines('   '+str(len(data2))+'\n')
    for i in data2:
       line1 = '%2s  %2d %2d %2d %4.1f' % (ant[i[0]], int(i[1]), int(i[2]), int(i[3]), i[4])
       line2 = ' %8.3f   %8.3f    %8.5f    %8.5f' % (i[5],i[6],i[7], i[8])
       f.writelines(line1+line2+'\n')

    f.close()

##############################################################################
#
def make_name_ionos(indata):
    ant=get_ant(indata)
    file='IONOS.FITS'
    data=loadtxt(file,skiprows=1)
    data2=[]
    n=int(max(data[:,0]))
    for j in range(1,n+1):
        for entry in data:
            if entry[0]==j:
                data2.append(entry)

    f=open('IONOS_NAME.FITS', 'w')
    f.writelines('   '+str(len(data2))+'\n')
    for i in data2:
       line1 = '%2s  %2d %2d %2d %4.1f' % (ant[i[0]], int(i[1]), int(i[2]), int(i[3]), i[4])
       line2 = ' %8.3f   %8.3f    %8.5f    %8.5f' % (i[5],i[6],i[7], i[8])
       f.writelines(line1+line2+'\n')

    f.close()

#############################################################################
#
def plotatmos(inter_flag, logfile):
    file='ATMOS.FITS'
    data=loadtxt(file,skiprows=1)

    ant=[]
    data2=[]
    avg=[]
    rms=[]

    for i in range(int(max(data[:,0]))):
        ant.append([])
        data2.append([])
        avg.append(-1)
        rms.append(-1)
    for row in data:
        if (row[0] in ant) == False:
            ant[int(row[0])-1]=int(row[0])
        time=row[1]*24.+(row[2]+row[3]/60.+row[4]/3600.)
        data2[int(row[0])-1].append([int(row[0]), time,row[5], row[6], row[7], row[8]])

    max_ant=len(data2)
    num_ant=0
    for i in ant:
        if i!=[]:
            num_ant+=1

    fig=figure(0)
    n=0
    start=100
    end=0

    for entry in data2:
        n=n+1
        if entry !=[]:
            ant_data=array(entry)
            if start>min(ant_data[:,1]):
                start=min(ant_data[:,1])
            if end<max(ant_data[:,1]):
                end=max(ant_data[:,1])
            sum=0
            for i in range(len(ant_data)):
                sum=sum+ant_data[i][2]
            avg[n-1]=(sum/len(ant_data))
            rms[n-1]=(ant_data[:,2].std())

    n=0
    plot=0
    span=2
    mprint('',logfile)
    mprint('Plotting ATMOS.FITS file.',logfile)
    mprint('Antenna       rms [cm]', logfile)
    for entry in data2:
        n+=1
        if entry !=[]:
            plot+=1
            ant_data=array(entry)
            ax=fig.add_subplot(num_ant,1,plot)
            line = ' %4d  %12.3f ' % (int(ant_data[0][0]), round(rms[n-1],3))

            if (max(ant_data[:,2])<avg[n-1]+span) and (min(ant_data[:,2])>avg[n-1]-span):
                ax.plot(ant_data[:,1], ant_data[:,2], 'gx')
                ax.set_ylim(avg[n-1]-span,avg[n-1]+span)
                line2 = ''
            elif (max(ant_data[:,2])<avg[n-1]+2*span) and (min(ant_data[:,2])>avg[n-1]-2*span):
                ax.plot(ant_data[:,1], ant_data[:,2], 'bx')
                ax.set_ylim(avg[n-1]-2*span,avg[n-1]+2*span)
                line2 = '  (variations > '+str(int(span))+' cm from mean)'
            elif (max(ant_data[:,2])<avg[n-1]+3*span) and (min(ant_data[:,2])>avg[n-1]-3*span):
                ax.plot(ant_data[:,1], ant_data[:,2], 'rx')
                ax.set_ylim(avg[n-1]-3*span,avg[n-1]+3*span)
                line2 = '  (variations > '+str(int(2*span))+' cm from mean)'
            else:
                ax.plot(ant_data[:,1], ant_data[:,2], 'ro')
                ax.set_ylim(avg[n-1]-4*span,avg[n-1]+4*span)
                line2 = '  (variations > '+str(int(3*span))+' cm from mean)'
            ax.set_xlim(start-1.,end+1.)
            yticks([int(avg[n-1])-2*span,int(avg[n-1]),int(avg[n-1])+2*span])

            mprint(line+line2, logfile)

            ax.text(0.03, 0.60, str(int(ant_data[0][0])), transform=ax.transAxes)

            if n==1:
                title('ATMOS.FITS zenith delays [cm]')
            if n<max_ant:
                ax.xaxis.set_major_locator(NullLocator())
            if n==max_ant:
                 xlabel('UT [hours]')

    mprint('',logfile)
    mprint('Green *: variations < '+str(int(span))+' cm from mean',logfile)
    mprint('Blue  x: variations between '+str(int(span))+' and '+str(int(2*span))+' cm from mean',logfile)
    mprint('Red   x: variations between '+str(int(2*span))+' and '+str(int(3*span))+' cm from mean',logfile)
    mprint('Red   o: variations > '+str(int(3*span))+' cm from mean',logfile)
    mprint('',logfile)

    draw()
    savefig('atmos.ps',bbox_inches='tight')

    if inter_flag==1:
        if int(matplotlib.__version__[0])==0:
            if int(matplotlib.__version__[2:4])<99:
                print 'Close figure to continue.'
                show()
                close()
                cont=raw_input('Continue with current ATMOS.FITS file? (y/n) ')
                if cont=='n' or cont=='N':
                    sys.exit()
                close()
            else:
                fig.show()
                cont=raw_input('Continue with current ATMOS.FITS file? (y/n) ')
                if cont=='n' or cont=='N':
                    sys.exit()
                close()
        else:
            print 'Close figure to continue.'
            show()
            close()
            cont=raw_input('Continue with current ATMOS.FITS file? (y/n) ')
            if cont=='n' or cont=='N':
                sys.exit()
            close()

#############################################################################
#
def plotionos(inter_flag, logfile):
    file='IONOS.FITS'
    data=loadtxt(file,skiprows=1)

    ant=[]
    data2=[]
    avg=[]
    rms=[]

    for i in range(int(max(data[:,0]))):
        ant.append([])
        data2.append([])
        avg.append(-1)
        rms.append(-1)
    for row in data:
        if (row[0] in ant) == False:
            ant[int(row[0])-1]=int(row[0])
        time=row[1]*24.+(row[2]+row[3]/60.+row[4]/3600.)
        data2[int(row[0])-1].append([int(row[0]), time,row[5], row[6], row[7], row[8]])

    max_ant=len(data2)
    num_ant=0
    for i in ant:
        if i!=[]:
            num_ant+=1

    fig=figure(0)
    n=0
    start=100
    end=0

    for entry in data2:
        n=n+1
        if entry !=[]:
            ant_data=array(entry)
            if start>min(ant_data[:,1]):
                start=min(ant_data[:,1])
            if end<max(ant_data[:,1]):
                end=max(ant_data[:,1])
            sum=0
            for i in range(len(ant_data)):
                sum=sum+ant_data[i][2]
            avg[n-1]=(sum/len(ant_data))
            rms[n-1]=(ant_data[:,2].std())

    n=0
    plot=0
    span=2
    mprint('',logfile)
    mprint('Plotting IONOS.FITS file.',logfile)
    mprint('Antenna       rms [cm]', logfile)
    for entry in data2:
        n+=1
        if entry !=[]:
            plot+=1
            ant_data=array(entry)
            ax=fig.add_subplot(num_ant,1,plot)
            line = ' %4d  %12.3f ' % (int(ant_data[0][0]), round(rms[n-1],3))

            if (max(ant_data[:,2])<avg[n-1]+span) and (min(ant_data[:,2])>avg[n-1]-span):
                ax.plot(ant_data[:,1], ant_data[:,2], 'gx')
                ax.set_ylim(avg[n-1]-span,avg[n-1]+span)
                line2 = ''
            elif (max(ant_data[:,2])<avg[n-1]+2*span) and (min(ant_data[:,2])>avg[n-1]-2*span):
                ax.plot(ant_data[:,1], ant_data[:,2], 'bx')
                ax.set_ylim(avg[n-1]-2*span,avg[n-1]+2*span)
                line2 = '  (variations > '+str(int(span))+' cm from mean)'
            elif (max(ant_data[:,2])<avg[n-1]+3*span) and (min(ant_data[:,2])>avg[n-1]-3*span):
                ax.plot(ant_data[:,1], ant_data[:,2], 'rx')
                ax.set_ylim(avg[n-1]-3*span,avg[n-1]+3*span)
                line2 = '  (variations > '+str(int(2*span))+' cm from mean)'
            else:
                ax.plot(ant_data[:,1], ant_data[:,2], 'ro')
                ax.set_ylim(avg[n-1]-4*span,avg[n-1]+4*span)
                line2 = '  (variations > '+str(int(3*span))+' cm from mean)'
            ax.set_xlim(start-1.,end+1.)
            yticks([int(avg[n-1])-2*span,int(avg[n-1]),int(avg[n-1])+2*span])

            mprint(line+line2, logfile)

            ax.text(0.03, 0.60, str(int(ant_data[0][0])), transform=ax.transAxes)

            if n==1:
                title('IONOS.FITS zenith delays [cm]')
            if n<max_ant:
                ax.xaxis.set_major_locator(NullLocator())
            if n==max_ant:
                 xlabel('UT [hours]')

    mprint('',logfile)
    mprint('Green *: variations < '+str(int(span))+' cm from mean',logfile)
    mprint('Blue  x: variations between '+str(int(span))+' and '+str(int(2*span))+' cm from mean',logfile)
    mprint('Red   x: variations between '+str(int(2*span))+' and '+str(int(3*span))+' cm from mean',logfile)
    mprint('Red   o: variations > '+str(int(3*span))+' cm from mean',logfile)
    mprint('',logfile)

    draw()
    savefig('ionos.ps')

    if inter_flag==1:
        if int(matplotlib.__version__[0])==0:
            if int(matplotlib.__version__[2:4])<99:
                print 'Close figure to continue.'
                show()
                close()
                cont=raw_input('Continue with current IONOS.FITS file? (y/n) ')
                if cont=='n' or cont=='N':
                    sys.exit()
                close()
            else:
                fig.show()
                cont=raw_input('Continue with current IONOS.FITS file? (y/n) ')
                if cont=='n' or cont=='N':
                    sys.exit()
                close()
        else:
            print 'Close figure to continue.'
            show()
            close()
            cont=raw_input('Continue with current IONOS.FITS file? (y/n) ')
            if cont=='n' or cont=='N':
                sys.exit()
            close()

##############################################################################
#
def runatmos(indata, atmos_file):
    atmos              = AIPSTask('CLCOR')
    atmos.indata       = indata
    atmos.gainver      = 3
    atmos.gainuse      = 4
    atmos.clcorprm[1:] = [1,0]
    atmos.opcode       = 'ATMO'
    atmos.infile       = 'PWD:'+atmos_file
    atmos()

def runionos(indata, ionos_file):
    atmos              = AIPSTask('CLCOR')
    atmos.indata       = indata
    atmos.gainver      = 4
    atmos.gainuse      = 4
    atmos.clcorprm[1:] = [1,0]
    atmos.opcode       = 'IONO'
    atmos.infile       = 'PWD:'+ionos_file
    atmos()

##############################################################################
#
def runpang(indata):
    pang              = AIPSTask('CLCOR')
    pang.indata       = indata
    pang.gainver      = 4
    pang.gainuse      = 4
    pang.opcode       = 'PANG'
    pang.clcorprm[1:] = [1,0]
    #antennas          = []
    # the next bit is to deal with HB,KE,YG being linear
    #for row in indata.table('AN', 0):
    #    if row['mntsta']==0:
    #        if not row['anname'].replace(' ','') in ('HB','KE'):
    #            antennas.append(row['nosta'])
    #pang.antennas[1:] = antennas
    pang()

##############################################################################
#
def runpang2(indata):
    pang        = AIPSTask('CLCOR')
    pang.indata = indata
    pang.gainver = 3
    pang.gainuse = 4
    pang.opcode  = 'PANG'
    pang.clcorprm[1:] = [1,0]
    #antennas          = []
    # the next bit is to deal with HB,KE,YG being linear
    #for row in indata.table('AN', 0):
    #    if row['mntsta']==0:
    #        if not row['anname'].replace(' ','') in ('HB','KE'):
    #            antennas.append(row['nosta'])
    #pang.antennas[1:] = antennas
    pang()

##############################################################################
#
def runaccor(indata):
    accor           = AIPSTask('ACCOR')
    accor.indata    = indata
    accor.timer[1:] = [0]
    accor.solint    = 2
    accor()

##############################################################################
#
def run_setjy(indata, source, flux):
    setjy = AIPSTask('SETJY')
    setjy.indata = indata
    setjy.source[1:] = [source]
    setjy.zerosp[1:] = flux
    setjy.optype = ''
    setjy.bif    = 0
    setjy.eif    = 0
    setjy.optype = 'VCAL'
    print setjy.optype
    setjy()

##############################################################################
#
def runantab(indata,antabfile):
    antab = AIPSTask('ANTAB')
    antab.indata = indata
    antab.calin  = antabfile
    antab.tyver  = 1
    antab.gcver  = 1
    antab.go()

##############################################################################
#
def runtysmo(indata, tywin, maxdev):
    while indata.table_highver('AIPS TY') > 1:
        indata.zap_table('TY', 0)
    tysmo           = AIPSTask('TYSMO')
    tysmo.indata    = indata
    tysmo.dobtween  = 0
    tysmo.cparm[1:] = [tywin,0,0,0,0,maxdev,0]
    tysmo.inext     = 'TY'
    tysmo.inver     = 1
    tysmo.outver    = 2
    tysmo()

##############################################################################
#
def runapcal(indata, tyver, gcver, snver, dofit):
    apcal            = AIPSTask('APCAL')
    apcal.indata     = indata
    apcal.tyver      = tyver
    apcal.gcver      = gcver
    apcal.snver      = snver
    ant              = get_ant(indata)
    for i in ant:
        apcal.dofit[i]  = dofit
    apcal()
##############################################################################
# runapcal added by Lucas Hyland Nov 2017 for LBA
def runapcal_lba(indata, snver, inver, outver):
    sqrtsefd        = {'AT': 9.48, 'CD': 28.3, 'HH': 22.36, 'HO': 31.6, 'MP': 25.10,
                       'PA': 7.42, 'WA': 25.5, 'KE': 54.77, 'YG': 60.0, 'HB': 54.77}
    ant             = get_ant(indata)
    check_sncl(indata, snver, inver, logfile)
    for i in ant:
        clcor             = AIPSTask('CLCOR')
        clcor.opcode      = 'GAIN'
        clcor.indata      = indata
        clcor.gainver     = 0
        clcor.gainuse     = outver
        clcor.antenna[1:] = [i]
        clcor.clcorprm[1] = sqrtsefd[ant[i]]
        clcor.go()
    return

##############################################################################
#
def man_pcal(indata, refant, mp_source, mp_timera, debug, logfile, dpfour):

    if mp_source == ['']:
        mp_source = []
        for source in indata.sources:
            if source[0]=='F':
                mp_source.append(source)
    fringe            = AIPSTask('FRING')
    fringe.indata     = indata
    fringe.refant     = refant
    fringe.docal      = 1
    fringe.solint     = 6
    fringe.bchan      = 0
    fringe.echan      = 0
    fringe.aparm[1:]  =[2,0]
    fringe.dparm[2]   = 250
    fringe.dparm[3]   = 50
    fringe.dparm[4]   = dpfour
    fringe.dparm[8]   = 1
    fringe.snver      = 0
    fringe.calso[1:]  = mp_source
#    fringe.inputs()
    if mp_timera==0:
        fringe.timer[1:]=[0]
        fringe()
    else:
        fringe.timer[1:] = mp_timera
        fringe()

    qualfile = indata.name+'.'+indata.klass+'-qual.dat'
    (source,timerange)=get_best_scan(indata,logfile, qualfile, 1)

    # Delete SN table from test fringe.
    check_sncl(indata, 2, 6, logfile)
    fringe.calsour[1]  = source
    fringe.timerang[1:] = timerange
    fringe()

    sn=indata.table('AIPS SN', 0)
    mprint('###########################################',logfile)
    mprint('Found solutions for '+str(len(sn))+' of '
           +str(len(indata.antennas))+' antennas.',logfile)
    mprint('###########################################',logfile)
    return source, timerange

##############################################################################
#
def fringecal(indata, fr_image, nmaps, refant, calsource,solint,smodel, doband, bpver, dpfour):
    fringe             = AIPSTask('FRING')
    if fr_image.exists():
        fringe.in2data = fr_image
        mprint('################################################',logfile)
        mprint('Using input model '+fringe.in2name+'.'+fringe.in2class+'.'+str(int(fringe.in2seq))+' on diks '+str(int(fringe.in2disk)), logfile)
        mprint('################################################',logfile)
    elif smodel!=[1,0]:
        fringe.smodel[1:] = smodel
        mprint('################################################',logfile)
        mprint('Using SMODEL='+str(smodel)+' for fringe.',logfile)
        mprint('################################################',logfile)
    else:
        mprint('################################################',logfile)
        mprint('Using point source as imput model for fringe.',logfile)
        mprint('################################################',logfile)

    if doband==1:
        mprint('################################################',logfile)
        mprint('Applying bandpass table '+str(bpver), logfile)
        mprint('################################################',logfile)
    else:
        mprint('################################################',logfile)
        mprint('Applying no bandpass table ', logfile)
        mprint('################################################',logfile)

    fringe.indata      = indata
    fringe.refant      = refant
    fringe.docal       = 1
    fringe.calsour[1:] = [calsource]
    fringe.solint      = solint
    fringe.aparm[1:]   = [2, 0, 0, 0, 1]
    fringe.dparm[1:]   = [1, 20, 50, 0]
    fringe.dparm[4]    = dpfour
    fringe.dparm[8]    = 0
    fringe.nmaps       = nmaps
    fringe.snver       = 0
    fringe.doband      = int(doband)
    fringe.bpver       = int(bpver)

    fringe()

##############################################################################
#
def runimagr(indata,source,niter,cz,iz,docal,imna,antennas,uvwtfn,robust,beam,baselines=[0],timer=[0]):

    if imna=='':
        outname=source
    else:
        outname         = source[:11-len(imna)]+imna
    mprint('#########################################################',logfile)
    mprint('Imaging '+source+' with imsize='+str(iz)+', cellsize='+str(cz)+
           ' and '+str(niter)+' iterations. Using antennas='+str(antennas)+
           '.',logfile)
    mprint('#########################################################',logfile)
    imagr              = AIPSTask('IMAGR')
    imagr.indata       = indata
    imagr.docal        = docal
    imagr.sourc[1:]    = [source]
    imagr.uvwtfn       = uvwtfn
    imagr.robust       = robust
    imagr.timer[1:]    = timer
    imagr.bif          = 0
    imagr.eif          = 0
    imagr.nchav        = 16
    imagr.bchan        = 0
    imagr.echan        = 0
    imagr.cellsize[1:] = [cz,cz]
    imagr.imsize[1:]   = [iz,iz]
    imagr.outna        = outname
    imagr.niter        = niter
    imagr.gain         = float(20/niter)
    imagr.outdisk      = indata.disk
    imagr.dotv         = -1
    imagr.antennas[1:] = antennas
    imagr.baseline[1:] = baselines
    imagr.bmaj         = beam[0]
    imagr.bmin         = beam[1]
    imagr.bpa          = beam[2]
    imagr()


##############################################################################
#
def rungridimagr(indata,source,niter,cz,iz,docal,uvwtfn,robust,beam):
    imagr           = AIPSTask('IMAGR')
    imagr.indata    = indata
    imagr.docal     = docal
    imagr.sourc[1:] = [source]
    imagr.uvwtfn       = uvwtfn
    imagr.robust       = robust
    imagr.bif       = 0
    imagr.eif       = 0
    imagr.nchav     = 16
    imagr.uvtaper[1:] = [100000,100000]
    imagr.uvwtfn    = 'N'
    imagr.bchan     = 0
    imagr.echan     = 0
    imagr.cellsize[1:] = [cz,cz]
    imagr.imsize[1:]   = [iz,iz]
    imagr.outna     = source
    imagr.niter     = niter
    imagr.outdisk   = indata.disk
    imagr.dotv      = -1
    imagr.bmaj         = beam[0]
    imagr.bmin         = beam[1]
    imagr.bpa          = beam[2]
    imagr()

##############################################################################
#
def runmaimagr(indata,source,niter,cz,iz,channel,docal,imna,uvwtfn,robust,beam,
    baselines=[0],timer=[0]):

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

    imagr.indata    = indata
    imagr.docal     = docal
    imagr.sourc[1:] = [source]
    imagr.uvwtfn       = uvwtfn
    imagr.robust       = robust
    imagr.nchav     = 1
    imagr.bchan     = channel
    imagr.echan     = channel
    imagr.cellsize[1:] = [cz,cz]
    imagr.imsize[1:]   = [iz,iz]
    imagr.outna     = outname
    imagr.niter     = niter
    imagr.outdisk   = indata.disk
    imagr.dotv      = -1
    imagr.bmaj         = beam[0]
    imagr.bmin         = beam[1]
    imagr.bpa          = beam[2]
    imagr.timer[1:]    = timer
    imagr.baseline[1:] = baselines
    imagr()

##############################################################################
#
def runcube(indata,source,niter,cz,iz,bch, ech, docal, ant, uvwtfn, robust,beam):

    imagr = AIPSTask('IMAGR')

    if indata.header['naxis'][3]>1:
        imagr.bif = input('Enter IF: ')
        imagr.eif = imagr.bif
    else:
        imagr.eif = 1
        imagr.bif = 1

    imagr.indata    = indata
    imagr.docal     = docal
    imagr.sourc[1:] = [source]
    imagr.antennas[1:] = ant
    imagr.nchav     = 1
    imagr.bchan     = bch
    imagr.echan     = ech
    imagr.uvwtfn    = uvwtfn
    imagr.robust    = robust
    imagr.gain      = float(20/niter)
    imagr.cellsize[1:] = [cz,cz]
    imagr.imsize[1:]   = [iz,iz]
    imagr.outna     = source
    imagr.niter     = niter
    imagr.outdisk   = indata.disk
    imagr.dotv      = -1
    imagr.bmaj         = beam[0]
    imagr.bmin         = beam[1]
    imagr.bpa          = beam[2]
    imagr()

##############################################################################
#
def shift_pos(indata, source, ra, dec, inver, outver):
    if source == '':
        source=findcal(indata, '')
    clcor             = AIPSTask('CLCOR')
    clcor.indata      = indata
    clcor.source[1]   = source
    clcor.opcode      = 'ANTC'
    clcor.clcorprm[5] = ra
    clcor.clcorprm[6] = dec
    clcor.gainver     = inver
    clcor.gainuse     = outver
    if ra!=0 or dec!=0:
        clcor()

##############################################################################
#
def findcal(indata, calsource):
    if calsource == '':
        n = 0
        for source in indata.sources:
            if source[0]=='G':
                calsource=source
                n=n+1
        if n>1:
            print 'More than one Maser source! Using '+calsource

    return calsource

##############################################################################
#
def findcvelsource(indata, cvelsource):
    if cvelsource == '' or cvelsource == ['']:
        cvelsource = []
        n = 0
        for source in indata.sources:
            if source[0]=='G' and source[1:3]!='RB':
                cvelsource.append(source)
                n=n+1
    return cvelsource

##############################################################################
#
def findtarget(indata, target):
    targets=[]
    for entry in target:
        targets.append(entry)

    if targets == ['']:
        targets = []
        n = 0
        for source in indata.sources:
            if source[0]=='J':
                targets.append(source)
                n=n+1
    return targets

##############################################################################
#
def get_split_sources(indata, target, cvelsource, calsource):

    split_sources=findtarget(indata, target)

    if calsource in cvelsource or calsource in target:
        pass
    else:
        split_sources.append(calsource)

    return split_sources

##############################################################################
#
def mafringe(indata, fr_image, calsource, channel, refant, outdisk, doband, bpver,dpfour):
    split                = AIPSTask('SPLIT')
    split.indata         = indata
    split.source[1]      = calsource
    split.bchan          = channel
    split.echan          = channel
    split.outcl          = 'CH'+str(channel)
    split.outdisk        = outdisk
    split.bif            = 0
    split.eif            = 0
    split.flagver        = 1
    split.docal          = 1
    split.aparm[1]       = 2
    split.aparm[6]       = 1
    split.doband         = doband
    split.bpver          = bpver

    splitdata = AIPSUVData(calsource, 'CH'+str(channel), outdisk, 1)
    if splitdata.exists():
        splitdata.clrstat()
        splitdata.zap()
    split()

    if len(calsource)>=4:
        reducedname = calsource[0:4]
    else:
        reducedname = calsource

    splitdata2 = AIPSUVData(reducedname, 'CH'+str(channel), outdisk, 1)
    if splitdata2.exists():
        splitdata2.clrstat()
        splitdata2.zap()

    splitdata.rename(reducedname, 'CH'+str(channel), 1)

    multi         = AIPSTask('MULTI')
    multi.indata  = splitdata
    multi.outname = reducedname
    multi.outcl   = '1IF'
    multi.outdisk = outdisk

    multidata = AIPSUVData(reducedname, multi.outcl, outdisk, 1)

    if multidata.exists():
        multidata.clrstat()
        multidata.zap()
    multi()

    indxr           = AIPSTask('INDXR')
    indxr.indata    = multidata
    indxr.cparm[1:] = [0,0./60.]
    indxr()

    vbglu         = AIPSTask('VBGLU')
    vbglu.indata  = multidata
    vbglu.in2data = multidata
    vbglu.outname = reducedname
    vbglu.outdisk = outdisk
    vbglu.outcl   = '2IF'

    vbgludata = AIPSUVData(reducedname, '2IF', outdisk, 1)

    if vbgludata.exists():
        vbgludata.clrstat()
        vbgludata.zap()
    vbglu()

    vbglu.indata  = vbgludata
    vbglu.in2data  = vbgludata
    vbglu.outcl   = '4IF'

    vbgludata2 = AIPSUVData(reducedname, '4IF', outdisk, 1)

    if vbgludata2.exists():
        vbgludata2.clrstat()
        vbgludata2.zap()
    vbglu()

    vbgludata.zap()

    vbglu.indata  = vbgludata2
    vbglu.in2data  = vbgludata2
    vbglu.outcl   = '8IF'

    vbgludata3 = AIPSUVData(reducedname, '8IF', outdisk, 1)

    if vbgludata3.exists():
        vbgludata3.clrstat()
        vbgludata3.zap()
    vbglu()

    vbgludata2.zap()

    indxr.indata    = vbgludata3
    indxr()

    fringe               = AIPSTask('FRING')

    if fr_image.exists():
        fringe.in2data = fr_image
        mprint('################################################',logfile)
        mprint('Using input model '+fringe.in2name+'.'+fringe.in2class+'.'+str(int(fringe.in2seq))+' on diks '+str(int(fringe.in2disk)), logfile)
        mprint('################################################',logfile)
    else:
        mprint('################################################',logfile)
        mprint('Using point source as imput model for fringe.',logfile)
        mprint('################################################',logfile)

    fringe.indata        = multidata
    fringe.refant        = refant
    fringe.docal         = 1
    fringe.calsour[1]    = ''
    fringe.solint        = 6
    fringe.aparm[1:]     = [2, 0, 0, 0, 0]
    fringe.dparm[1:]     = [1, -1, 0, 0]
    fringe.dparm[4]      = dpfour
    fringe.snver         = 0
    fringe()

    fringe.indata        = vbgludata3
    fringe()

    return multidata, vbgludata3

##############################################################################
#
def mafringe2(indata, calsour, channel, refant, outdiks, doband, bpver, dpfour):

    fringe               = AIPSTask('FRING')

    if fr_image.exists():
        fringe.in2data = fr_image
        mprint('################################################',logfile)
        mprint('Using input model '+fringe.in2name+'.'+fringe.in2class+'.'+str(int(fringe.in2seq))+' on disk '+str(int(fringe.in2disk)), logfile)
        mprint('################################################',logfile)
    else:
        mprint('################################################',logfile)
        mprint('Using point source as imput model fro fringe.',logfile)
        mprint('################################################',logfile)

    fringe.indata        = indata
    fringe.refant        = refant
    fringe.bchan         = channel
    fringe.echan         = channel
    fringe.docal         = 1
    fringe.calsour[1:]   = [calsour]
    fringe.solint        = 6
    fringe.aparm[1:]     = [2, 0, 1, 0, 0]
    fringe.dparm[1:]     = [1, -1, 0, 0]
    fringe.dparm[4]      = dpfour
    fringe.snver         = 4
    fringe.doband        = doband
    fringe.bpver         = bpver
    fringe()

##############################################################################
#
def sncor(indata):

    bif = input('Enter IF with the maser line: ')
    sncor              = AIPSTask('SNCOR')
    sncor.indata       = indata
    sncor.opcode       = 'CPSN'
    sncor.bif         = 0
    sncor.eif         = 0
    sncor.sncorprm[1:] = [bif]
    sncor.snver        = 4
    sncor()

##############################################################################
#
def check_sn_ver(indata):
        sn=indata.table('AIPS SN',0)
        return len(sn[0])

##############################################################################
#
def runprtmasn(indata,channel):
    prtab           = AIPSTask('PRTAB')
    prtab.indata    = indata
    prtab.inext     = 'SN'
    prtab.invers    = 0
    prtab.docrt     = -1
    prtab.box[1][1] = 1
    prtab.box[1][2] = 3
    prtab.box[1][3] = 4
    prtab.box[1][4] = 9
    if check_sn_ver(indata)>15:
        prtab.box[2][1] = 15
        prtab.box[2][2] = 17
    else:
        prtab.box[2][1] = 13
        prtab.box[2][2] = 15
    prtab.dohms     = -1
    (year, month, day)=get_observation_year_month_day(indata)
    monthlist=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

    if day<10:
        strday='0'+str(day)
    else:
        strday=str(day)

    namma=str(year)+monthlist[int(month)-1]+strday

    prtab.outprint='PWD:'+namma+'_'+prtab.inname+'_CH'+str(channel)+'.RATE'
    prtab()

##############################################################################
#
def runprtmasu(indata,channel):
    prtab           = AIPSTask('PRTAB')
    prtab.indata    = indata
    prtab.inext     = 'SU'
    prtab.invers    = 0
    prtab.docrt     = -1
    prtab.box[1][1] = 1
    prtab.box[1][2] = 2
    prtab.box[1][3] = 11
    prtab.box[1][4] = 12
    prtab.box[2][1] = 13
    prtab.box[2][2] = 14
    prtab.box[2][3] = 15
    prtab.dohms     = -1
    prtab.ndig      = 4
    (year, month, day)=get_observation_year_month_day(indata)
    monthlist=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

    if day<10:
        strday='0'+str(day)
    else:
        strday=str(day)

    namma=str(year)+monthlist[int(month)-1]+strday

    prtab.outprint='PWD:'+namma+'_'+prtab.inname+'_CH'+str(channel)+'.SU'
    prtab()

##############################################################################
#
def runcvel(indata, cvelsource, vel, inter_flag, doband, bpver):

    print 'Running CVEL.'
    if inter_flag==1:
        cvelsource=check_calsource(indata, cvelsource)

    naxis = indata.header['naxis']
    crpix = indata.header['crpix']

    if isinstance(vel, float) or isinstance(vel, int):
        vel = [vel]

    if naxis[3]!=len(vel):
        print 'You have '+str(naxis[3])+' IFs, and '+str(len(vel))+' velocities.'
        sys.exit()
    (linename,restfreq) = get_line_name(indata)

    setjy = AIPSTask('SETJY')
    setjy.indata = indata
    setjy.source[1:] = cvelsource
    setjy.restfreq[1:] = restfreq
    setjy.optype = ''
    setjy.veltyp = 'LSR'
    setjy.veldef = 'RADIO'
    channum = indata.header['naxis'][2]
    setjy.aparm[1:] = [channum/2.+crpix[2],0]

    for i in range(naxis[3]):
        setjy.sysvel     = vel[i]
        setjy.bif = i+1
        setjy.eif = i+1
        setjy()

    cvel = AIPSTask('CVEL')
    cvel.indata = indata
    cvel.source[1:] = cvelsource
    cvel.outna = cvel.inna
    cvel.outcl = cvel.incl
    cvel.gainuse = 7
    cvel.freqid = 1
    cvel.outseq = cvel.inseq+1
    cvel.outdisk = cvel.indisk
    cvel.doband  = doband
    cvel.bpver   = bpver

    cveldata = AIPSUVData(cvel.inna, cvel.incl, int(cvel.indisk), 2)
    if cveldata.exists():
        cveldata.clrstat()
        cveldata.zap()
    cvel()

    indxr = AIPSTask('indxr')
    indxr.indata = cveldata
    indxr()

def runcvel_lba(indata, cvelsource, vel, inter_flag, doband, bpver, channel):
    #uvsort(indata)
    print 'Running CVEL.'
    if inter_flag==1:
        cvelsource=check_calsource(indata, cvelsource)

    naxis = indata.header['naxis']
    crpix = indata.header['crpix']

    if isinstance(vel, float) or isinstance(vel, int):
        vel = [vel]

    if naxis[3]!=len(vel):
        print 'You have '+str(naxis[3])+' IFs, and '+str(len(vel))+' velocities.'
        sys.exit()
    (linename,restfreq) = get_line_name(indata)

    setjy = AIPSTask('SETJY')
    setjy.indata = indata
    setjy.source[1:] = cvelsource
    setjy.restfreq[1:] = restfreq
    #setjy.optype = 'VCAL'
    setjy.veltyp = 'LSR'
    setjy.veldef = 'RADIO'
    channum = indata.header['naxis'][2]
    if velchan==0: setjy.aparm[1:] = [int(channum/2.)+1, 0]
    else:          setjy.aparm[1:] = [velchan, 0]
    #[channum/2.+crpix[2],0]

    for i in range(naxis[3]):
        setjy.sysvel     = vel[i]
        setjy.bif = i+1
        setjy.eif = i+1
        setjy()

    cvel = AIPSTask('CVEL')
    cvel.indata = indata
    cvel.source[1:] = cvelsource
    cvel.outna = cvel.inna
    cvel.outcl = cvel.incl
    cvel.gainuse = 7
    cvel.freqid = 1
    cvel.outseq = cvel.inseq+1
    cvel.outdisk = cvel.indisk
    cvel.doband  = doband
    cvel.bpver   = bpver
    #cvel.aparm[10] = 1

    cveldata = AIPSUVData(cvel.inna, cvel.incl, int(cvel.indisk), 2)
    if cveldata.exists():
        cveldata.clrstat()
        cveldata.zap()
    cvel()

    #uvsort(cveldata)
    indxr = AIPSTask('indxr')
    indxr.indata = cveldata
    indxr()

##############################################################################
#
def uvsort(indata):
    uvsrt = AIPSTask('uvsrt')
    uvsrt.default
    uvsrt.indata = indata
    uvsrt.outdata = indata
    uvsrt.sort = 'TB'
    uvsrt.go()

##############################################################################
#
def get_velocity(indata,cvelsource):
    cvelsource=findcvelsource(indata, cvelsource)
    light = 299792458
    nchan = indata.header['naxis'][2]
    res   = indata.header['cdelt'][2]
    su    = indata.table('AIPS SU',0)
    nif   = indata.header['naxis'][3]
    for entry in su:
        if entry['source'].strip() in cvelsource:
            vel=entry['lsrvel']
            restfreq=entry['restfreq']

    if nif>1:
        restfreq=restfreq[0]
    spacing=res/restfreq*light
    return vel, restfreq, nchan, res, spacing

##############################################################################
#
def get_real_sources(indata):
    nx=indata.table('AIPS NX',0)
    su=indata.table('AIPS SU',0)

    real_sources=[]
    for entry in nx:
        if (su[entry['source_id']-1]['source'].rstrip() in real_sources)==False:
            real_sources.append(su[entry['source_id']-1]['source'].rstrip())

    return real_sources

##############################################################################
#
def check_calsource(indata, calsource):

    if isinstance(calsource, str):
        sour = raw_input('Using source '+calsource+'? (y/n) ')
    else:
        cals_str = ''
        for i in calsource:
            cals_str=cals_str+i+' '
        sour = raw_input('Using source '+cals_str+'? (y/n) ')

    if sour=='n' or sour=='N':
        print 'Searching for sources with data.'
        real_sources=get_real_sources(indata)
        print 'Available sources: '
        for i in range(len(real_sources)):
             print str(i+1)+': '+real_sources[i]

        calsource=raw_input('Which source? ')
        if (calsource in real_sources)==False:
           try:
               k=int(calsource)
               calsource=real_sources[k-1]
           except:
               print 'No such source.'
               sys.exit()
        else:
           print 'No such source.'
           sys.exit()

    return calsource

##############################################################################
#
def get_antname(an_table,n):
    name=''
    for entry in an_table:
        if n==entry['nosta']:
            name=entry['anname']
    return name

##############################################################################
#
def get_phasecal_sources(indata,mp_source,logfile):
    if mp_source == ['']:
        mp_source = []
        for source in indata.sources:
            if source[0]=='F':
                mp_source.append(source)

    newdata = AIPSUVData(indata.name, 'UVCOP', indata.disk, 1)
    if newdata.exists():
        newdata.zap()

    uvcop        = AIPSTask('UVCOP')
    uvcop.indata = indata
    uvcop.source[1:] = mp_source
    uvcop.go()

    return newdata

def check_RDBE(indata,logfile,inter_flag,dtype):

    mprint('################################################',logfile)
    mprint('### Checking geoblock data for RDBE errors #####',logfile)
    mprint('################################################',logfile)
    avspc             = AIPSTask('AVSPC')
    nchan          = indata.header['naxis'][2]

    data1 = AIPSUVData(indata.name, 'AVSPC', indata.disk, 1)

    if data1.exists():
        data1.zap()

    avspc.indata   = indata
    avspc.outdisk  = indata.disk
    avspc.outclass = 'AVSPC'
    avspc.channel  = nchan
    avspc.go()

    data2 = AIPSUVData(indata.name, 'UVAVG', indata.disk, 1)
    if data2.exists():
        data2.zap()

    uvavg             = AIPSTask('UVAVG')
    uvavg.indata      = AIPSUVData(indata.name, 'AVSPC', indata.disk, 1)
    uvavg.yinc        = 300
    uvavg.outclass    = 'UVAVG'
    uvavg.outdisk     = indata.disk
    uvavg.go()

    block_nr = make_check_RDBE(indata,logfile,inter_flag,dtype)

    data1.zap()
    data2.zap()
    return block_nr

def make_check_RDBE(data,logfile,inter_flag,dtype):

    wdata = WAIPSUVData(data.name,'UVAVG',data.disk,data.seq)
    nif     = data.header['naxis'][3]

    sources = get_sources(data)
    su_table = data.table('AIPS SU', 0)
    an_table = data.table('AIPS AN', 0)

    times=[]
    for visibility in wdata:
        there = False
        if visibility.time in times:
            pass
        else:
            if len(times)>0:
                for time in times:
                    if visibility.time-time<6.94e-4:
                        there = True
            else: pass

            if there: pass
            else: times.append(visibility.time)

    blocks=[]
    for entry in times:
        if len(blocks)==0:
            blocks.append(entry)
        else:
            new_block = True
            for time in blocks:
                if abs(entry-time)<4.17e-2:
                    new_block = False
            if new_block:
                blocks.append(entry)

    antennas = []
    for ant in an_table:
        antennas.append(ant['nosta'])

    if len(antennas)<max(antennas):
        for i in range(1,max(antennas)):
            if i in antennas: pass
            else: antennas.append(i)

    antennas=sort(antennas)

    block_nr=0
    block_data=[]
    bad = 0
    flagline=[]
    bt=check_geo(data)

    for block in blocks:
        block_nr+=1

        if dtype=='GEO':
            mprint('#####################################',logfile)
            mprint('###### Checking geoblock Nr.:'+str(block_nr)+' ######', logfile)
            mprint('#####################################',logfile)
        elif dtype=='CONT':
            mprint('######################################',logfile)
            mprint('###### Checking contblock Nr.:'+str(block_nr)+' ######', logfile)
            mprint('######################################',logfile)
        antdata = []
        for ant in antennas:
            if_data=[]
            for IF in range(nif):
               if_data.append([])
            antdata.append(if_data)

        avgantdata = []
        for ant in antennas:
            if_data=[]
            for IF in range(nif):
                if_data.append([])
            avgantdata.append(if_data)

        n = 0
        for visibility in wdata:
            if abs(block-visibility.time)<4.17e-2:
                n+=1
                for IF in range(nif):
                   vis = visibility.visibility[IF][0][0]
                   amp=(sqrt(vis[0]**2+vis[1]**2))
                   antdata[visibility.baseline[0]-1][IF].append(amp)
                   antdata[visibility.baseline[1]-1][IF].append(amp)

        for ant in antennas:
            for IF in range(nif):
                if len(antdata[ant-1][IF])>0:
                    avg=average(antdata[ant-1][IF])
                else:
                    avg=0.0
                avgantdata[ant-1][IF]=avg

        plot=1
        an_table = data.table('AIPS AN', 0)
        for entry in antennas:
            myantdata=array(avgantdata[entry-1])
            ymax=(1.3*max(myantdata)+0.000001)
            avg = average(myantdata)*1000
            rms = std(myantdata)*1000
    #            print get_antname(an_table,entry),average(myantdata)/min(myantdata)
            if min(myantdata)>0 and average(myantdata)/min(myantdata)>2:
                if 'MK' in get_antname(an_table,entry) or 'SC' in get_antname(an_table,entry):
                    if average(myantdata)/min(myantdata)>2.2:
                        bad_IF=[]
                        for IF in range(nif):
                            if average(myantdata)/myantdata[IF]>2.2:
                                bad_IF.append(IF+1)
                        mprint('Probable RDBE Error: '+get_antname(an_table,entry).strip()+', bad IFs:'+str(bad_IF),logfile)
                        bad+=1
                else:
                    bad_IF=[]
                    for IF in range(nif):
                        if average(myantdata)/myantdata[IF]>2:
                            bad_IF.append(IF+1)
                    if len(bad_IF)==1 and 7 in bad_IF:
                        mprint('Probable IF 7 Error: '+get_antname(an_table,entry).strip()+', bad IFs:'+str(bad_IF),logfile)
                    else:
                        mprint('Probable RDBE Error: '+get_antname(an_table,entry).strip()+', bad IFs:'+str(bad_IF),logfile)
                        flagline.append([get_antname(an_table,entry).strip(),bt[block_nr-1],bt[block_nr]])
                    bad+=1
            fig=figure(block_nr,figsize=(8, 12))
            ax=fig.add_subplot(len(antennas),1,plot)
            ax.set_xlim(-0.8,nif-0.5)
            ax.set_ylim(0.0, ymax*1000)
            ax.text(0.03, 0.60, get_antname(an_table,entry), transform=ax.transAxes)
            if round(ymax*900,1)<0.01:
                yticks([0,0.1])
            else:
                yticks([0,round(ymax*900,1)])
            ax.plot(myantdata[:]*1000, 'bx')
            ax.plot([0,nif-1],[avg,avg],'k-')

            if entry==min(antennas):
                title('Aplitude in arbitrary units across the IFs\n Block:'+str(block_nr))
            if entry<max(antennas):
                ax.xaxis.set_major_locator(NullLocator())
            if entry==max(antennas):
                xlabel('No. IF - 1')
            plot+=1

        #        print 'test'
        draw()
        if dtype=='GEO':
            savefig('RDBE_check_geoblock'+str(block_nr)+'.ps')
        elif dtype=='CONT':
            savefig('RDBE_check_contblock'+str(block_nr)+'.ps')

        if inter_flag==1:
            if int(matplotlib.__version__[0])==0:
                if int(matplotlib.__version__[2:4])<99:
                    print 'Close figure to continue.'
                    show()
                    close()
                else:
                    fig.show()
                    raw_input('Press enter to close figure and continue. ')
                    close()
            else:
                print 'Close figure to continue.'
                show()
                close()

    mprint('########################################',logfile)
    mprint('#### '+str(bad)+' probable RDBE errors detected ###',logfile)
    mprint('########################################',logfile)

    date=get_observation_year_month_day(data)
    doy = get_day_of_year(date[0],date[1],date[2])

    if bad>0:
        flagfile=data.name+'.flagfile'
        f=open(flagfile,'w')
        for entry in flagline:
            d1,h1,m1,s1=time_to_hhmmss(entry[1])
            d2,h2,m2,s2=time_to_hhmmss(entry[2])
            f.writelines('ANT_NAME=\''+entry[0]+'\' TIMERANG='+str(d1+doy)+','+str(h1)+','+str(m1)+','+str(s1)+','+str(d2+doy)+','+str(h2)+','+str(m2)+','+str(s2)+' REASON=\'RDBE problem\'/\n')
        f.close()


        if inter_flag==1:
            cont=raw_input('Apply flags? (y/n) ')
            if cont=='y' or cont=='Y':
                runuvflg(data,flagfile,logfile)
        else:
            runuvflg(data,flagfile,logfile)

    return len(blocks)

#    if data.exists():
#        print 'Yes'
#    else:
#        print 'No'

##############################################################################
#
def runpossm(indata, calsource, refant, tv, doband, bpver):

    indata.zap_table('PL', -1)

    calsource=check_calsource(indata, calsource)

    bchan = input('Enter bchan: ')
    echan = input('Enter echan: ')

    if bchan==echan and echan!=0:
        print 'Only one channel selected.'

    elif bchan>echan:
        print 'bchan > echan'

    else:
        possm             = AIPSTask('POSSM')
        if indata.header['naxis'][3]>1:
            possm.bif = input('Enter IF: ')
            possm.eif = possm.bif
        else:
            possm.eif = 1
            possm.bif = 1
        possm.indata      = indata
        possm.source[1:]  = [calsource]
        possm.antenna[1:] = [refant]
        possm.stokes      = 'RR'
        possm.docal       = 1
        possm.gainuse     = 0
        possm.nplots      = 9
        possm.bchan       = bchan
        possm.echan       = echan
        possm.doband      = doband
        possm.bpver       = bpver

        if AIPSTV.AIPSTV().exists():
            possm.dotv        = 1
            gv_flag           = 0
        else:
            possm.dotv        = -1
            gv_flag           = 1
        possm()


        if gv_flag==1:
            if os.path.exists('possm.ps'):
                os.popen('rm possm.ps')
            lwpla         = AIPSTask('LWPLA')
            lwpla.indata  = indata
            lwpla.inver   = 0
            lwpla.outfile = 'PWD:possm.ps'
            lwpla()

            indata.zap_table('PL', -1)
            os.popen('gv possm.ps')
            os.popen('rm possm.ps')

##############################################################################
#
def run_snplt(indata, inter_flag):

    indata.zap_table('PL', -1)
    n_ant         = len(get_ant(indata))
    snplt         = AIPSTask('SNPLT')
    snplt.indata  = indata
    snplt.stokes  = 'RR'
    snplt.inver   = 4
    snplt.inext   = 'SN'
    snplt.optype  = 'PHAS'
    snplt.nplots  = n_ant
    snplt.pixrange[1:] = [-200,200]
    snplt.bif     = 1
    snplt.eif     = 1
    snplt.dotv    = -1
    snplt()

    #snplt.stokes  = 'DIFF'
    #snplt()

    name1=indata.name+'_sn4.ps'
    #name2=indata.name+'_diffsn4.ps'

    if os.path.exists(name1):os.popen('rm '+name1)
    #if os.path.exists(name2):os.popen('rm '+name2)

    lwpla         = AIPSTask('LWPLA')
    lwpla.indata  = indata
    lwpla.inver   = 1
    lwpla.outfile = 'PWD:'+name1
    lwpla()

    #lwpla.inver   = 2
    #lwpla.outfile = 'PWD:'+name2
    #lwpla()

    if inter_flag==1:
        tv=AIPSTV.AIPSTV()
        if tv.exists()==False:
            tv.start()
        if tv.exists():
            tv.clear()

        if AIPSTV.AIPSTV().exists():
            snplt.dotv        = 1
            snplt()
        else:
            os.popen('gv '+name1)

    indata.zap_table('PL', -1)

def run_snplt_diff(indata, inter_flag):

    indata.zap_table('PL', -1)
    n_ant         = len(get_ant(indata))
    snplt         = AIPSTask('SNPLT')
    snplt.indata  = indata
    snplt.stokes  = 'DIFF'
    snplt.inver   = 4
    snplt.inext   = 'SN'
    snplt.optype  = 'PHAS'
    snplt.nplots  = n_ant
    snplt.pixrange[1:] = [-200,200]
    snplt.bif     = 1
    snplt.eif     = 1
    snplt.dotv    = -1
    snplt()

    name2=indata.name+'_diffsn4.ps'

    if os.path.exists(name2):os.popen('rm '+name2)

    lwpla         = AIPSTask('LWPLA')
    lwpla.indata  = indata
    lwpla.inver   = 1
    lwpla.outfile = 'PWD:'+name2
    lwpla()

    if inter_flag==1:
        tv=AIPSTV.AIPSTV()
        if tv.exists()==False:
            tv.start()
        if tv.exists():
            tv.clear()

        if AIPSTV.AIPSTV().exists():
            snplt.dotv        = 1
            snplt()
        else:
            os.popen('gv '+name1)

    indata.zap_table('PL', -1)

def setup_plotfiles(indata):
    (year, month, day)=get_observation_year_month_day(indata)
    ant=get_ant(indata)
    code=''
    mon_dir={1:'jan',2:'feb',3:'mar',4:'apr',5:'may',6:'jun',7:'jul',8:'aug'
             ,9:'sep',10:'oct',11:'nov',12:'dec'}
    path0=('http://www.vlba.nrao.edu/astro/VOBS/astronomy/')
    if code=='':
        code_str         = indata.header.observer.lower()
    else:
        code_str      = code

    naxis    = indata.header['naxis']
    frq=get_center_freq(indata)/1e9
    if frq>1.3 and frq< 1.8:  band='L'
    elif frq>2.1 and frq<2.4: band='S'
    elif frq>4.5 and frq<8:
        if check_sx(indata, logfile): band='S/X'
        else: band='C'
    elif frq>8.0 and frq<10.0:  band='X'
    elif frq>11.5 and frq<16.0: band='U'
    elif frq>21.5 and frq<24.5: band='K'
    elif frq>40.5 and frq<45.5: band='Q'
    else:
        band='unknown'

    dinfo    = band+'-Band: '+str(naxis[1])+' Stokes, '+str(naxis[3])+' IFs with '+str(naxis[2])+' channels'


    path1=path0+mon_dir[month]+str(year)[2:4]+'/'+code_str+'/'
    path2=path1+code_str+'log.vlba'

    sources = get_sources(indata)

    f = open('./index.html','w')
    f.writelines('<html>\n')
    f.writelines('<A NAME="TOP">')
    f.writelines('Datafile: '+indata.name+'<br>\n')
    f.writelines('Observed on '+str(year)+'/'+str(month)+'/'+str(day)+'<br>\n')
    f.writelines('<A HREF = "'+path1+'">NRAO VOBS page</A>, <A HREF="'+path2+'">Logfile</A>')
    f.writelines('<hr>\n')
    f.writelines(dinfo+'<br>\n')
    f.writelines('Antennas: '+str(ant)+'<br>\n')
    f.writelines('Tsys plots: ')
    for entry in ant:
        f.writelines('<A HREF = "'+path1+code_str+'tsm.'+ant[entry]+'.ps.gz">'+ant[entry]+'</A> ')
    f.writelines('<br><br>\n')
    f.writelines('Sources: '+str(sources)+'\n')
    f.writelines('<hr>\n')
    f.writelines('<A HREF = "#SN1">SN1</A> <A HREF = "#SN2">SN2</A> <A HREF = "#SN3">SN3</A> <A HREF = "#SN4">SN4</A> <A HREF = "#RDBE">RDBE Check</A><br>')
    f.writelines('<hr>\n')

    if (os.path.exists('delay-rate.ps')):
        f.writelines('<TABLE BORDER="3" CELLPADDING="10" CELLSPACING="10">\n')
        f.writelines('<TD><A HREF="plotfiles/delay-rate.png"><img SRC="plotfiles/delay-rate.png" width=800></A>\n')
    if (os.path.exists('ATMOS_NAME.FITS')):
        f2=open('./ATMOS_NAME.FITS', 'r')
        data=f2.readlines()
        f2.close()
        f.writelines('<TD><table border=1>\n')
        f.writelines('<tr><th>Antenna<th>UT Time<th>zenith delay<th>clock offset<th>clock rate\n')

        for entry in data:
            if len(entry.split())==1:
                last_ant=''
            else:
                time=entry.split()[1]+':'+entry.split()[2]+':'+entry.split()[3]+':'+entry.split()[4]
                print last_ant
                if entry.split()[0]!=last_ant:
                    f.writelines('<tr><td><td><td><td><td>\n')
                    last_ant=entry.split()[0]
                else:
                    last_ant=entry.split()[0]
                f.writelines('<tr><td align="center">'+entry.split()[0]+'<td align="center">'+time+'<td align="center">'+entry.split()[5]+'<td align="center">'+entry.split()[6]+'<td align="center">'+entry.split()[8]+'\n')
        f.writelines('</table></table>\n')

    if (os.path.exists('delay-rate-ionos.ps')):
        f.writelines('<TABLE BORDER="3" CELLPADDING="10" CELLSPACING="10">\n')
        f.writelines('<TD><A HREF="plotfiles/delay-rate-ionos.png"><img SRC="plotfiles/delay-rate-ionos.png" width=800></A>\n')
    if (os.path.exists('IONOS_NAME.FITS')):
        f2=open('./IONOS_NAME.FITS', 'r')
        data=f2.readlines()
        f2.close()
        f.writelines('<TD><table border=1>\n')
        f.writelines('<tr><th>Antenna<th>UT Time<th>zenith delay<th>clock offset<th>clock rate\n')

        for entry in data:
            if len(entry.split())==1:
                last_ant=''
            else:
                time=entry.split()[1]+':'+entry.split()[2]+':'+entry.split()[3]+':'+entry.split()[4]
                print last_ant
                if entry.split()[0]!=last_ant:
                    f.writelines('<tr><td><td><td><td><td>\n')
                    last_ant=entry.split()[0]
                else:
                    last_ant=entry.split()[0]
                f.writelines('<tr><td align="center">'+entry.split()[0]+'<td align="center">'+time+'<td align="center">'+entry.split()[5]+'<td align="center">'+entry.split()[6]+'<td align="center">'+entry.split()[8]+'\n')
        f.writelines('</table></table>\n')

    f.close()

##############################################################################
#
def run_snplt_2(indata, inver, optype, snplt_name):

    indata.zap_table('PL', -1)
    nif   = indata.header['naxis'][3]
    n_ant         = len(get_ant(indata))
    snplt         = AIPSTask('SNPLT')
    snplt.indata  = indata
    snplt.stokes  = ''
    snplt.inver   = inver
    snplt.inext   = 'SN'
    snplt.optype  = optype
    if nif==1:
        snplt.nplots  = n_ant
    else:
        snplt.nplots  = 8
    snplt.bif     = 1
    snplt.eif     = 0
    snplt.dotv    = -1
    snplt()

    if (os.path.exists('plotfiles')==False):
        os.mkdir('plotfiles')

    f = open('./index.html','a')
    f.writelines('<A NAME="'+snplt_name+'">\n')
    f.writelines(snplt_name+' <A HREF = "#TOP">TOP</A><br>\n')
    for n in range(1,indata.table_highver('AIPS PL')+1):
        name=indata.name.strip()+'_'+snplt_name+'_'+str(n)+'.ps'
        name2=indata.name.strip()+'_'+snplt_name+'_'+str(n)+'.png'

        if os.path.exists(name):
            os.popen('rm '+name)

        lwpla         = AIPSTask('LWPLA')
        lwpla.indata  = indata
        lwpla.plver   = n
        lwpla.inver   = n
        lwpla.dparm[5] = 1
        lwpla.outfile = 'PWD:'+name
        lwpla()
        f.writelines('<A HREF="plotfiles/'+name+'"><img SRC="plotfiles/'+name2+'" width=500></A>\n')
        os.popen(r'convert '+name+' '+name2)
        os.popen(r'mv '+name+' plotfiles/')
        os.popen(r'mv '+name2+' plotfiles/')

    f.writelines('<br><hr>\n')
    f.close()
    indata.zap_table('PL', -1)

def rdbe_plot(data,logfile,dtype):
    block_nr=check_RDBE(data,logfile,0,dtype)
    if dtype=='GEO':
        os.popen(r'pstoimg RDBE_check_geoblock* plotfiles/')
        os.popen(r'cp  RDBE_check_geoblock* plotfiles/')


    f = open('./index.html','a')
    f.writelines('<A NAME="RDBE">\n')
    f.writelines('RDBE check <A HREF = "#TOP">TOP</A><br>\n')
    for i in range(1,block_nr+1):
        rdbplot='RDBE_check_geoblock'+str(i)+'.ps'
        rdbplot2='RDBE_check_geoblock'+str(i)+'.png'
        f.writelines('<A HREF="plotfiles/'+rdbplot+'"><img SRC="plotfiles/'+rdbplot2+'" width=500></A>\n')

    f.writelines('<br><hr>\n')
    f.close()

##############################################################################
#
def run_split(indata, source, outclass, doband, bpver):

    channels = indata.header['naxis'][2]
    if channels==16:
        bad=1                     # remove 1 channel from each side
    elif channels==32:
        bad=2                     # remove 2 channels from each side
    elif channels==64:
        bad=4                     # remove 2 channels from each side
    elif channels==128:
        bad=6                     # remove 6 channels from each side
    elif channels==256:
        bad=8                     # remove 8 channels from each side
    elif channels==512:
        bad=10                    # remove 10 channels from each side
    elif channels==1024:
        bad=12                    # remove 12 channels from each side
    else:
        bad=0

    [bchan,echan]=[1+bad,channels-bad]

    print 'Averaging channels '+str(bchan)+' - '+str(echan)+'.'

    target           = findtarget(indata, source)
    if isinstance(target, str):
        target=[target]

    if target!=[]:
        for source in target:
            split_data=AIPSUVData(source,outclass,indata.disk,1)
            if split_data.exists():
                split_data.clrstat()
                split_data.zap()

        split            = AIPSTask('SPLIT')
        split.indata     = indata
        split.bchan      = bchan
        split.echan      = echan
        split.docalib    = 1
        split.flagver    = 0
        split.source[1:] = target
        split.outclass   = outclass
        split.aparm[1:]  = [3,0]
        split.aparm[6]   = 1
        split.outdisk    = indata.disk
        split.doband     = doband
        split.bpver      = bpver
        split.smooth[1:] = smooth

        split()

def run_fittp_data(source, outcl, disk,logfile):
    fittp         = AIPSTask('FITTP')
    data          = AIPSUVData(source, outcl, disk, 1)
    fittp.indata  = data
    if os.path.exists(fittp.inname+'.'+fittp.inclass+'.fits'):
        os.popen(r'rm '+fittp.inname+'.'+fittp.inclass+'.fits')
    fittp.dataout = 'PWD:'+fittp.inname+'.'+fittp.inclass+'.fits'

    if data.exists():
        mprint('Writing out calibrated and splitted uv-data for '+source,logfile)
        fittp.go()
    else:
        mprint('No calibrated and splitted uv-data for '+source,logfile)

def run_grid(indata, source, cellsize, imsize, n, m, grid_offset,uvwtfn, robust,beam):

    grid=[]
    if n%2!=0:
        for i in range(n):
            if m%2!=0:
                for j in range(m):
                    grid.append([i-n/2,j-m/2])

    for shift in grid:
        shift_pos(indata, source, shift[0]*grid_offset, shift[1]*grid_offset, 8, 9)

        channels = indata.header['naxis'][2]
        if channels==16:
            bad=1                     # remove 1 channel from each side
        elif channels==32:
            bad=2                     # remove 2 channels from each side
        elif channels==64:
            bad=4                     # remove 2 channels from each side
        elif channels==128:
            bad=6                     # remove 6 channels from each side
        elif channels==256:
            bad=8                     # remove 8 channels from each side
        elif channels==512:
            bad=10                    # remove 10 channels from each side
        elif channels==1024:
            bad=12                    # remove 12 channels from each side
        else:
            bad=0

        [bchan,echan]=[1+bad,channels-bad]

        split_data=AIPSUVData(source,'GRIDS',indata.disk,1)
        if split_data.exists():
            split_data.clrstat()
            split_data.zap()

        split            = AIPSTask('SPLIT')
        split.indata     = indata
        split.bchan      = bchan
        split.echan      = echan
        split.docalib    = 1
        split.flagver    = 0
        split.source[1:] = [source]
        split.outclass   = 'GRIDS'
        split.aparm[1:]  = [3,0]
        split.aparm[6]   = 1
        split.outdisk    = indata.disk
        split()

        print 'Imaging '+str(shift[0])+' '+str(shift[1])
        rungridimagr(split_data, source, 10, cellsize, imsize, -1,uvwtfn,robust,beam)
        split_data.zap()

        restore_su(indata, logfile)
        check_sncl(indata, 4, 8,logfile)

def run_masplit(indata, source, outclass, doband, bpver, smooth, channel):
    source           = findcal(indata, source)
    if isinstance(source, str):
        source=[source]

    if source!=[]:
        for target in source:
            split_data=AIPSUVData(target,outclass,indata.disk,1)
            if split_data.exists():
                split_data.clrstat()
                split_data.zap()

        split            = AIPSTask('SPLIT')
        split.indata     = indata
        split.bchan      = 0
        split.echan      = 0
        split.flagver    = 0
        split.docalib    = 1
        split.source[1:] = source
        split.outclass   = outclass
        split.aparm[1:]  = [0,0]
        split.aparm[6]   = 1
        split.outdisk    = indata.disk
        split.doband     = doband
        split.bpver      = bpver
        split.smooth[1:] = smooth
        split()

##############################################################################
# phase_selfcal
def phase_selfcal(indata, source, solint, outdisk, niter, cellsize, imsize,
                  logfile, imna, cal_ants, im_ants, refant, startmod, beam):

    mprint('#############################################################',logfile)
    mprint('### Phase Self-cal on: '+source+(35-len(source))*' '+'###',logfile)

    if imna=='':
        outname = source
    else:
        outname = source[:11-len(imna)]+imna
    nimage = 1
    while AIPSImage(outname,'ICL001',line_data.disk, nimage+1).exists():
        nimage+=1

    model  = AIPSImage(outname,'ICL001',line_data.disk, nimage)


    if (indata.table_highver('AIPS SN'))>0:
        while indata.table_highver('AIPS SN')>0:
            indata.zap_table('AIPS SN', 0)

    if (indata.table_highver('AIPS NX'))>0:
        while indata.table_highver('AIPS NX')>0:
            indata.zap_table('AIPS NX', 0)

    calib              = AIPSTask('CALIB')
    calib.indata       = indata
    calib.calsour[1:]  = [source]
    calib.aparm[1:]    = [0,0,1,0,1,0]
    calib.soltype      = 'L1R'
    calib.solmode      = 'P'
    calib.solint       = solint
    calib.outdisk      = outdisk
    calib.outname      = outname
    calib.refant       = refant
    calib.antennas[1:] = cal_ants

    if model.exists():
        calib.in2data     = model
        mprint('### Using Model: '+str(model)+'###', logfile)
        mprint('### Solint = '+str(int(solint*60))+' seconds ###', logfile)
        mprint('#############################################################'
               ,logfile)
        xpixel = model.header['naxis'][0]
        ypixel = model.header['naxis'][1]

        ccedt           = AIPSTask('CCEDT')
        ccedt.indata    = model
        ccedt.nboxes    = 1
        ccedt.clbox[1][1:] = [1,1,xpixel,ypixel]
        ccedt.cutoff    = -2*model.header['datamin']
        ccedt()

        calib()
    elif startmod.exists():
        calib.in2data = startmod
        mprint('### Using input model '+calib.in2name+'.'+calib.in2class+'.'+str(int(calib.in2seq))+' on diks '+str(int(calib.in2disk))+' ###', logfile)
        calib()
    else:
        mprint('### Using Point source model                              ###'
               , logfile)
        mprint('### Solint = '+str(int(solint*60))+' seconds ###', logfile)
        mprint('#############################################################'
               ,logfile)
        calib.smodel[1:] = [1,0]
        calib()

    ncal = 1
    while AIPSUVData(outname,'CALIB',line_data.disk,ncal+1).exists():
        ncal+=1
    cal_data=AIPSUVData(outname,'CALIB',line_data.disk,ncal)
    runimagr(cal_data, source, niter, cellsize, imsize, -1, imna, im_ants,uvwtfn, robust,beam,
        baselines=ant_bls,timer=imgr_timer)

    return cal_data


##############################################################################
# amp_selfcal
def amp_selfcal(indata, source, solint, outdisk, niter, cellsize,
                imsize, logfile, imna, antennas, refant, dofit,beam):

    mprint('#############################################################',logfile)
    ants=''
    for j in range(len(dofit)):
        if dofit[j]!=0:
            if j!=len(dofit):
                ants=ants+str(get_ant(indata)[entry])+', '
            else:
                ants=ants+str(get_ant(indata)[entry])
    mprint('### Ampl. Self-cal on: '+source+(35-len(source))*' '+'###',logfile)
    mprint('### Using dofit = '+str(dofit)+' ('+str(ants)+') ###',logfile)

    if imna=='':
        outname = source
    else:
        outname = source[:11-len(imna)]+imna
    nimage = 1
    while AIPSImage(outname,'ICL001',line_data.disk, nimage+1).exists():
        nimage+=1

    model  = AIPSImage(outname,'ICL001',line_data.disk, nimage)

    if (indata.table_highver('AIPS SN'))>0:
        while indata.table_highver('AIPS SN')>0:
            indata.zap_table('AIPS SN', 0)

    if (indata.table_highver('AIPS NX'))>0:
        while indata.table_highver('AIPS NX')>0:
            indata.zap_table('AIPS NX', 0)

    calib              = AIPSTask('CALIB')
    calib.indata       = indata
    calib.calsour[1:]  = [source]
    calib.aparm[1:]    = [0,0,1,0,1,0]
    calib.soltype      = 'L1R'
    calib.solmode      = 'A&P'
    calib.dofit[1:]    = dofit
    calib.solint       = solint
    calib.outdisk      = outdisk
    calib.outname      = outname
    calib.refant       = refant
    calib.cparm[2]     = 1

    if model.exists():
        calib.in2data     = model
        mprint('### Using Model: '+str(model)+'###', logfile)
        mprint('### Solint = '+str(int(solint*60))+' seconds ###', logfile)
        mprint('#############################################################'
               ,logfile)
        xpixel = model.header['naxis'][0]
        ypixel = model.header['naxis'][1]

        ccedt           = AIPSTask('CCEDT')
        ccedt.indata    = model
        ccedt.nboxes    = 1
        ccedt.clbox[1][1:] = [1,1,xpixel,ypixel]
        ccedt.cutoff    = -2.0*model.header['datamin']
        ccedt()

        calib()
    else:
        mprint('### Using Point source model                              ###'
               , logfile)
        mprint('### Solint = '+str(int(solint*60))+' seconds ###', logfile)
        mprint('#############################################################'
               ,logfile)
        calib.smodel[1:] = [1,0]
        calib()

    ncal = 1
    while AIPSUVData(outname,'CALIB',line_data.disk,ncal+1).exists():
        ncal+=1
    cal_data=AIPSUVData(outname,'CALIB',line_data.disk,ncal)
    runimagr(cal_data, source, niter, cellsize, imsize, -1, imna, antennas,uvwtfn,robust,beam,
        baselines=ant_bls,timer=imgr_timer)

    return cal_data

#############################################################################
# apply_selfcal
def do_apply_selfcal(source,calsource,split_outcl,outdisk,refant):

    single_data  = AIPSUVData(source,split_outcl,outdisk,0)

    multi          = AIPSTask('MULTI')
    multi.indata   = single_data
    multi.outclass = split_outcl
    multi.outseq   = 2
    multi()

    multi_data = AIPSUVData(source,split_outcl,outdisk,2)

    cal_data  = AIPSUVData(calsource,split_outcl,outdisk,1)

    tacop         = AIPSTask('TACOP')
    tacop.indata  = cal_data
    tacop.outdata = multi_data
    tacop.inext   = 'SN'
    tacop.inver   = 0
    if source==calsource: pass
    else: tacop()

    runclcal(multi_data,1,1,2,'',1,refant)

    ncal = 0
    while AIPSUVData(calsource,'CALIB',outdisk,ncal+2).exists():
        ncal+=1
        cal_data=AIPSUVData(calsource,'CALIB',outdisk,ncal)
        tacop.indata = cal_data
        tacop()
        runclcal(multi_data,ncal+1,ncal+1,ncal+2,'',1,refant)

################################################################################
# Extract line name from a uvdata
#
def get_line_name(indata):
    freq = get_center_freq(indata)/1e9
    if freq>12 and freq<13:
        restfreq = [1.2178E+10,597000]
        linename='CH3OH_12GHz'
        print 'Assuming 12.2 GHz methanol maser.'
    elif freq>22 and freq<23:
        restfreq = [2.2235E+10,80000]
        linename='H2O'
        print 'Assuming 22.2 GHz water maser.'
    elif freq>6 and freq<7:
        restfreq = [6.668e+09, 519200]
        linename='CH3OH_6.7GHz'
        print 'Assuming 6.7 GHz methanol maser.'
    elif freq>43.1 and freq<43.2:
        restfreq = [43.122e+09,80000]
        linename='SiO_43.1GHz'
        print 'Assuming 43.122 GHz SiO maser.'
    else:
        print 'Unknown maser line.'
        exit()

    return (linename, restfreq)

##############################################################################
# run possm to make an ampscalar spectrum using only the inner-5 antennae
#
def runrpossm(indata, cvelsource, tv, interflag,antennas):

    indata.zap_table('PL', -1)

    if indata.table_highver('AIPS SN') ==0:
        docal=0
    else:
        docal=1

    if interflag==1:
#        nchvel = input('Enter 1 for vel or 2 for ch: ')
#        chvel_str=['vel','ch']
#        chvel=chvel_str[nchvel-1]
        chvel='vel'

        if indata.header['naxis'][3]>1:
            bif = input('Enter IF: ')
            eif = bif
        else:
            eif=1
            bif=1

        possm             = AIPSTask('POSSM')
        possm.indata      = indata
        #possm.source[1:]  = cvelsource[0]
        possm.source[1:]  = cvelsource
        #possm.antenna[1:] = [2,4,5,8,9]
        possm.antenna[1:] = antennas
        possm.stokes      = 'I'
        possm.docal       = docal
        possm.gainuse     = 0
        possm.nplots      = 0
        possm.bchan       = 0
        possm.echan       = 0
        possm.bif         = bif
        possm.eif         = eif
        possm.aparm[8]    = 0
        possm.aparm[1]    = -1
        possm.codetype    = 'AMP'
        if AIPSTV.AIPSTV().exists():
            possm.dotv        = 1
            gv_flag           = 0
        else:
            possm.dotv        = -1
            gv_flag           = 1
        possm()

        bchan = input('Enter bchan: ')
        echan = input('Enter echan: ')

    else:
        chvel ='vel'
        bchan = 0
        echan = 0
        bif   = 1
        eif   = 1


    if bchan==echan and echan!=0:
        print 'Only one channel selected.'

    elif bchan>echan:
        print 'bchan > echan'

    else:
        possm             = AIPSTask('POSSM')
        possm.indata      = indata
        #possm.source[1:]  = cvelsource[0]
        possm.source[1:]  = cvelsource
        #possm.antenna[1:] = [2,4,5,8,9]
        possm.antenna[1:] = [3,4,7,8]
        possm.stokes      = 'I'
        possm.docal       = docal
        possm.gainuse     = 0
        possm.nplots      = 0
        possm.bchan       = bchan
        possm.echan       = echan
        possm.bif         = bif
        possm.eif         = eif
        possm.aparm[8]    = 0
        possm.aparm[1]    = -1
        possm.codetype    = 'AMP'

        (linename,restfreq) = get_line_name(indata)

#--- x-axis labeled with channel number
        if chvel == 'ch':
            if bchan==0 and echan==0:
                spectfile='PWD:'+cvelsource[0]+'.'+linename+'.spectrum.txt'
                filename=cvelsource[0]+'.'+linename+'.POSSM.ps'
            else:
                spectfile =''
                filename=cvelsource[0]+'.'+linename+'.POSSM_CH'+str(bchan)+'-'+str(echan)+'.ps'

            if gv_flag == 0:
                possm.dotv = -1
                possm.outtext = spectfile
                possm()

            if os.path.exists(filename):
                os.popen('rm '+filename)
            lwpla         = AIPSTask('LWPLA')
            lwpla.indata  = indata
            lwpla.inver   = 0
            lwpla.outfile = 'PWD:'+filename
            lwpla()

#--- x-axis labeled with velocity
        if chvel == 'vel':
            possm.aparm[7]     = 2
            possm.aparm[10]    = 1
            possm.dotv = -1
            possm()

            filename=cvelsource[0]+'.'+linename+'.POSSM.ps'

            if os.path.exists(filename):
                os.popen('rm '+filename)
            lwpla         = AIPSTask('LWPLA')
            lwpla.indata  = indata
            lwpla.inver   = 0
            lwpla.outfile = 'PWD:'+filename
            lwpla()
        os.popen(r'ps2pdf '+filename)

        indata.zap_table('PL', -1)
        #os.popen('gv '+filename)

###############################################################################
# Get peak and rms from an image
#
def runimean(imgdata, blc=[0,0,0], trc=[0,0,0]):
    '''Must set imgdata'''
    if os.path.exists('imean.txt'):
        os.popen('rm imean.txt')
    imean = AIPSTask('imean')
    imean.indata = imgdata
    imean.blc[1:] = blc
    imean.trc[1:] = trc
    if int(aipsver[6])>1:
        imean.doprint = 1
    imean.outtext = 'PWD:imean.txt'
    imean()
    #datamax = imgdata.header.datamax
    datamax = get_image_peak()
    return (datamax, imean.pixstd)

##############################################################################
#
def get_image_peak(idir='./'):
    infile = idir+"imean.txt"
    if os.path.exists(infile):
        myfile = open(infile)

        for line in myfile.readlines():
            strtmp = 'Maximum'
            if cmp(strtmp,line[0:7]) == 0:
                peak = line[9:20]
        myfile.close()

        peak=float(peak)
        return peak
    else:
        print('###--- '+infile)
        print('###--- Input file dose not exsit')

################################################################################
# run SAD for maser source
#
def run_ma_sad(inimg,indata,cut,dyna):
    if inimg.exists():
        inimg.clrstat()
        num_ch = inimg.header.naxis[3-1]
        srcname = inimg.name
        obs=inimg.header['observer'];
        date_obs=inimg.header['date_obs']
        (linename,restfreq) = get_line_name(indata)
        for line in inimg.history:
            if line[0:11] == 'IMAGR BCHAN':
                bchan = int(line[15:22])
            if line[0:11] == 'IMAGR ECHAN':
                echan = int(line[15:22])
                break

        bch_str = str(bchan+1000)
        bch_str = bch_str[1:4]
        ech_str = str(echan+1000)
        ech_str = ech_str[1:4]

        sad_file = srcname+'_'+linename+'_'+obs+'_'+\
                date_obs+'_CH'+bch_str+'-'+ech_str+'_sad.txt'
        if os.path.exists(sad_file):
            os.popen('rm '+sad_file)
        if os.path.exists('sad.txt'):
            os.popen('rm sad.txt')

        sad           =  AIPSTask('SAD')
        sad.indata    = inimg
        sad.doresid   = -1
        sad.ngauss    = 3
        #--- for weak source, we may need to increase gain
        sad.gain      = 0.3
        sad.sort      = 'S'
        sad.docrt     = -4
        sad.fitout    = 'PWD:'+'sad.txt'
        mpeak=[]
        mrms=[]
        for i in range(1,num_ch+1):
            (peak, rms) = runimean(inimg,[0,0,i],[0,0,i])
            mpeak.append(peak)
            mrms.append(rms)

            c1 = peak*0.5
            c3 = max(rms*cut,peak*dyna)
            c2 = (c1 + c3)*0.5

            d1 = c3
            d2 = c3
            d3 = 0
#--- pixel (depended on the compactness)
#--- cellsize = 0.1 mas, then 40 means 4 mas
            d4 = 40
            d5 = 40
            cparm=[c1,c2,c3]
            cparm.sort()
            cparm.reverse()
            sad.cparm[1:] = cparm
            sad.icut      = c3
            sad.dparm[1:] = [d1,d2,d3,d4,d5]
            sad.blc[3]    = i
            sad.trc[3]    = i
            if peak>rms*cut:
                sad()

        cmd='cp sad.txt '+sad_file
        os.popen(cmd)
        mprint ('Channel:   Peak     rms      SNR',logfile)
        mprint ('            (Jy)   (Jy)',logfile)
        mprint ('---------------------------------',logfile)
        for i in range(len(mpeak)):
            mprint(('%3d     %7.4f %7.4f %8.4f' % (i+bchan-1,
                  mpeak[i], mrms[i], mpeak[i]/mrms[i])),logfile)
    else:
        print 'Image does not exist, please check it.'
        exit()

def orfit_to_plot_sorted(inimg,indata,bchan,vel):
    f=open('sad.txt')
    content=f.readlines()
    f.close()

    nif     = indata.header['naxis'][3]
    altrval = inimg.header['altrval']
    altrpix = inimg.header['altrpix']

    srcname = inimg.name
    obs=inimg.header['observer'];
    date_obs=inimg.header['date_obs']
    (linename,restfreq) = get_line_name(indata)

    cvelsource=inimg.name
    (vel2, restfreq, nchan, res, spacing) = get_velocity(line_data,cvelsource)

    BW = indata.header['cdelt'][2]
    channum = indata.header['naxis'][2]
    n_center  = (channum/2.) + 1
    dV        = -spacing/1000
    vel2 = inimg.header['altrval']

    sinfo={}
    ichan0    = bchan
    for entry in content:
        if entry[0]=='!':
            pass
        else:
            name=entry.split()[0]
            sname=entry[0:10].strip()
#            schan=(int(entry[12:16].strip()))+ichan0-1
            schan=(int(entry[12:16].strip()))-altrpix+1
            sffla=entry[31]
            if sffla==' ':
                factor=1
            elif sffla=='m':
                factor=0.001
            else:
                factor=1000
            sflux=(float(entry[16:24].strip())*factor)
            sfler=(float(entry[24:31].strip())*factor)
            sraof=(float(entry[32:44].strip())/1000)
            sraer=(float(entry[44:52].strip())/1000)
            sdeof=(float(entry[52:64].strip())/1000)
            sdeer=(float(entry[64:72].strip())/1000)
            srera=(entry[73:84].strip())
            srede=(entry[84:96].strip())
#            svelo=vel[0]+(schan-n_center)*dV
            svelo=(vel2-(schan-1)*spacing)/1000
            sinfo[sflux]=[schan,sfler,svelo,sraof,sraer,sdeof,sdeer]


    filename = srcname+'.'+linename+'.spotmap.txt'
    f2=open(filename,'w')
    f2.writelines('!  Results for '+obs+' '+srcname+' observations\n')
    f2.writelines('!  of '+str(round(restfreq/1e9,1))+' GHz '+linename
                  +' masers conducted on '+date_obs+'\n')
    f2.writelines('!  Reference position: '+srera+' '+srede+'\n')
    f2.writelines('!  Channel  Vlsr     S_peak       dX           +-'+
                  '         dY           +-\n')
    f2.writelines('!    #     (km/s)    (Jy/bm)      (")         (")'+
                  '         (")         (")\n')

    order=sorted(sinfo)
    order.reverse()
    for i in order:
#            print '   %4i    %6.2f    %6.2f  %10.6f  %10.6f  %10.6f  %10.6f' % (sinfo[i][0],sinfo[i][2],i,sinfo[i][3],sinfo[i][4],sinfo[i][5],sinfo[i][6])
            f2.writelines('   %4i    %6.2f    %6.2f  %10.6f  %10.6f  %10.6f  %10.6f\n' % (sinfo[i][0],sinfo[i][2],i,sinfo[i][3],sinfo[i][4],sinfo[i][5],sinfo[i][6]))

    f2.close()

def orfit_to_plot_sorted2(inimg,indata):
    f=open('sad.txt')
    content=f.readlines()
    f.close()

    srcname = inimg.name
    obs=inimg.header['observer'];
    date_obs=inimg.header['date_obs']
    (linename,restfreq) = get_line_name(indata)

    cvelsource=inimg.name
    (vel, restfreq, nchan, res, spacing) = get_velocity(line_data,cvelsource)

    altrval=inimg.header['altrval']
    altrpix = inimg.header['altrpix']

    n_center  = (nchan/2.) + 1
    dV        = -spacing

    sinfo={}
    for entry in content:
        if entry[0]=='!':
            pass
        else:
            name=entry.split()[0]
            sname=entry[0:10].strip()
            schan=(int(entry[12:16].strip()))-altrpix+1
            sffla=entry[31]
            if sffla==' ':
                factor=1
            elif sffla=='m':
                factor=0.001
            else:
                factor=1000
            sflux=(float(entry[16:24].strip())*factor)
            sfler=(float(entry[24:31].strip())*factor)
            sraof=(float(entry[32:44].strip())/1000)
            sraer=(float(entry[44:52].strip())/1000)
            sdeof=(float(entry[52:64].strip())/1000)
            sdeer=(float(entry[64:72].strip())/1000)
            srera=(entry[73:84].strip())
            srede=(entry[84:96].strip())
            svelo=(altrval+(schan-1)*dV)/1000
            sinfo[sflux]=[schan,sfler,svelo,sraof,sraer,sdeof,sdeer]



    filename = srcname+'.'+linename+'.spotmap.txt'
    f2=open(filename,'w')
    f2.writelines('!  Results for '+obs+' '+srcname+' observations\n')
    f2.writelines('!  of '+str(round(restfreq/1e9,1))+' GHz '+linename
                  +' masers conducted on '+date_obs+'\n')
    f2.writelines('!  Reference position: '+srera+' '+srede+'\n')
    f2.writelines('!  Channel  Vlsr     S_peak       dX           +-'+
                  '         dY           +-\n')
    f2.writelines('!    #     (km/s)    (Jy/bm)      (")         (")'+
                  '         (")         (")\n')

    order=sorted(sinfo)
    order.reverse()
    for i in order:
#            print '   %4i    %6.2f    %6.2f  %10.6f  %10.6f  %10.6f  %10.6f' % (sinfo[i][0],sinfo[i][2],i,sinfo[i][3],sinfo[i][4],sinfo[i][5],sinfo[i][6])
            f2.writelines('   %4i    %6.2f    %6.2f  %10.6f  %10.6f  %10.6f  %10.6f\n' % (sinfo[i][0],sinfo[i][2],i,sinfo[i][3],sinfo[i][4],sinfo[i][5],sinfo[i][6]))

    f2.close()

def plot_spot_map(inimg,indata):
    srcname = inimg.name
    obs=inimg.header['observer']
    xsize=abs(inimg.header['naxis'][0]*inimg.header['cdelt'][0])*3600
    ysize=inimg.header['naxis'][1]*inimg.header['cdelt'][1]*3600
    date_obs=inimg.header['date_obs']
    (linename,restfreq) = get_line_name(indata)
    filename = srcname+'.'+linename+'.spotmap.txt'
    f=open(filename)
    content=f.readlines()
    f.close()
    x=[]
    y=[]
    z=[]
    f=[]
    f_max=0
    n_poi=0
    for entry in content:
        if entry[0]=='!':
            pass
        else:
            n_poi+=1
            x_off = float(entry.split()[3])
            y_off = float(entry.split()[5])
            if abs(x_off)<0.48*xsize and abs(y_off)<0.48*ysize:
                x.append(x_off)
                y.append(y_off)
                z.append(float(entry.split()[1]))
                flux=float(entry.split()[2])
                f.append(flux**0.5)
                if float(entry.split()[2])>f_max:
                    f_max=float(entry.split()[2])
    min_x=min(x)
    max_x=max(x)
    dx=max_x-min_x
    min_y=min(y)
    max_y=max(y)
    size=max(max_y-min_y,max_x-min_x,0.2)/2
    x_cent = (max_x+min_x)/2
    y_cent = (max_y+min_y)/2
    dy=max_y-min_y
    f2=array(f)
    scale = 300/max(f2)
    f1=figure()

    ra,dec= deg_to_radec([inimg.header['crval'][0],inimg.header['crval'][1]])
    title(srcname+' '+linename+' on '+date_obs+'\n (0,0)='+ra+','+
          dec+' (J2000)')
    ax=axes()
    s=ax.scatter(x,y,c=z,s=f2*scale,marker='o',edgecolors='',cmap='plasma')
    ax.scatter(x_cent+1.05*size,y_cent-1.15*size,c='0.8',s=max(f2)*scale)
    cb = plt.colorbar(s)
    cb.set_label('v$_{LSR}$ [km s$^{-1}$]')
    xlabel('East Offset [arcseconds]')
    ylabel('North Offset [arcseconds]')
    ax.set_xlim(x_cent+1.2*size,x_cent-1.2*size)
    ax.set_ylim(y_cent-1.3*size,y_cent+1.1*size)
    xy=(x_cent+0.93*size,y_cent-1.17*size)
    ax.annotate('= '+str(round(f_max,1))+
            ' Jy bm$^{-1}$', xy, xytext=None, xycoords='data',
         textcoords='data', arrowprops=None,
         bbox=dict(boxstyle="round", fc="0.8"),)
    draw()
    imname = srcname+'.'+linename+'.spotmap.ps'
    #show()
    savefig(imname)
    os.popen('ps2pdf '+imname)

def deg_to_radec(pos):
    if pos[0]<0:
        pos[0]=360+pos[0]
    ra_deg=pos[0]/15.
    dec_deg=pos[1]
    hour = int(ra_deg)
    min  = int(60*(ra_deg-hour))
    sec  = (60*(60*(ra_deg-hour)-min))
    if abs(60-round(sec,5))<1e-5:
        sec=0
        min+=1
    if min==60:
        min=0
        hour+=1
    if hour<0: hour=hour+24
    if hour>=24: hour=hour-24
    hp=''
    mp=''
    sp=''
    if hour<10: hp='0'
    if min<10: mp='0'
    if sec<10: sp='0'
    ra=hp+str(hour)+':'+mp+str(min)+':'+sp+str(round(sec,5))
    deg  = abs(int(dec_deg))
    amin  = int(60*(abs(dec_deg)-deg))
    asec  = (60*(60*(abs(dec_deg)-deg)-amin))
    if 60-abs(round(asec,4))<1e-4:
        asec=0
        amin+=1
    if amin==60:
        amin = 0
        deg+=1
    if dec_deg<0:
        sign='-'
    else:
        sign='+'
    dp=''
    amp=''
    asp=''
    if deg<10: dp='0'
    if amin<10: amp='0'
    if asec<10: asp='0'
    dec = sign+dp+str(abs(deg))+':'+amp+str(abs(amin))+':'+asp+str(round(abs(asec),4))
    return ra, dec

##############################################################################
# check_cal

def check_cal(indata, source, outdisk, logfile, imna):
    if imna=='':
        outname = source
    else:
        outname = source[:11-len(imna)]+imna

    ncal = 1
    while AIPSUVData(outname,'CALIB',line_data.disk,ncal+1).exists():
        ncal+=1
    cal_data=AIPSUVData(outname,'CALIB',line_data.disk,ncal)
    if cal_data.exists():
        return cal_data
    else:
        return indata

##############################################################################
# Get download options from mail

def read_mail(mail_path, ou, op, of, inter_flag):
    f=open(mail_path)
    content=f.readlines()
    f.close()

    users=[]
    passs=[]
    dates=[]
    pubfiles=[]
    date_pubfiles_nr={}
    pubfile_nr=[]
    k = 0
    kk = 0
    j = 0

    for n in range(len(content)):
        if 'From: VLA/VLBA' in content[n]:
            j=j+1
            date=content[n-1].rstrip()
            dates.append(date)
        pfa = 'Public File available :'
        if pfa in content[n] and (pfa in content[n-1])==False:
            while pfa in content[n+kk]:
                kk=kk+1
                k=k+1
                pubfile_nr.append(k)
            date_pubfiles_nr[j]=pubfile_nr
        pubfile_nr=[]
        kk=0


    for entry in content:
        if 'Proprietary File Dir :' in entry:
            prop_dir = entry.split(' ')[4].rstrip()
        elif 'Your access username, password :' in entry:
            user  = entry.split(' ')[5]
            passw = entry.split(' ')[7].rstrip()
            users.append(user)
            passs.append(passw)
        elif 'Public File available :' in entry:
            user = 'nopass'
            passw = ''
            users.append(user)
            passs.append(passw)
            pubfiles.append(entry.split('/')[4].rstrip())


    if len(dates)>1:
        print 'Found VLA/VLBA archive emails from:'
        for i in range(len(dates)):
            print str(i+1)+': '+dates[i]

        if inter_flag==1:
            date_nr=input('Which one? (1-'+str(len(dates))+') ')-1
            if (date_nr in range(len(dates)))==False:
                print 'No such date.'
                sys.exit()
        else:
            date_nr=len(dates)-1

        print 'Using: '+date
        print ''
        user=users[date_nr]
        passw=passs[date_nr]
        print user, passw


    if user=='nopass':
        file_names = []
        for entry in date_pubfiles_nr[date_nr+1]:
            file_names.append(pubfiles[entry-1])
    else:
        os.popen(r'wget --no-check-certificate --http-user '+user+
                 ' --http-passwd '+passw+' -t45 -O files.txt'+
                 ' https://archive.nrao.edu/secured/'+user+'/', 'w')
        print ' https://archive.nrao.edu/secured/'+user+'/'

        f=open('files.txt')
        content=f.read()
        f.close()

        t=content.split('"')
        file_names=[]
        for entry in t:
            if 'fits' in entry:
                file_names.append(entry)

        os.popen(r'rm files.txt')

    print_download_options(user,passw, file_names)

    if inter_flag==1:
        cont=raw_input('Use these filenames? (y/n) ')
        if cont=='n' or cont=='N':
            print 'Using other filenames:'
            print_download_options(ou,op,of)
            return ou, op, '', of, len(of)

        else:
            return user, passw, prop_dir, file_names, len(file_names)
    else:
        return user, passw, prop_dir, file_names, len(file_names)

def get_download_names(ou, op, of):
    user=ou
    passw=op
    if user=='nopass':
        file_names = []
        file_sizes = []
        for entry in date_pubfiles_nr[date_nr+1]:
            file_names.append(pubfiles[entry-1])
            file_sizes.append('')
    else:
        os.popen(r'wget --no-check-certificate --http-user '+user+
                 ' --http-passwd '+passw+' -t45 -O files.txt'+
                 ' https://archive.nrao.edu/secured/'+user+'/', 'w')

        f=open('files.txt')
        content=f.read()
        f.close()

        t=content.split('"')
        t=content.split('\n')
        file_names=[]
        file_sizes=[]
        for entry in t:
            if 'fits' in entry:
                t2=entry.split('"')
                file_names.append(t2[1])
                file_sizes.append(entry[-8:])
        os.popen(r'rm files.txt')

    if inter_flag==1 and of[0]==0:
        print_download_options(user,passw, file_names,file_sizes)
        cont=raw_input('Use these filenames? (y/n) ')
        if cont=='n' or cont=='N':
            order=[]
            for i in range(0,len(file_names)):
                entry = file_names[i]
                print entry,file_sizes[i]
                cont=raw_input('File number? ')
                if int(cont)>len(file_names)-1:
                    print 'Only Number 0 to '+str(len(file_names)-1)+' allowed.'
                    cont=raw_input('File number? ')
                order.append(cont)


            new_file_names=range(0,len(file_names))
            for i in range(0,len(file_names)):
                new_file_names[int(order[i])]=file_names[i]
            file_names = new_file_names
            print_download_options(user,passw,file_names,file_sizes)
            return file_names, len(file_names)

        else:
            return file_names, len(file_names)

    elif inter_flag==1 and of[0]!=0:
        print_download_options(user,passw, of, file_sizes)
        cont=raw_input('Use these filenames? (y/n) ')
        if cont=='n' or cont=='N':
            order=[]
            for entry in file_names:
                print entry
                cont=raw_input('File number? ')
                if int(cont)>len(file_names)-1:
                    print 'Only Number 0 to '+str(len(file_names)-1)+' allowed.'
                    cont=raw_input('File number? ')
                order.append(cont)

            new_file_names=range(0,len(file_names))
            for i in range(0,len(file_names)):
                new_file_names[int(order[i])]=file_names[i]
            file_names = new_file_names
            print_download_options(user,passw,file_names, file_sizes)
            return file_names, len(file_names)

        else:
            return of, len(of)

    else:
        if of[0]!=0:
            print_download_options(ou,op,of,file_sizes)
            return of, len(of)
        else:
            print_download_options(user,passw, file_names,file_sizes)
            return file_names, len(file_names)

def print_download_options(user,passw,file_names,file_sizes):
    print '    file        = range(n)'
    print '    arch_user   = \''+user+'\''
    print '    arch_pass   = \''+passw+'\''
    print ''

    n=0
    for i in range(0,len(file_names)):
        print '    file['+str(n)+'] = \''+file_names[i]+'\' ('+file_sizes[i]+')'
        n=n+1
    print ''

def check_data(data, n, geo, cont, line, logfile):
    count = 0
    for i in range(len(data)):
        if data[i].exists():
            count=count+1
            data_info(data[i],i, geo, cont, line, logfile)
    if count==n: mprint('Found '+str(count)+' data files on disk.',logfile)
    else: mprint('Expected '+str(n)+' files, but found '+str(count)+' data files on disk.',logfile)

def data_info(indata, i, geo, cont, line, logfile):
    if indata.exists():
        frq=get_center_freq(indata)/1e9
        if frq>1.3 and frq< 1.8:  band='L'
        elif frq>2.1 and frq<2.4: band='S'
        elif frq>4.5 and frq<8:
            if check_sx(indata, logfile): band='S/X'
            else: band='C'
        elif frq>8.0 and frq<10.0:  band='X'
        elif frq>11.5 and frq<16.0: band='U'
        elif frq>21.5 and frq<24.5: band='K'
        elif frq>40.5 and frq<45.5: band='Q'
        else:
            band='unknown'
        if i==geo: add = ' (geoblock data)'
        elif i==cont and i!=line: add = ' (continuum data)'
        elif i==line and i!=cont: add = ' (line data)'
        elif i==line and i==cont: add = ' (line or continuum data)'
        else: add = ' (data not used)'
        naxis    = indata.header['naxis']
        mprint('File '+str(i)+': '+indata.name+' '+band+' band, '+str(naxis[1])
              +' Stokes, '+str(naxis[2])+' channels '+str(naxis[3])
              +' IFs'+add, logfile)
        return band, naxis[1], naxis[2], naxis[3]

##############################################################################
# Added by LJH

# deletes beam
def _zapbeam(source,inseq=1):
    beam=AIPSImage(source,'IBM001',1,inseq)
    if beam.exists():
        beam.zap()



# loads files
def get_file(path):
    #opens and external file and makes it into a list
    fopen = path
    f=open(fopen, 'r+')
    g=list(f)
    g=map(lambda s: s.strip(), g)
    return g

#splits strings into lists/numpy arrays
def splitt(old_list):
    #splits the list entries into sublists
    new_list=[]
    for i in old_list:
        new_list+=[i.split()]
    return np.array(new_list)

# extra AIPS tasks
def runcalib(indata,sources=[''],gainuse=0,docal=-1,snver=0,solmode='',soltype='',aparm7=0,chan=0,refant=0):
    calib             = AIPSTask('CALIB')
    calib.indata      = indata
    calib.calsour[1:] = sources
    calib.gainu       = gainuse
    calib.snver       = snver
    calib.solmode     = solmode
    calib.soltype     = soltype
    calib.aparm[3]    = 1
    calib.aparm[7]    = aparm7
    calib.refant      = refant
    calib()

def runsplat(indata,outdata,sources=[''],stokes='HALF',docal=-1):
    splat = AIPSTask('splat')
    splat.indata      = indata
    splat.outdata     = outdata
    splat.sources[1:] = sources
    splat.stokes      = stokes
    splat.docalib     = docal
    splat.gainu       = 0
    splat.aparm[1:]   = [2,0]
    splat()

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

# make multiview control file from indata=calibrator in2data=linedata
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

# python version of sed            
def replace(file, searchExp, replaceExp):
    for line in fileinput.input(file, inplace=1):
        line = line.replace(searchExp, replaceExp)
        sys.stdout.write(line)

# END defs
##############################################################################


#############################################################################
####  Do not change beyond this line unless you know what you are doing! ####
#############################################################################

AIPS.log    = open(logfile, 'a')

mprint('#############################################',logfile)
mprint('### Using definition file from '+version_date+' ###',logfile)
mprint('### Using AIPS Version '+aipsver+(19-len(aipsver))*' '+'###', logfile)
mprint('#############################################',logfile)

try:    debug = debug
except: debug = 0

try:
    if inter_flag==0:
        print 'Running in non-interactive mode.'
    else:
        print 'Running in interactive mode.'
except:
    inter_flag=1
    print 'Running in interactive mode.'

try:
    if split_outcl=='':
        split_outcl = 'SPLIT'
    else:
        if len(split_outcl)>6:
            split_outcl=split_outcl[0:6]
            mprint('################################################',logfile)
            mprint('split_outcl longer than 6 characters. Truncating',logfile)
            mprint('it to: '+split_outcl, logfile)
            mprint('################################################',logfile)
except:
    split_outcl = 'SPLIT'

##############################################################################
#####                 Set default parameters, if not set                 #####
##############################################################################

if 'delzn_flag' in locals() and globals(): pass
else: delzn_flag = 0
if 'restore_fg_flag' in locals() and globals(): pass
else: restore_fg_flag = 0
if 'restore_su_flag' in locals() and globals(): pass
else: restore_su_flag = 0
if 'split_flag' in locals() and globals(): pass
else: split_flag = 0
if 'ma_imagr_flag' in locals() and globals(): pass
else: ma_imagr_flag=0
if 'co_imagr_flag' in locals() and globals(): pass
else: co_imagr_flag   = 0
if 'cube_imagr_flag' in locals() and globals(): pass
else: cube_imagr_flag = 0
if 'fr_n' in locals() and globals(): pass
else: fr_n = ''
if 'fr_c' in locals() and globals(): pass
else: fr_c = ''
if 'fr_d' in locals() and globals(): pass
else: fr_d = defdisk
if 'fr_s' in locals() and globals(): pass
else: fr_s = 1
if 'nmaps' in locals() and globals(): pass
else: nmaps = 1
if 'flux' in locals() and globals(): pass
else: flux = {'' : [0,0,0,0]}
if 'niter' in locals() and globals(): pass
else: niter = 100
if 'grid_flag' in locals() and globals(): pass
else: grid_flag = 0
if 'gridsource' in locals() and globals(): pass
else: gridsource = ''
if 'n_grid' in locals() and globals(): pass
else: n_grid = 0
if 'm_grid' in locals() and globals(): pass
else: m_grid = 0
if 'grid_offset' in locals() and globals(): pass
else: grid_offset = 0
if 'dual_geo' in locals() and globals(): pass
else: dual_geo = 0
if 'arch_user' in locals() and globals(): pass
else: arch_user = ''
if 'arch_pass' in locals() and globals(): pass
else: arch_pass = ''
if 'file' in locals() and globals(): pass
else: file = []
if 'kntr_flag' in locals() and globals(): pass
else: kntr_flag = 0
if 'fittp_flag' in locals() and globals(): pass
else: fittp_flag = 0
if 'get_key_flag' in locals() and globals(): pass
else: get_key_flag = 0
if 'code' in locals() and globals(): pass
else: code = ''
if 'max_ant' in locals() and globals(): pass
else: max_ant = 12
if 'phase_cal_flag' in locals() and globals(): pass
else: phase_cal_flag = 0
if 'amp_cal_flag' in locals() and globals(): pass
else: amp_cal_flag = 0
if 'imna' in locals() and globals(): pass
else: imna = ''
if 'phase_target_flag' in locals() and globals(): pass
else: phase_target_flag=''
if 'amp_target_flag' in locals() and globals(): pass
else: amp_target_flag=''
if 'antennas' in locals() and globals(): pass
else: antennas = [0]
if 'refeed_flag' in locals() and globals(): pass
else: refeed_flag = 0
if 'plot_tables' in locals() and globals(): pass
else: plot_tables = -1
if 'dofit' in locals() and globals(): pass
else: dofit = [0]
if 'apply_selfcal' in locals() and globals(): pass
else: apply_selfcal = 0
if 'tysmo_flag' in locals() and globals(): pass
else: tysmo_flag = 0
if 'solint' in locals() and globals(): pass
else: solint = 0
if 'smodel' in locals() and globals(): pass
else: smodel = [1,0]
if 'uvwtfn' in locals() and globals(): pass
else: uvwtfn = ''
if 'robust' in locals() and globals(): pass
else: robust = 0
if 'bandcal' in locals() and globals(): pass
else: bandcal = ['']
if 'do_band_flag' in locals() and globals(): pass
else: do_band_flag = 0
if 'dpfour' in locals() and globals(): pass
else: dpfour = 0
if 'min_elv' in locals() and globals(): pass
else: min_elv = 0
if 'rpossm_flag' in locals() and globals(): pass
else: rpossm_flag = 0
if 'ma_sad_flag' in locals() and globals(): pass
else: ma_sad_flag = 0
if 'plot_map' in locals() and globals(): pass
else: plot_map = 0
if 'min_snr' in locals() and globals(): pass
else: min_snr = 7
if 'smooth' in locals() and globals(): pass
else: smooth = [0]
if 'beam' in locals() and globals(): pass
else: beam = [0,0,0]
if 'RDBE_check' in locals() and globals():pass
else: RDBE_check = 0
if 'TECU_model' in locals() and globals():pass
else: TECU_model = 'jplg'
if 'imultiv_flag' in locals() and globals(): pass
else: imultiv_flag = 0
if 'imv_imagr_flag' in locals() and globals(): pass
else: imv_imagr_flag = 0
if 'mvwin' in locals() and globals():pass
else: mvwin = 30.0
if 'ant_bls' in locals() and globals(): pass
else: ant_bls = [0]
if 'imv_prep_flag' in locals() and globals(): pass
else: imv_prep_flag = 0
if 'imv_app_flag' in locals() and globals(): pass
else: imv_app_flag = 0
if 'imgr_timer' in locals() and globals(): pass
else: imgr_timer = [0,0,0,0,0,0,0,0]

##############################################################################
# Start main script

mprint('######################',logfile)
mprint(get_time(),logfile)
mprint('######################',logfile)

refant_flag=0
if refant==0: refant_flag=1


if download_flag==1:
    try:
        if os.path.exists(mail_path):
            (arch_user, arch_pass, prop_dir, file, m)=read_mail(mail_path, arch_user, arch_pass, file, inter_flag)
    except:
        if arch_user!='nopass':
            (file,m)=get_download_names(arch_user, arch_pass, file)

            mprint('########################################',logfile)
            mprint('Found '+str(m)+' files to download from archive.', logfile)
            if m<n:
                mprint('Not what expected! Using only '+str(m)
                      +' files.',logfile)
                n=m
            elif m>n:
                mprint('Not what expected! Making '+str(m-n)
                    +' dummy files',logfile)
                for i in range(m-n):
                    filename.append('Dummy'+str(i))
                    outname.append('Dummy'+str(i))
                    outclass.append('UVDATA')
                    nfiles.append(0)
                    ncount.append(1)
                    doconcat.append(-1)
                    flagfile.append('')
                    antabfile.append('')
                    outdisk.append(defdisk)
                n=m
            mprint('########################################',logfile)
        elif arch_user=='nopass' and file[0]==0:
            mprint('###########################################',logfile)
            mprint('Need to specify filenames for public data.',logfile)
            mprint('###########################################',logfile)
            sys.exit()

    for i in range(n):
        mprint('##################################',logfile)
        mprint('Download file '+str(i)+' from archive.', logfile)
        download(arch_user, arch_pass, file[i], file_path, filename[i])
        mprint(get_time(),logfile)
        mprint('##################################',logfile)

if load_flag==1:
    for i in range(n):
        loadindx(file_path,filename[i],outname[i],outclass[i],outdisk[i],
                 nfiles[i],ncount[i],doconcat[i],antname,logfile)

data = range(n)

mprint('##################################',logfile)
for i in range(n):
    data[i]=AIPSUVData(outname[i], outclass[i], int(outdisk[i]), int(1))
    if data[i].exists():
        data[i].clrstat()
        if listr_flag==1:
            runlistr(data[i])
        if get_key_flag==1:
                get_key_file(data[i], code)
mprint('##################################',logfile)

###########################
# Data Preparation

if geo_data_nr < 0: do_geo_block=0
else: do_geo_block=1

pr_data_nr  = []

if cont!=-1:
    pr_data_nr.append(cont)
    if line != cont:
        pr_data_nr.append(line)

check_data(data, n, geo_data_nr, cont, line, logfile)

# Download TEC maps and EOPs

if pr_prep_flag==1 or geo_prep_flag==1:
    if do_geo_block==1:
        (year, month, day)=get_observation_year_month_day(data[geo_data_nr])
        num_days=get_num_days(data[geo_data_nr])
    else:
        (year, month, day)=get_observation_year_month_day(data[pr_data_nr[0]])
        num_days=get_num_days(data[pr_data_nr[0]])

    doy=get_day_of_year(year, month, day)
    
    get_TEC(year,doy,TECU_model)
    if not os.path.exists(eop_path):
        os.mkdir(eop_path)
    get_eop(eop_path)

    if num_days==2: get_TEC(year,doy+1,TECU_model)

mprint('######################',logfile)
mprint(get_time(),logfile)
mprint('######################',logfile)

###########################b########################################
# Geodetic block analysis

if RDBE_check==1:
    blocks = check_RDBE(data[geo_data_nr],logfile,inter_flag,'GEO')

if geo_prep_flag>0:
    geo_data= data[geo_data_nr]
    runuvflg(geo_data,flagfile[geo_data_nr],logfile)
    check_sncl(geo_data, 0, 1,logfile)
    if geo_data.header['telescop']=='EVN':
        if geo_prep_flag==1:
            runTECOR(geo_data,year,doy,num_days,3,TECU_model)
        else:
            runtacop(geo_data, geo_data, 'CL', 1, 3, 0)
    else:
        if geo_prep_flag==1:
            runTECOR(geo_data,year,doy,num_days,2,TECU_model)
        else:
            runtacop(geo_data, geo_data, 'CL', 1, 2, 0)
        runeops(geo_data, eop_path)

    geo_data= data[geo_data_nr]
    sx_geo = False
    if check_sx(geo_data,logfile):
        mprint('####################################',logfile)
        mprint('Geoblock observation in SX mode',logfile)
        mprint('####################################',logfile)
        (geo_data_l, geo_data_h)=split_sx(geo_data, 1)
        if dual_geo>0:
            sx_geo = True
            cw_geo = False
        geo_data = geo_data_h
    elif check_clh(geo_data,logfile):
        mprint('#########################################',logfile)
        mprint('Geoblock observation in wide C-band mode',logfile)
        mprint('########################################',logfile)
        (geo_data_l, geo_data_h)=split_clh(geo_data, 1)
        if dual_geo>0:
            sx_geo = True
            cw_geo = True
        geo_data = geo_data_h

if do_geo_block==1:
    geo_data= data[geo_data_nr]
    sx_geo = False
    if check_sx(geo_data,logfile):
        mprint('####################################',logfile)
        mprint('Geoblock observation in SX mode',logfile)
        mprint('####################################',logfile)
        (geo_data_l, geo_data_h)=split_sx(geo_data, 0)
        if dual_geo>0:
            sx_geo = True
            cw_geo = False
        if geo_data_h.exists(): geo_data = geo_data_h
    elif check_clh(geo_data,logfile):
        mprint('#########################################',logfile)
        mprint('Geoblock observation in wide C-band mode',logfile)
        mprint('########################################',logfile)
        (geo_data_l, geo_data_h)=split_clh(geo_data, 0)
        if dual_geo>0:
            sx_geo = True
            cw_geo = True
        if geo_data_h.exists(): geo_data = geo_data_h

    if delzn_flag==0:
        if (os.path.exists('geoblock')==False):
            os.mkdir('geoblock')

if geo_fringe_flag==1 and do_geo_block==1 and sx_geo==False:
    mprint('##########################################',logfile)
    mprint('Processing geoblock file: '+geo_data.name+'.'+geo_data.klass,logfile)
    mprint('##########################################',logfile)

    check_sncl(geo_data, 0, 3, logfile)
    runquack(geo_data, [0], 2./60.)
    if refant_flag==1: refant=select_refant(geo_data)
    testfringe(geo_data, refant, 0, logfile)
    if refant_flag==1: refant=select_refant2(geo_data,logfile)
    geo_man_pcal(geo_data, refant, logfile)
    runclcal(geo_data, 1, 3, 4, '', 1, refant)
    fringegeo(geo_data, refant)
    get_geosource_stats(geo_data)

if geo_fringe_flag==1 and do_geo_block==1 and sx_geo==True:
    mprint('###################################################',logfile)
    mprint('Processing geoblock files: '+geo_data_l.name+'.'+geo_data_l.klass+
           ' and '+geo_data_h.name+'.'+geo_data_h.klass,logfile)
    mprint('###################################################',logfile)
    check_sncl(geo_data_l, 0, 3, logfile)
    check_sncl(geo_data_h, 0, 3, logfile)
    runquack(geo_data_l, [0], 2./60.)
    runquack(geo_data_h, [0], 2./60.)
    if refant_flag==1: refant=select_refant(geo_data_h)
    testfringe(geo_data_h, refant, 0, logfile)
    if refant_flag==1: refant=select_refant2(geo_data_h,logfile)
    geo_man_pcal_sx(geo_data_l, geo_data_h, refant, logfile)
    runclcal(geo_data_l, 1, 3, 4, '', 1, refant)
    runclcal(geo_data_h, 1, 3, 4, '', 1, refant)
    fringegeo(geo_data_l, refant)
    fringegeo(geo_data_h, refant)
    get_geosource_stats(geo_data_h)

    mprint('######################',logfile)
    mprint(get_time(),logfile)
    mprint('######################',logfile)

if doprt_flag==1 and do_geo_block==1 and delzn_flag==0:
    check_sn_ver(geo_data)
    # Creates directories
    if (os.path.exists('geoblock')==False):
        os.mkdir('geoblock')
    else:
        os.popen(r'rm -rf geoblock/*')


    make_station_file(geo_data)

    if sx_geo==False:
        make_calibrator_file(geo_data)
        make_control_file(geo_data,refant,max_ant)
    else:
        make_calibrator_file_sx(geo_data_h)
        make_control_file_sx(geo_data_h,data[line],refant,max_ant)

    if sx_geo:
        runprtsn_sx(geo_data_l,geo_data_h,data[line])
        runprtsu(geo_data_h)
        os.popen(r'mv 20* geoblock/')
    else:
        runprtsn(geo_data)
        runprtsu(geo_data)
        os.popen(r'mv 20* geoblock/')

if dofit_flag==1 and do_geo_block==1:
    if refant==0: refant=select_refant(geo_data)
    if delzn_flag==1:
        run_delzn(geo_data)
    else:
        if sx_geo==False:
            os.chdir('geoblock')
            if os.path.exists('tropos'):
                os.popen('rm -rf tropos')
            os.mkdir('tropos')
            os.popen('rm -rf fit_geodetic_bl*')
            os.popen(fit_path+'fit_geoblocks_tropos > fit_geoblocks.prt','w')
            os.popen('cat fit_geoblocks.prt','w')
            os.popen('mv fit_geodetic_bl* tropos/')
            os.chdir('../')
            os.popen(r' cp geoblock/ATMOS.FITS .')
            make_name_atmos(geo_data)
        else:
            os.chdir('geoblock')
            if os.path.exists('tropos'):
                os.popen('rm -rf tropos')
            if os.path.exists('ionos'):
                os.popen('rm -rf ionos')
            os.mkdir('tropos')
            os.mkdir('ionos')
            os.popen('rm -rf fit_geodetic_bl*')
            os.popen(fit_path+'diff_dual_freq_delays','w')
            os.popen('cp calibrator_file_tropos.inp calibrator_file.inp')
            os.popen('cp control_file_tropos.inp control_file.inp')
            os.popen(fit_path+'fit_geoblocks_tropos > fit_geoblocks_tropos.prt','w')
            os.popen('cat fit_geoblocks_tropos.prt','w')
            os.popen('mv fit_geodetic_bl* tropos/')
            os.chdir('../')
            os.popen(r' cp geoblock/ATMOS.FITS .')

            os.chdir('geoblock')
            os.popen('cp calibrator_file_ionos.inp calibrator_file.inp')
            os.popen('cp control_file_ionos.inp control_file.inp')
            os.popen(fit_path+'fit_geoblocks_ionos > fit_geoblocks_ionos.prt','w')
            os.popen('cat fit_geoblocks_ionos.prt','w')
            os.popen('mv fit_geodetic_bl* ionos/')
            os.chdir('../')
            os.popen(r' cp geoblock/IONOS.FITS .')
            make_name_atmos(geo_data)
            make_name_ionos(geo_data)

if doplot_flag>0 and do_geo_block==1 and delzn_flag==0:
    if refant==0: refant=select_refant(geo_data)
    plot_baseline(geo_data,refant, inter_flag,'A',doplot_flag,logfile)
    plotatmos(inter_flag, logfile)
    if sx_geo:
        plot_baseline(geo_data,refant, inter_flag,'I',doplot_flag,logfile)
        plotionos(inter_flag, logfile)

mprint('######################',logfile)
mprint(get_time(),logfile)
mprint('######################',logfile)

###################################################################
# Phase referencing analysis

line_data  = data[line]
cont_data  = data[cont]

if RDBE_check==2:
    newdata=get_phasecal_sources(cont_data,mp_source,logfile)
    blocks = check_RDBE(newdata,logfile,inter_flag,'CONT')

if not antname == 'LBA':
    if line_data.sources==cont_data.sources:
       pass
    else:
       mprint('#########################################################',logfile)
       mprint('### Sources in line and continuum files not identical ###',logfile)
       mprint('#########################################################',logfile)
       sys.exit()
else:
    ''

if do_geo_block==1 and pr_prep_flag==1 and delzn_flag==0:
    for i in pr_data_nr:
        if (get_ant(geo_data)==get_ant(data[i])):
            mprint('Antennas identical in geoblock data and file '+str(i),logfile)
        else:
            mprint('## ATTENTION ## ATTENTION ## ATTENTION ##',logfile)
            mprint('Antennas different for geoblock data and file '+str(i),logfile)
            mprint('## ATTENTION ## ATTENTION ## ATTENTION ##',logfile)
    mprint('####################################',logfile)
    checkatmos(inter_flag, logfile)
    mprint('####################################',logfile)

n=0
for i in pr_data_nr:
    n=n+1
    pr_data = data[i]
    mprint('#######################################',logfile)
    mprint('Processing phase-ref file: '+outname[i],logfile)
    mprint('#######################################',logfile)

    if tasav_flag==1:
        runtasav(pr_data, i, logfile)

    if restore_su_flag==1 or pr_prep_flag==1:
        restore_su(pr_data, logfile)

    if restore_fg_flag==1:
        restore_fg(pr_data, logfile)

    if pr_prep_flag>0:
        runuvflg(pr_data,flagfile[i],logfile)
        check_sncl(pr_data, 0, 1,logfile)
        if pr_data.header['telescop']=='EVN':
            if pr_prep_flag==1:
                runTECOR(pr_data,year,doy,num_days,3,TECU_model)
            else:
                runtacop(pr_data, pr_data, 'CL', 1, 3, 0)
        else:
            if pr_prep_flag==1:
                runTECOR(pr_data,year,doy,num_days,2,TECU_model)
            else:
                runtacop(pr_data, pr_data, 'CL', 1, 2, 0)
            runeops(pr_data, eop_path)

        if refant_flag==1:
            refant=select_refant(pr_data)
        if do_geo_block==1:
            if delzn_flag==1:
                atmos_file='DELZN.FITS'
            else:
                geo_data=data[geo_data_nr]
                make_name_atmos(geo_data)
                atmos_file='ATMOS_NAME.FITS'
            runatmos(pr_data, atmos_file)
            if sx_geo==True and dual_geo==1:
                ionos_file='IONOS_NAME.FITS'
                make_name_ionos(geo_data)
                runionos(pr_data, ionos_file)
            runpang(pr_data)
            for source in pos_shift:
                [ra, dec] = [pos_shift[source][0],pos_shift[source][1]]
                if not source in pr_data.sources:
                    continue
                if ra!=0 or dec!=0:
                    if source == '':
                        source=findcal(pr_data, '')
                    mprint('###############################'+
                           '###########################',logfile)
                    mprint('Shift '+source+' by '+str(ra)+
                    ' arcsec in RA and '+str(dec)+
                    ' arcsec in DEC', logfile)
                    mprint('###############################'+
                           '###########################',logfile)
                    shift_pos(pr_data, source, ra, dec, 4, 4)
        else:
            mprint('####################################',logfile)
            mprint('Using no ATMOS.FITS file',logfile)
            mprint('####################################',logfile)
            runpang2(pr_data)
            for source in pos_shift:
                [ra, dec] = [pos_shift[source][0],pos_shift[source][1]]
                if ra!=0 or dec!=0:
                    if source == '':
                        source=findcal(pr_data, '')
                    mprint('####################################'+
                           '################################',logfile)
                    mprint('Shift '+source+' by '+str(ra)+
                    ' arcsec in RA and '+str(dec)+
                    ' arcsec in DEC', logfile)
                    mprint('####################################'+
                           '################################',logfile)
                    shift_pos(pr_data, source, ra, dec, 4, 4)


    if apcal_flag==1:
        vla_pres=False
        check_sncl(pr_data, 0, 4,logfile)
        if flux.items()[0][0]!='':
            ants = get_ant(pr_data)
            for ant in ants:
                if ants[ant]=='Y ':
                    vla=ant
                    vla_pres=True
                    runquack(pr_data, [ant], 10./60.)
            sources  = get_sources(pr_data)
            for entry in sources:
                run_setjy(pr_data, entry, [1,0,0,0])
        for entry in flux:
            if entry!='':
                run_setjy(pr_data, entry, flux[entry])
        if refant_flag==1:
            refant=select_refant(pr_data)
        if antabfile[i]!='':
            runantab(pr_data,antabfile[i])
        if tysmo_flag==1:
            runtysmo(pr_data,90,10)
        if pr_data.header['telescop']=='EVN':
            runapcal(pr_data, 0, 1, 2, 0)
            runclcal(pr_data, 2, 4, 6, '',  1, refant)
        else:
            if antname=='LBA':
                #if not pr_data.table_highver('TY')>=1:
                # no antenna antables pre-applied, run defaults
                runaccor(pr_data)
                runclcal(pr_data, 1, 4, 5, '',  1, refant)
                runapcal_lba(pr_data, 1, 5, 6)
                runtacop(pr_data, pr_data, 'SN', 1, 2, 0)
                mprint('#####################################',logfile)
                mprint('Now run amplitude_cal external script',logfile)
                mprint('#####################################',logfile)
                #else:
                #    #found TY tables, applying
                #    LBA antab is unreliable, don't bother - LJH
                #    runaccor(pr_data)
                #    runclcal(pr_data, 1, 4, 5, '',  1, refant)
                #    runapcal(pr_data, 0, 1, 2, 0)
                #    runclcal(pr_data, 2, 5, 6, '',  1, refant)
            else:
                runaccor(pr_data)
                runclcal(pr_data, 1, 4, 5, '',  1, refant)
                runapcal(pr_data, 0, 1, 2, 0)
                runclcal(pr_data, 2, 5, 6, '',  1, refant)

        mprint('######################',logfile)
        mprint(get_time(),logfile)
        mprint('######################',logfile)

    if pr_fringe_flag==1:
        check_sncl(pr_data, 2, 6,logfile)
        if refant_flag==1:
            refant=select_refant2(pr_data, logfile)
        (so,ti)=man_pcal(pr_data, refant, mp_source, mp_timera,debug,logfile, dpfour)
        if n==1:
            so_ti=[so,ti]
        if n==2:
            if so_ti[0]==so and so_ti[1]==ti:
                mprint('#############################################',logfile)
                mprint('### Both manual phasecal scans identical. ###',logfile)
                mprint('#############################################',logfile)
            else:
                mprint('#############################################',logfile)
                mprint('### Manual phasecal scans different.      ###',logfile)
                mprint('### Select one manually or close enough   ###',logfile)
                mprint('### '+data[1].name+' '+str(so_ti[0])+' '+str(so_ti[1])+' ###',logfile)
                mprint('### '+data[2].name+' '+str(so)+' '+str(ti)+' ###',logfile)
                mprint('#############################################',logfile)
                #sys.exit()
        runclcal(pr_data, 3, 6, 7, '', 1, refant)

    if do_band_flag==1:
        check_sncl(pr_data, 3, 7,logfile)
        do_band(pr_data, bandcal, logfile)

line_data  = data[line]
cont_data  = data[cont]
line_data2 = AIPSUVData(line_data.name,line_data.klass,line_data.disk,2)
cont_data2 = AIPSUVData(cont_data.name,cont_data.klass,cont_data.disk,2)

if bandcal==['']:
    doband = -1
    bpver  = -1
else:
    doband = 1
    bpver  = 1

if cvel_flag==1:
    check_sncl(line_data, 3, 7,logfile)
    check_sncl(cont_data, 3, 7,logfile)
    calsource = findcal(line_data, calsource)
    if isinstance(cvelsource, str):
        cvelsource=[cvelsource]
    if cvelsource[0]=='':
        cvelsource[0]=calsource
    if not antname == 'LBA':
        runcvel(line_data,cvelsource,vel,inter_flag, doband, bpver)
    else:
        runcvel_lba(line_data,cvelsource,vel,inter_flag, doband, bpver, channel)

if possm_flag==1:
    if refant_flag==1:
        refant=select_refant(line_data)
    if line_data2.exists():
        line_data2.clrstat()
    calsource = findcal(line_data, cvelsource[0])

    tv=AIPSTV.AIPSTV()
    if tv.exists()==False: tv.start()

    rep='Y'
    while (rep=='Y' or rep=='y' or rep=='yes'):
        if line_data2.exists():
            mprint('Using shifted data (CVEL).',logfile)
            runpossm(line_data2, calsource, refant, tv, doband, bpver)
        else:
            mprint('Using unshifted data (CVEL).',logfile)
            runpossm(line_data, calsource, refant, tv, doband, bpver)
        rep=raw_input('Repeat POSSM with different channels? (y/n) ')

    channel=input('Enter channel for fringe fit: ')

###################################################################

fr_image=AIPSImage(fr_n, fr_c, fr_d, fr_s)

if ma_fringe_flag==1 and line != cont:
    if line_data2.exists():
        line_data2.clrstat()
        check_sncl(line_data2, 3, 7,logfile)

    check_sncl(line_data, 3, 7,logfile)
    check_sncl(cont_data, 3, 7,logfile)

    calsource = findcal(line_data, calsource)

    mprint('######################',logfile)
    mprint('Selected channel: '+str(channel),logfile)
    mprint('######################',logfile)

    if line_data2.exists():
        mprint('Using shifted data (CVEL).',logfile)
        line_used = line_data2
    else:
        mprint('Using unshifted data (CVEL).',logfile)
        line_used = line_data

    (outdata1,outdata2)=mafringe(line_used, fr_image, calsource, channel, refant, line_data2.disk, doband, bpver, dpfour)
    runtacop(outdata1, line_used, 'SN', 1, 4, 1)
    runtacop(outdata2, cont_data, 'SN', 1, 4, 1)
    runclcal(line_used, 4, 7, 8, '', 1, refant)
    if snflg_flag==1:
        runsnflg(cont_data, 4, calsource)
    if min_elv>0:
        run_elvflag(cont_data,min_elv,logfile)
    runclcal(cont_data, 4, 7, 8, '', 1, refant)
    run_snplt(line_used, inter_flag)
    run_snplt_diff(line_used, inter_flag)

    mprint('######################',logfile)
    mprint(get_time(),logfile)
    mprint('######################',logfile)

if co_fringe_flag==1 and line!=cont:

    check_sncl(cont_data, 3, 7,logfile)
    fringecal(cont_data,fr_image,nmaps,refant,calsource,solint,smodel,doband,bpver,dpfour)
    runclcal(cont_data, 4, 7, 8, '', 1, refant)
    run_snplt(cont_data, inter_flag)

    if line_data2.exists():
        line_used=line_data2
    else:
        line_used=line_data

    line_used.clrstat()
    check_sncl(line_used, 3, 7,logfile)
    runtacop(cont_data, line_used, 'SN', 4, 4, 1)
    if snflg_flag==1:
        runsnflg(line_used, 4, calsource)
    if min_elv>0:
        run_elvflag(line_used,min_elv,logfile)
    runclcal(line_used, 4, 7, 8, '', 1, refant)

    mprint('######################',logfile)
    mprint(get_time(),logfile)
    mprint('######################',logfile)

###################################################################

if ma_fringe_flag==1 and line == cont:

    line_data.clrstat()
    check_sncl(line_data, 3, 7,logfile)

    mprint('######################',logfile)
    mprint('Selected channel: '+str(channel),logfile)
    mprint('######################',logfile)

    calsource = findcal(line_data, calsource)

    if line_data2.exists():
        line_data2.clrstat()
        line_used = line_data2
        mprint('Using shifted data (CVEL).',logfile)
    else:
        mprint('Using unshifted data (CVEL).',logfile)
        line_used = line_data

    check_sncl(line_used, 3, 7,logfile)
    mafringe2(line_used, fr_image,calsource, channel, refant, line_used.disk, doband, bpver,dpfour)
    sncor(line_used)
    if snflg_flag==1:
        runsnflg(line_used, 4, calsource)
    if min_elv>0:
        run_elvflag(line_used,min_elv,logfile)
    runclcal(line_used, 4, 7, 8, '', 1, refant)
    run_snplt(line_used, inter_flag)

if co_fringe_flag==1 and line==cont:

    if cont_data2.exists():
        mprint('Using shifted data (CVEL).',logfile)
        cont_used=cont_data2
    else:
        mprint('Using unshifted data (CVEL).',logfile)
        cont_used=cont_data

    check_sncl(cont_used, 3, 7,logfile)
    fringecal(cont_used,fr_image,nmaps,refant,calsource,solint,smodel,doband,bpver,dpfour)
    if snflg_flag==1:
        runsnflg(cont_used, 4, calsource)
    if min_elv>0:
        run_elvflag(cont_used,min_elv,logfile)
    runclcal(cont_used, 4, 7, 8, '', 1, refant)
    run_snplt(cont_used, inter_flag)

    mprint('######################',logfile)
    mprint(get_time(),logfile)
    mprint('######################',logfile)

###################################################################

if split_flag==1:

    split_sources=get_split_sources(cont_data, target, cvelsource, calsource)

    if cont_data2.exists() and line==cont:
        check_sncl(cont_data2, 4, 8,logfile)
        run_split(cont_data2, split_sources, split_outcl, doband, bpver)
    else:
        check_sncl(cont_data, 4, 8,logfile)
        run_split(cont_data, split_sources, split_outcl, doband, bpver)

    cvelsource = findcvelsource(line_data, cvelsource)

    if line_data2.exists():
        check_sncl(line_data2, 4, 8,logfile)
        run_masplit(line_data2, cvelsource, split_outcl, doband, bpver,smooth,channel)
    else:
        check_sncl(line_data, 4, 8,logfile)
        run_masplit(line_data, cvelsource, split_outcl, doband, bpver,smooth,channel)

if print_sn_flag==1:

    if line_data2.exists():
        outdata1 = line_data2
    else:
        outdata1 = line_data

    if outdata1.exists():
        outdata1.clrstat()
        runprtsn(outdata1)
        runprtsu(outdata1)

        mprint('######################',logfile)
        mprint(get_time(),logfile)
        mprint('######################',logfile)

if co_imagr_flag==1:   

    split_sources=get_split_sources(cont_data, target, cvelsource, calsource)

    for source in split_sources:
         split_data=AIPSUVData(source,split_outcl,cont_data.disk,1)
         if split_data.exists():
             runimagr(split_data, source, niter, cellsize, imsize, -1, imna,
                      antennas, uvwtfn, robust,beam,baselines=ant_bls,timer=imgr_timer)
             _zapbeam(source)

if ma_imagr_flag==1:

    cvelsource=findcvelsource(line_data, cvelsource)
    for source in cvelsource:
         split_data=AIPSUVData(source,split_outcl,line_data.disk,1)
         if split_data.exists():
             if len(split_data.name)>=12:
                 print split_data.name
                 split_data.rename(split_data.name[0:8],split_data.klass,split_data.seq)
             runmaimagr(split_data, source, niter, cellsize, imsize,channel,
                        -1, imna, uvwtfn, robust, beam,baselines=ant_bls)
             _zapbeam(split_data.name[0:8])

if cube_imagr_flag==1:

    cvelsource=findcvelsource(line_data, cvelsource)
    for source in cvelsource:
         split_data=AIPSUVData(source,split_outcl,line_data.disk,1)
         if split_data.exists():
             runcube(split_data, source, niter, cellsize, imsize,bchan,echan,
                     -1, antennas, uvwtfn, robust, beam)

if grid_flag==1:
    check_sncl(cont_data, 4, 8,logfile)
    run_grid(cont_data, gridsource, cellsize, imsize, n_grid, m_grid, grid_offset,uvwtfn, robust, beam)

if fittp_flag==1:

    split_sources=get_split_sources(cont_data, target, cvelsource, calsource)

    mprint('########################################################', logfile)
    for source in split_sources:
        run_fittp_data(source, split_outcl, defdisk, logfile)

    cvelsource = findcvelsource(line_data, cvelsource)
    for source in cvelsource:
        run_fittp_data(source, split_outcl, defdisk, logfile)
    mprint('########################################################', logfile)

###########################################################################

if imv_prep_flag==1:
    if line_data2.exists(): 
        line_data2.clrstat()
        check_sncl(line_data2, 4, 8,logfile)

    check_sncl(cont_data, 4, 8,logfile)
    calsource = findcal(line_data2, calsource)

    if not os.path.exists('./multiview'):
        os.makedirs('./multiview')

    if cont_data2.exists() and line!=cont: 
        cont_data2.zap()
    elif cont_data2.exists() and line==cont:
        mprint('########################################################', logfile)
        mprint('Deleting cont_data2, hope you did not want it',logfile)
        mprint('########################################################', logfile)
        cont_data2.zap()

    runsplat(cont_data,cont_data2,sources=target,stokes='I',docal=1)

    runcalib(cont_data2,docal=1,snver=1,solmode='P',soltype='L1R',aparm7=1,refant=refant)
    #runcalib(line_data2,docal=1,snver=5,solmode='P',soltype='L1R',aparm7=1,chan=channel)

    make_multiviewcontrol(cont_data2,line_data2,mv_window=mvwin)

    if os.path.exists('./multiview/multiview.TBOUT'): os.remove('./multiview/multiview.TBOUT')
    runtbout(cont_data2,'SN',1,'./multiview/multiview.TBOUT')
    replace('./multiview/multiview.TBOUT',"'INDE'",'-999.9')

    if os.path.exists('./multiview/target.TBOUT'): os.remove('./multiview/target.TBOUT')
    runtbout(line_data2,'SN',4,'./multiview/target.TBOUT')
    replace('./multiview/target.TBOUT',"'INDE'",'-999.9')

if imultiv_flag==1:
    if line_data2.exists(): 
        line_data2.clrstat()
        check_sncl(line_data2, 4, 8,logfile)

    calsource = findcal(line_data2, calsource)

    if not os.path.exists('./multiview'):
        sys.exit('./multiview directory does not exist. Run with imv_prep_flag first.')

    os.chdir('./multiview/')
    mv_task = os.system(mv_path+'multiview_4x > multiview_4x.prt')
    if mv_task==0:
        print 'MULTV: Appears to have ended successfully'
    else: 
        print 'MUTLV: Purports to die of UNNATURAL causes'
        print 'Please check input files in ./multiview/'
        raise ValueError

    os.chdir('../')

if imv_app_flag==1:
    if line_data2.exists(): 
        line_data2.clrstat()
        check_sncl(line_data2, 4, 8,logfile)
        
    runtbin(line_data2,'./multiview/target.TBIN')
    runclcal(line_data2,5,8,9,'',1,refant)

    mprint('######################',logfile)
    mprint(get_time(),logfile)
    mprint('######################',logfile)

if imv_imagr_flag==1:
    if line_data2.exists(): 
        line_data2.clrstat()
        check_sncl(line_data2, 5, 9,logfile)

    mprint('######################',logfile)
    mprint('Imaging channel: '+str(channel),logfile)
    mprint('######################',logfile)

    runmaimagr(line_data2,calsource,niter,cellsize,imsize,channel,1,'',
        uvwtfn,robust,beam,baselines=ant_bls,timer=imgr_timer)
    _zapbeam(calsource,inseq=1)
    _zapbeam(calsource,inseq=2)
    mprint('######################',logfile)
    mprint(get_time(),logfile)
    mprint('######################',logfile)

###########################################################################

if phase_cal_flag==1:
    source = calsource
    indata = AIPSUVData(source,split_outcl,line_data.disk,1)
    cal_ants = antennas
    im_ants = []
    for entry in dofit[0]:
        im_ants.append(-entry)

    for i in range(len(phase_loop)):
        cal_data = phase_selfcal(indata, source, phase_loop[i], line_data.disk,
                                 niter, cellsize,imsize,logfile, imna,
                                 cal_ants, im_ants, refant, fr_image, beam)
        indata = cal_data
    mprint('########################################################', logfile)

if amp_cal_flag==1:
    source = calsource
    if len(dofit)!=len(amp_loop):
        dofit=range(amp_loop)
        mprint('dofit and amp_loop not equal length, solving for all antennas', logfile)
        for i in range(len(dofit)):
            dofit[i]=0

    indata = AIPSUVData(source,split_outcl,line_data.disk,1)
    indata = check_cal(indata, source, outdisk, logfile, imna)
    for i in range(len(amp_loop)):
        cal_data = amp_selfcal(indata, source, amp_loop[i], line_data.disk,
                               niter, cellsize,imsize,logfile, imna, antennas,
                               refant, dofit[i], beam)
        indata = cal_data
    mprint('########################################################', logfile)

if refeed_flag==1:

    source=calsource
    if imna=='':
        outname = source
    else:
        outname = source[:11-len(imna)]+imna
    nimage = 1
    while AIPSImage(outname,'ICL001',line_data.disk, nimage+1).exists():
        nimage+=1

    model  = AIPSImage(outname,'ICL001',line_data.disk, nimage)


    if cont_data2.exists():
        mprint('Using shifted data (CVEL).',logfile)
        cont_used=cont_data2
    else:
        mprint('Using unshifted data (CVEL).',logfile)
        cont_used=cont_data

    check_sncl(cont_used, 3, 7,logfile)
    fringecal(cont_used,model,nmaps,refant,calsource,solint,smodel,doband,bpver,dpfour)
    if snflg_flag==1:
        runsnflg(cont_used, 4, calsource)
    if min_elv>0:
        run_elvflag(cont_used,min_elv,logfile)
    runclcal(cont_used, 4, 7, 8, '', 1, refant)
    run_snplt(cont_used, inter_flag)

    mprint('######################',logfile)
    mprint(get_time(),logfile)
    mprint('######################',logfile)

    split_sources=get_split_sources(cont_data, target, cvelsource, calsource)

    if cont_data2.exists():
        check_sncl(cont_data2, 4, 8,logfile)
        run_split(cont_data2, split_sources, split_outcl, doband, bpver)
    else:
        check_sncl(cont_data, 4, 8,logfile)
        run_split(cont_data, split_sources, split_outcl, doband, bpver)

    cvelsource = findcvelsource(line_data, cvelsource)

    if line_data2.exists():
        check_sncl(line_data2, 4, 8,logfile)
        run_masplit(line_data2, cvelsource, split_outcl, doband, bpver,smooth,channel)
    else:
        check_sncl(line_data, 4, 8,logfile)
        run_masplit(line_data, cvelsource, split_outcl,smooth,channel)


    split_sources=get_split_sources(cont_data, target, cvelsource, calsource)

    for source in split_sources:
         split_data=AIPSUVData(source,split_outcl,cont_data.disk,1)
         if split_data.exists():
             runimagr(split_data, source, niter, cellsize, imsize, -1, imna,
                      antennas,uvwtfn, robust, beam, baselines=ant_bls,timer=imgr_timer)

if phase_target_flag!='':
    source = phase_target_flag
    indata = AIPSUVData(source,split_outcl,line_data.disk,1)
    for entry in phase_loop:
        cal_data = phase_selfcal(indata, source, entry, line_data.disk, niter,
                                 cellsize,imsize,logfile, imna, antennas,
                                 refant, beam)
        indata = cal_data

if amp_target_flag!='':
    source = phase_target_flag
    indata = AIPSUVData(source,split_outcl,line_data.disk,1)
    for entry in amp_loop:
        cal_data = amp_selfcal(indata, source, entry, line_data.disk, niter,
                               cellsize,imsize,logfile, imna, antennas,
                               beam, refant)
        indata = cal_data

if apply_selfcal==1:
    do_apply_selfcal(calsource,calsource,split_outcl,defdisk, refant)
    for entry in target:
        do_apply_selfcal(entry,calsource,split_outcl,defdisk, refant)

if plot_tables!=-1:

    if plot_tables==line:
        if line_data2.exists():
            mprint('Using shifted data (CVEL).',logfile)
            line_used = line_data2
        else:
            mprint('Using unshifted data (CVEL).',logfile)
            line_used = line_data
        plot_data=line_used
    else:
        plot_data=data[plot_tables]

    setup_plotfiles(plot_data)
    run_snplt_2(plot_data, 1,  'AMP', 'SN1')
    run_snplt_2(plot_data, 2,  'AMP', 'SN2')
    run_snplt_2(plot_data, 3, 'DELA', 'SN3')
    run_snplt_2(plot_data, 4, 'PHAS', 'SN4')
    if (os.path.exists('delay-rate.ps')):
        os.popen(r'convert delay-rate.ps delay-rate.png')
        os.popen(r'mv delay-rate.png plotfiles/')
    if (os.path.exists('delay-rate-ionos.ps')):
        os.popen(r'convert delay-rate-ionos.ps delay-rate-ionos.png')
        os.popen(r'mv delay-rate-ionos.png plotfiles/')
#    rdbe_plot(data[geo_data_nr],logfile,'GEO')

if rpossm_flag==1:
    cvelsource=findcvelsource(line_data, cvelsource)
    for source in cvelsource:
         split_data=AIPSUVData(source,split_outcl,line_data.disk,1)
         if split_data.exists():
            tv=AIPSTV.AIPSTV()
            if inter_flag==1:
                if tv.exists()==False: tv.start()
            runrpossm(split_data,cvelsource,tv,inter_flag,antennas)

if ma_sad_flag==1:
    cvelsource=findcvelsource(line_data, cvelsource)
    for source in cvelsource:
        ma_cube=AIPSImage(source,'ICL001',line_data.disk,1)
        run_ma_sad(ma_cube,line_data,min_snr,dyna)
#        orfit_to_plot_sorted(ma_cube, line_data)
        orfit_to_plot_sorted(ma_cube, line_data,bchan,vel)

if plot_map==1:
    cvelsource=findcvelsource(line_data, cvelsource)
    for source in cvelsource:
        ma_cube=AIPSImage(source,'ICL001',line_data.disk,1)
        plot_spot_map(ma_cube,line_data)

tv=AIPSTV.AIPSTV()
if tv.exists():
    raw_input('Press Enter to close window. ')
    tv.kill()

# End main script
##############################################################################
