from AIPSData import AIPSUVData as AUV
from AIPS import AIPS
from AIPSTask import AIPSTask
from Wizardry.AIPSData import AIPSUVData as WAUV
from Wizardry.AIPSData import AIPSImage  as WAIm
from pylab import *

import peakutils

import sys, os, random, copy, difflib, re, pdb, time
import numpy as np
from math import exp, log, pi, atan2

##############################################################################
##############################################################################

''' 
Version = 3.0

This script scales antenna gains from AIPS off maser autocorrections.
A single antenna (refant) is assumed to be perfect in pol[0].
Polynomial-<n> baselines are fit to each antenna individually in n.
Corrections are reapplied to principal data back to most recent CL table.

This script assumes (less than v2):
    - single IF data and a STRONG maser correlated at high spectral
           resolution;
    - that there has been accurate pre-amp calibration applied to refant
    - Data comes from klass = UVDATA;

If these assumptions aren't true fix them yourself or write a different script
Lucas J. Hyland 2019
Lucas.Hyland@utas.edu.au

Usage:

parseltongue ./amplitude_autocorrect_v3.py source refant

e.g.

parseltongue ./amplitude_autocorrect_v3.py G123.56+0.7 2 

'''

# G123.56+0.7  2       
# argv[1]      argv[2]                                      

##############################################################################
##############################################################################

kla = 'UVDATA'
seq = 1
dsk = 1

##############################################################################
##############################################################################

def main():
    ''' ####################
        ### Help section ###
        #################### '''
    if True in [sys.argv[1] in ['help','-help','--help','h','-h','--h'],len(sys.argv)<3]:
        #######1234567890123456789012345678901234567891234567
        print '#################### HELP ####################' 
        print '#                                            #'
        print '#  Format: parseltongue <script>.py source   #'
        print '#                                   refant   #'
        print '#                                            #'
        print '################# PARAMETERS #################'     
        print '#                                            #'
        print '#  source: Name of maser in AIPS             #'
        print '#                                            #'
        print '#  refant: AIPS antenna number  of station   #'
        print '#          which  has a good detection of    #'
        print '#          maser in POL1 and accurate pre-   #'
        print '#          calibration applied               #'
        print '#                                            #'
        print '################## END HELP ##################'  
        sys.exit()
    ''' ############################################
        ### Tentatively check arguments are okay ###
        ############################################ '''

    ### Get reasonable experiment name ###
    logfile = open(sys.argv[0].replace('.py','.log'),'a')
    AIPS.userno,exp=get_experiment()
    indata = AUV(exp,kla,dsk,seq)
    if not indata.exists():
        sys.exit(indata+' in AIPSID=\
            '+str(AIPS.userno)+' does not exist')
    else: print str(indata)+' does exist'
    ### Make sure we operate on CL1 ###
    CL = indata.table_highver('CL')
    if CL==1: sys.exit('CL1 !! Need to pre-calibrate!')
    ### Check if source is valid ###
    source   = sys.argv[1]
    sources  = indata.sources
    msources = difflib.get_close_matches(source, sources)
    try:
        if not msources[0]==source:
            read = raw_input('Did you mean '+msources[0]+'? [y|n]\n')
            if not read in ['Y','y']: sys.exit('Unknown source')
            else: source = msources[0]
    except IndexError:
        sys.exit('No suggested sources, unknown source '+source)
    ### Check if ref antenna is valid ###
    reference_antenna = int(sys.argv[2])
    antable = indata.table('AN',1)
    an = {}
    for row in antable:
        an.update({row['nosta']:row['anname'].replace(' ','')})
    if not reference_antenna in an:
        print 'Antenna table is '+str(an)
        reference_antenna = raw_input('Please enter\
    antenna number not name\n')
    ''' ###################################
        ### Get relevant visbility data ###
        ################################### '''
    if len(source)>12: sname = source[:12] #AIPS innames can't handle>12
    else: sname=source
    spl_data = AUV(sname,'TMPL',1,69)
    if spl_data.exists(): spl_data.zap()
    #sys.exit()
    _split(indata,source)
    wdata = WAUV(spl_data.name, spl_data.klass, 1, 69)
    t = []
    for vis in wdata:
        t += [vis.time]
    nchans = wdata.header['naxis'][2]
    npolas = wdata.header['naxis'][1]
    ncmplx = wdata.header['naxis'][0]
    data = np.ones(shape=(len(t), nchans, npolas, ncmplx))
    # create key dict to map time -> ant
    antkey = {}
    for i in range(len(wdata.antennas)):
        antkey.update({i:[]})
    i=0
    for vis in wdata:
        data[i] = vis.visibility
        antenna = vis.baseline[0]
        antkey[antenna-1]+=[i]
        i=i+1
    ''' #######################################
        ### Calculate and apply corrections ###
        ####################################### '''  
    #find peaks from ACOR spectrum (avg time,ant,pol)
    #no weighting or scaling applied
    peaks = find_peaks(data[:,:,:,0].mean(0).mean(1))
    S_peak = {}
    # calculating peak for each ant/polar (avg time)
    for i in range(len(indata.antennas)):
        S_peak.update({i:[]})
        S_debaselined = np.ones(shape=data[antkey[i],:,:,0].mean(0).shape)
        for j in range(data[antkey[i],:,:,0].mean(0).shape[1]):
            S_debaselined[:,j] = baselined_spectrum(data[antkey[i],:,j,0].mean(0))
        S_peak[i] = S_debaselined[peaks,:].T
    S_ref = S_peak[reference_antenna-1][0]
    correctionsq = np.zeros(shape=(len(indata.antennas),npolas,len(peaks)))
    for i in range(len(indata.antennas)):
        correctionsq[i] = S_ref/S_peak[i]
    corrections = correctionsq**0.5
    polz = {}
    for p in range(len(indata.polarizations)):
        polz.update({p:indata.polarizations[p]})
    for ant in range(len(indata.antennas)):
        for pol in range(len(indata.polarizations)):
            cor = round(median(corrections[ant,pol]),2)
            if True in [cor>=1.0, cor<=1.0]:
                ''
            else:
                cor=0.0
            print >> logfile, "Correction for "+an[ant+1]+' '+polz[pol]+' is '+str(cor)
            _clcor(indata,ant+1,polz[pol],cor,CL)
    spl_data.zap()
    print 'CLCOR2: localhos 31DEC16 TST: Cpu=      0.0  Real=      0  IO=         1'
    print '                                        '
    print '########################################################################'
    print 'AMPL: Appears to have ended successfully'
    print '########################################################################'
    print '                                        '
    print 'If you want please check '+sys.argv[0].replace('.py','.log')+' for list of applied corrections'
    print 'And <antenna>_amplitude_tmp.png diagnostic plots (coming soon)'
    ''' ### DONE ### '''

##############################################################################
# Spectra fitting functions
##############################################################################

def calculate_ideal_baseline(np_array):
    p = []
    for i in range(20):
        b = peakutils.baseline(np_array,i)
        res = np_array - b
        p += [calcProbability_conservative_np(res)]
    best = np.array(p).argmax()+1
    baseline = peakutils.baseline(np_array,best)
    return baseline

def calcProbability_conservative_np(res):
    probability = []
    for i in range(res.shape[0]):
        residual = res[i]
        probability_list = np.log ( ( 1 - np.exp(-1.0/2*(residual**2 + 2e-15) ) ) ) - np.log ( residual**2 + 1e-15)
        probability += [probability_list.sum()]
    PR = np.ma.masked_invalid(probability).sum()
    return  PR

def calculate_peaks(spectrum):
    try:
        baseline  = calculate_ideal_baseline(spectrum)
    except ValueError:
        return 0.0
    return (spectrum-baseline)[peakutils.indexes(spectrum-baseline,\
     thres = 20.0/((spectrum-baseline)/(std_approx(spectrum\
        - baseline))).max())]

def find_peaks(spectrum):
    try:
        baseline  = calculate_ideal_baseline(spectrum)
    except ValueError:
        return 0.0
    return peakutils.indexes(spectrum-baseline,\
     thres = 20.0/((spectrum-baseline)/(std_approx(spectrum\
        - baseline))).max())

def plot_peaks(spectrum):
    try:
        baseline  = calculate_ideal_baseline(spectrum)
    except ValueError:
        return np.array([0])
    x = np.array(range(len(spectrum)))
    index = peakutils.indexes(spectrum-baseline, thres = 10.0/((spectrum-baseline)/(std_approx(spectrum-baseline))))
    print index, spectrum[index]
    plt.plot(x, spectrum-baseline, x[index], (spectrum-baseline)[index], 'ro')
    return

def baselined_spectrum(spectrum):
    try:
        baseline  = calculate_ideal_baseline(spectrum)
    except ValueError:
        return np.array([0])
    return spectrum-baseline

##############################################################################
# AIPS-PARSELTONGUE FUNCTIONS
##############################################################################

def get_vis(split_data,chan):
    ### Print and load in data ###
    ants = range(len(split_data.antennas))
    npols = split_data.header['naxis'][1]
    nchan = split_data.header['naxis'][2]    
    nifs  = split_data.header['naxis'][3]
    if nifs>1:
        sys.exit('Can only deal with one IF')
    if npols==4:
        print 'Assuming full-stokes RR,LL,RL,LR'
    elif npols==1:
        print 'Assuming single-pol '+str(pols[0])
    else:
        sys.exit('Unknown polarizations, tried full/single')
    _prtuv(split_data,chan)
    f = splitt(get_file('/tmp/prtuv.tmp'))
    #pdb.set_trace()
    F = np.array(refine_fits(f,18,'-',3,None))
    t = F[:,1]
    ### Make data structure to fit data into ###
    L = np.zeros(shape=(len(t),len(ants),npols))
    ### Buffer values to take care of empties ###
    ant_buff = {}
    for i in ants:
        ant_p = {}
        for p in range(npols):
            ant_p.update({p:0.0})
        ant_buff.update({i:ant_p})
    for k in range(len(ant_buff)):
        for j in range(F.shape[0]):
            if F[j,3]==str(int(k+1)):
                ant_buff[k].update({0:float(F[j,6]),\
                    1:float(F[j,9]),2:float(F[j,12]),3:float(F[j,15])})
                break
    ### Format data into data structure ###
    for i in range(len(t)):
        for k in range(len(ant_buff)):
            if F[i,3]==str(int(k+1)):
                ant_buff[k].update({0:float(F[i,6]),\
                    1:float(F[i,9]),2:float(F[i,12]),3:float(F[i,15])})
                break
        for k in range(len(ant_buff)):
            for p in range(npols):
                L[i,k,p]=ant_buff[k][p]*100
    return t, L

def _clcor(data,ant,pol,clcorparm,cltable):
    clcor=AIPSTask('clcor')
    clcor.default
    clcor.indata=data
    clcor.stokes=pol
    clcor.gainu=cltable
    clcor.gainv=cltable
    clcor.clcorprm[1]=round(clcorparm,2)
    clcor.opcode='GAIN'
    clcor.antenna[1]=ant
    clcor.bif=0
    clcor.eif=0
    #clcor.log=open('./amplitude_autocorrect_v3.log', 'a')
    clcor.go()

def _split(data,source):
    spli=AIPSTask('split')
    spli.default
    spli.indata=data
    spli.source[1]=source
    spli.docalib=1
    spli.gainu=0
    spli.aparm[5]=2
    spli.aparm[6]=1
    spli.outclass='TMPL'
    spli.outseq=69
    spli.flagver=1    
    spli.go()

def _prtuv(split_data,channel):
    prtuv=AIPSTask('prtuv')
    prtuv.default
    prtuv.channel=channel
    prtuv.indata=split_data
    prtuv.dparm[7]=1
    prtuv.cparm[10]=0.01
    prtuv.docrt=-1
    removefile('/tmp/prtuv.tmp')
    prtuv.outprint='/tmp/prtuv.tmp'  
    prtuv.go()

def get_experiment():
    files = os.listdir('./')
    reffile = []
    for i in range(len(files)):
        if "LST" in files[i]:
            reffile += [files[i]]
    if reffile==[]:
        print 'No AIPS <EXP>.LST  file to read from'
        print "You should make them, they're useful"
        inp = raw_input('USERID INNAME ?\n')
        aipsid     = inp.replace(',',' ').split()
        experiment = inp.replace(',',' ').split()
    else:
        if len(reffile)>1:
            print ' '
            print 'Multiple experiments detected, please choose N'
            print 'n : EXP'
            for i in range(len(reffile)): print str(i+1)+' : '+reffile[i]
            choice = raw_input('')
            try:
                index  = int(choice)-1
            except IndexError:
                sys.exit('What')
            readfile = get_file(reffile[index])
        else:
            readfile = get_file(reffile[0])
        aipsid   = readfile[0].split()[2]
        experiment = readfile[1].split()[2]
    return int(aipsid), experiment

def refine_fits(old_list, length_cat, symbol_cat,\
             column_search_number, column_retrival_number):
    #to refine an imported list after the use of get_file()
    #to remove header information or filter the list
    # only works for 2-dim type objects
    #symbol_cat is a single character
    tmp_list=[]
    if not length_cat==None:
        for i in old_list:
            if len(i)==length_cat:
                tmp_list+=[i]
    else:
        tmp_list=old_list
    tmp_list_2=[]
    if not symbol_cat==None:
        for i in tmp_list:
            if not column_search_number==None:
                column=i[column_search_number-1] #searches specific column for match
            else:
                column=i #searches all columns.
            for row_element in column: #searching rows in columns
                if row_element.count(symbol_cat)>0:
                    tmp_list_2+=[i]
                    break #ends if it finds it to prevent line repeats
                else:
                    continue #continues to look if it doesn't
    else:
        tmp_list_2=tmp_list
    tmp_list_3=[]
    if column_search_number==None:
        if not column_retrival_number==None:
            for i in tmp_list_2:
                tmp_list_3+=[i[column_retrival_number-1]]
        else:
            tmp_list_3=tmp_list_2
    else:
        tmp_list_3=tmp_list_2
    tmp_list_4 = []
    for k in range(len(tmp_list_3)):
        if 'localhos' not in tmp_list_3[k]:
            tmp_list_4+=[tmp_list_3[k]]
    return tmp_list_4

##############################################################################
# GENERAL-use FUNCTIONS
##############################################################################

def axis_median3D(np_array3D,axis):
    s = np_array3D.shape
    np2Dp = np.zeros(shape=s[1:])
    for i in range(np2Dp.shape[0]):
        for j in range(np2Dp.shape[1]):
            np2Dp[i,j]=np.percentile(np_array3D[:,i,j],50)
    return np2Dp

def axis_median2D(np_array2D):
    # only does median along axis 1
    s = np_array2D.shape
    np1Dp = np.zeros(shape=s[0])
    for i in range(np1Dp.shape[0]):
        np1Dp[i]=np.percentile(np_array2D[i,0],50)
    return np1Dp

def median_smooth(x,M,solint):
    if False in [type(x)==np.ndarray,type(M)==np.ndarray,len(x)==len(M)]:
        sys.exit('Requires np.numpy arrays of same primary length')
    # use lists due to unknown output length
    Y = []
    X = []
    count = 0
    j = 0
    while j<len(x):
        j = j + count
        count = 0
        y     = []
        for i in range(len(x)-j):
            if x[j+i]-x[j]<solint:
                count=count+1
                y+=[M[i+j]]
            else:
                Y+=[axis_median(np.array(y))]
                X+=[0.5*(x[j]+x[j+i])]
                break
    # reshape into numpy objects along 0th axis
    m = np.zeros(shape=(len(Y),Y[0].shape[0],Y[0].shape[1]))
    t = np.zeros(len(X))
    for k in range(len(Y)):
        m[k] = Y[k]
        t[k] = X[k]
    return t, m

def average_repeats(x, M):
    '''
    Finds repeats in x, replaces values 
    then averages values in M
    '''
    if False in [type(x)==np.ndarray,type(M)==np.ndarray,len(x)==len(M)]:
        sys.exit('Requires np.numpy arrays of same primary length')
    # use lists due to unknown output length
    Y = []
    X = []
    count = 0
    j = 0
    while j<len(x):
        j = j + count
        count = 0
        y     = []
        for i in range(len(x)-j):
            if x[j]==x[j+i]:
                count=count+1
                y+=[M[i+j]]
            else:
                Y+=[axis_median(np.array(y))]
                X+=[x[j]]
                break
    # reshape into numpy objects along 0th axis
    m = np.zeros(shape=(len(Y),Y[0].shape[0],Y[0].shape[1]))
    t = np.zeros(len(X))
    for k in range(len(Y)):
        m[k] = Y[k]
        t[k] = X[k]
    return t, m

def locate_index(reference_t,T):
    ind = abs(T - reference_t).argmin()
    return ind

def median(data):
    average=np.percentile(data,50)
    return average

def std_approx(data):
    approx_std = 0.5*(np.percentile(data,84)-np.percentile(data,16))
    return approx_std

def make_int(list):
    for i in range(len(list)):
        try:
            list[i]=int(list[i])
        except ValueError:
            try:
                list[i]=list[i-1]
            except ValueError:
                list[i]=list[i+1]
    return list

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
    return new_list

def convertdhms2fday(dhms):
    return float('%1.3f'%(dhms[0]+dhms[1]/24.0+
            dhms[2]/(24.0*60.0)+
            dhms[3]/(24.0*60.0*60.0)))

def convertfday2dhms(frac_day):
    D = int(frac_day)
    H = int(24*(frac_day - D))
    M = int(60*(24*(frac_day - D) - H))
    S = int(60*(60*(24*(frac_day - D) - H) - M))
    return D, H, M, S

def convertstrdate2dhms(str_date):
    return make_int(str_date.replace('/',':').split(':'))

def convertdhms2strdate(dhms):
    return '%s/%2.2d:%2.2d:%2.2d' % (dhms[0],dhms[1],dhms[2],dhms[3])

def removefile(dirpath):
    if os.path.exists(dirpath):
        os.remove(dirpath)

##############################################################################
# MAIN
##############################################################################

if __name__=='__main__':
    main()

##############################################################################
##############################################################################