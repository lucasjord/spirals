#!/home/observer/anaconda2/bin/ParselTongue

from AIPSData import AIPSUVData
from AIPS import AIPS
from AIPSTask import AIPSTask
import AIPS, os, math, time, sys, os, random
import copy, difflib, re, numpy as np, argparse
from math import exp, log, pi, atan2
from numpy.linalg import inv

import pdb

##############################################################################
##############################################################################

''' 
Version 3.0
Written by Lucas Jordan Hyland
University of Tasmania
Date 19/10/21
'''

def main():
    print '\n'

    parser = argparse.ArgumentParser()
    parser.add_argument("epoch",
                        help="Experiment code e.g. s006i",
                        type=str)
    parser.add_argument("aipsID",
                        help="AIPS ID, e.g 666",
                        type=int)
    parser.add_argument("-c","--cltable",
                        help="CL table to save corrections to. Default 6",
                        type=int,default=6)
    parser.add_argument("-D","--disk",
                        help="Disk data is located. Default 1.",
                        type=int,default=1)
    parser.add_argument("-q","--seq",
                        help="Sequence of continuum data. Default 1",
                        type=int,default=1)
    parser.add_argument('-n','--niter',
                        help='Number of iterations to do. Default 3',
                        type=int,default=3)
    parser.add_argument('-s','--source',
                        help='Calibrator source to match. Default determine automatically.',
                        type=str, default=None)
    parser.add_argument('-f','--flux',
                        help='Flux of the calibrator if known/suspected (Jy). Default 1.0',
                        type=float,default=1.0)
    parser.add_argument('-v','--verbose',
                        help='Increase verbosity.',
                        action='count',default=0)

    args = parser.parse_args()

    kla = 'UVDATA'
    seq = args.seq
    dsk = args.disk 

    AIPS.userno = args.aipsID
    exp         = args.epoch.upper() + '_C'    #get_experiment()
    indata      = AIPSUVData(exp,kla,dsk,seq)
    npols       = len(indata.polarizations)
    if not indata.exists():
        sys.exit(indata+' in AIPSID='+str(AIPS.userno)+' does not exist')
    else:
        print str(indata)+' does exist'
    while indata.table_highver('CL')>args.cltable:
        tab = indata.table_highver('CL')
        indata.zap_table('CL',tab)
    # redirect stdout to null while doing aips stuff
    old_stdout = sys.stdout
    null = open(os.devnull,'w')
    n = 0
    print ''
    print '{:3s} {:20s}'.format('Ant:','Correction factors')
    while n<=args.niter:
        sys.stdout = null
        get_plots(indata)
        indata.zap_table('PL',-1)
        if not 'indata' in globals(): globals().update({'indata':indata})
        if    npols==1: W, baselines_key = read_in_data_1pol(indata)
        elif  npols==2: W, baselines_key = read_in_data_2pol(indata)
        else: sys.exit('Unknown polarization configuration')
        polz = {}
        for i in range(npols):
            polz.update({i:indata.polarizations[i]})
        ifs = indata.header['naxis'][3]
        ant_key, epsilon = calculate_offsets(indata,W,baselines_key,args)
        clcorprm = (epsilon + 1).tolist()
        for ant in range(len(ant_key)):
            for freq in range(ifs):
                for pol in range(npols):
                    if not kla=='SPLIT':
                        clcor(data=indata,bif=freq+1,ant=int(ant_key[ant]+1),
                            clcorparm=clcorprm[ant][freq+pol*ifs],polar=polz[pol],cltable=args.cltable)
                    else:
                        sncor(indata,freq+1,int(ant_key[ant]+1),clcorprm[ant][freq+pol*ifs],polz[pol])
        #
        sys.stdout = old_stdout
        #
        for ant in range(len(ant_key)):
            print '{:3.0f}:'.format(ant+1),
            for freq in range(ifs):
                for pol in range(npols):
                    print '{:4.2f}'.format(clcorprm[ant][freq+pol*ifs]),
            print ''    
        print ''
        n = n + 1
    #
    s, tm = get_best_sources(indata)
    print '#################################'
    print ' Please check possm on '+s
    print ' for timer ',
    for i in tm: print str(i),
    print ''
    print '    If amp is not flat, rerun    '
    print ' Else, now run manual phase cal  '
    print '#################################'

##############################################################################
##############################################################################

def get_plots(data):
    source,time=get_best_sources(data)
    possm=AIPSTask('possm')
    possm.default
    possm.indata=data
    data.clrstat()
    possm.nplots=0
    possm.stokes='HALF'
    possm.aparm[1:]=[0,1,0,0,0,0,0,0,1, 0]
    possm.flagver=0
    possm.sources[1]=source
    possm.timer[1:]=time
    possm.solint=0
    possm.codetype='AMP'
    possm.gainu=6
    possm.docalib=1
    possm.dotv=-1
    possm.bchan=8
    possm.echan=24
    #possm.smooth[1:]=[16,16,16]
    possm.outtext='/tmp/possm.tmp'
    removefile(possm.outtext)
    perm = list(product([0,1],repeat=len(data.antennas)))
    #possm.grchan=1
    final_perm=[]
    for i in range(len(perm)):
        if perm[i].count(1)==2:
            final_perm+=[list(perm[i])]
    for j in range(len(final_perm)):
        possm.antennas[1:]=make_int(list(np.array(final_perm[j])*(np.array(range(len(data.antennas)))+1)))
        possm.baseline    = possm.antennas
        try:
            possm()
        except RuntimeError:
            print 'No '+str(possm.baseline)+' data'
            continue

def read_in_data_1pol(indata):
    possmdata = split_block(splitt(get_file('/tmp/possm.tmp')),'Header')
    baselines = {}
    averaged_data = range(len(possmdata))
    for k in range(len(possmdata)):
        baselinefits = possmdata[k]
        header = baselinefits[0:14]
        dataR   = np.matrix(refine_fits(baselinefits,7,'R',3,None))
        dataL   = np.matrix(refine_fits(baselinefits,7,'L',3,None))
        if dataR.size==0 and dataL.size!=0:
            data = copy.copy(dataL)
        elif dataR.size!=0 and dataL.size==0:
            data = copy.copy(dataR)
        else:
            raise IndexError('Problem with data')
            sys.exit()
        number_of_if = indata.header['naxis'][3]
        number_chans = int(header[5][2])
        a = header[7][2]  # antenna 1
        b = header[8][2]  # antenna 2
        averaged_data[k] = list(np.ones(number_of_if)*99999)
        for n in range(number_of_if):
            if_data = data[n*number_chans:(n+1)*number_chans]
            try:
                t = int(if_data[:,1][0].tostring().replace('\x00',''))-1
            except IndexError:
                continue
            averaged_data[k][t] = if_data[:,5].astype(np.float).mean()
        baselines.update({k:[a,b]})
    W = np.matrix(averaged_data)    # flux matrix
    return W, baselines

def read_in_data_2pol(indata):
    possmdata = split_block(splitt(get_file('/tmp/possm.tmp')),'Header')
    baselines = {}
    averaged_data = range(len(possmdata)/2)
    for k in range(len(averaged_data)):
        baselinefits1 = possmdata[2*k]
        baselinefits2 = possmdata[2*k+1]
        header = baselinefits1[0:14]
        header2 = baselinefits2[0:14]
        if not header==header2:
            sys.exit('Consecutive headers not equal')
        dataR   = refine_fits(baselinefits1,7,'R',3,None)+refine_fits(
            baselinefits2,7,'R',3,None)
        dataL   = refine_fits(baselinefits1,7,'L',3,None)+refine_fits(
            baselinefits2,7,'L',3,None)
        number_of_if = indata.header['naxis'][3]*2
        number_chans = int(header[5][2])
        a = header[7][2]  # antenna 1
        b = header[8][2]  # antenna 2
        averaged_data[k] = list(np.ones(number_of_if)*99999)
        data = np.matrix(dataR + dataL)
        for n in range(number_of_if):
            if_data = data[n*number_chans:(n+1)*number_chans]
            try:
                t = int(if_data[:,1][0].tostring().replace('\x00',''))-1
            except IndexError:
                continue
            averaged_data[k][n] = if_data[:,5].astype(np.float).mean()
        baselines.update({k:[a,b]})
    W = np.matrix(averaged_data)    # flux matrix
    return W, baselines

def calculate_offsets(data,W_matrix,antenna_key,args):
    ant_map = np.array(range(len(data.antennas)))
    m = np.zeros(shape=(len(antenna_key),len(data.antennas))).tolist()
    for i in range(len(antenna_key)):
        for a in antenna_key[i]:
            m[i][int(a)-1]=1
    # removing antennas with no data from array and mapping
    M = np.array(m)
    s = copy.copy(M)
    for i in range(s.shape[1]):
        j=s.shape[1]-1-i
        if s.sum(0)[j]==0.0:
            M = np.delete(M,(j),axis=1)
            ant_map = np.delete(ant_map,(j),axis=0)
    M = np.matrix(M)
    P = M.T*M
    try: 
        tmpP = inv(P)
    except LinAlgError:
        sys.exit('Offset Matrix is singular, cannot invert')
    S = args.flux #match_source(data)
    reweight = 0
    #look for missing or flagged data via IF (doesn't work right now)
    if reweight==1:
        for i in range(W_matrix.shape[1]):
            #pdb.set_trace()
            for j in range(W_matrix.shape[0]):   #over baselines
                b=-1
                if W_matrix[j,i]>99998.0:
                    diagnostic = M.T*W_matrix[:,i]
                    b = diagnostic.argmax()
                    break
                else:
                    continue
            if b!=-1:
                ants = range(4)
                ants.remove(b)
                for a in ants:
                    p_aj=[]
                    for j in range(W_matrix.shape[0]):
                        if m[j][a]==0 and m[j][b]==0:
                            p_cd = j
                        elif m[j][a]==1 and m[j][b]==0:
                            p_aj += [j]
                        elif m[j][a]==1 and m[j][b]==1:
                            p_ab = j
                    W_matrix[p_ab,i]=-3/(2/W_matrix[p_aj[0],i]+2/W_matrix[p_aj[1],i]-7/W_matrix[p_cd,i])

    D  = S/W_matrix -1 
    x  = (inv(P)*M.T)*D
    return ant_map, x

def clcor(data=None,bif=0,ant=0,clcorparm=0,polar='R',cltable=0):
    clcor=AIPSTask('clcor')
    clcor.default
    clcor.indata=data
    clcor.stokes=polar
    clcor.gainu=cltable
    clcor.gainv=cltable
    clcor.clcorprm[1]=clcorparm
    clcor.opcode='GAIN'
    clcor.antenna[1]=ant
    clcor.bif=bif
    clcor.eif=bif
    clcor.go()

def sncor(data,bif,ant,clcorparm,polar):
    sncor=AIPSTask('sncor')
    sncor.default
    sncor.indata=data
    sncor.snver=1
    sncor.stokes=polar
    sncor.sncorprm[1]=clcorparm
    sncor.opcode='MULA'
    sncor.antenna[1]=ant
    sncor.bif=bif
    sncor.eif=bif
    sncor.go()

##############################################################################
##############################################################################

def remove_repeats(list):
    track_count=[]
    for entry in list:
        if entry in track_count:
            continue
        else:
            track_count+=[entry]
    return track_count

def permutations(iterable, r=None):
    # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
    # permutations(range(3)) --> 012 021 102 120 201 210
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    if r > n:
        return
    indices = list(range(n))
    cycles = list(range(n, n-r, -1))
    yield tuple(pool[i] for i in indices[:r])
    while n:
        for i in reversed(range(r)):
            cycles[i] -= 1
            if cycles[i] == 0:
                indices[i:] = indices[i+1:] + indices[i:i+1]
                cycles[i] = n - i
            else:
                j = cycles[i]
                indices[i], indices[-j] = indices[-j], indices[i]
                yield tuple(pool[i] for i in indices[:r])
                break
        else:
            return

def product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

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

def get_experiment():
    files = os.listdir('./')
    ref = []
    for i in range(len(files)):
        if "LST" in files[i]:
            ref += [files[i]]
    if ref==[]:
        sys.exit('No AIPS output file to read from')
    reffile = [s for s in ref if "_C" in s]
    if len(reffile)>1:
        print 'Multiple possibilities detected, please indicate which one'
        print reffile
        for i in range(len(reffile)): print str(i+1),
        choice = raw_input('??\n')
        try:
            index  = int(choice)-1
        except IndexError:
            sys.exit('What')
        readfile = get_file(reffile[index])
    else:
        readfile = get_file(reffile[0])
    aipsid   = readfile[1].split()[-1].strip('=')
    experiment = readfile[1].split()[2]
    return int(aipsid), experiment

def removefile(dirpath):
    if os.path.exists(dirpath):
        os.remove(dirpath)

def get_best_sources(data):
    try: 
        best_src = get_file('%s.%s-ampcal.dat' % (data.name,data.klass))[0].split()[0]
        best_time= make_int(re.split(']',get_file('%s.%s-ampcal.dat' % (data.name,data.klass))[0].replace('[',']'))[-2].split(','))
    except IOError:
        #print 'Using BeSSeL'
        #let the bessel script output do the hard lifting here
        best_src = get_file('%s.%s-qual.dat' % (data.name,data.klass))[0].split()[0]
        best_time= make_int(re.split(']',get_file('%s.%s-qual.dat' % (data.name,data.klass))[0].replace('[',']'))[-2].split(','))
    return best_src,best_time

def split_block(total_block, key):
    #this arbitrarily splits up a large list into smaller lists via a keyword, 
    #in this case most likely the word 'Header'
    #this will allow the printed files to be separated into different times
    ind=[]
    for i in range(len(total_block)):
        cell=total_block[i]
        if cell.count(key)>0:
            ind+=[i]
    n_blocks=len(ind)
    spl_block=[]
    for i in range(n_blocks-1):
        tmp=[]
        tmp=total_block[ind[i]:ind[i+1]]
        spl_block+=[tmp]
    if n_blocks>1:
        tmp=[]
        tmp=total_block[ind[len(ind)-1]:len(total_block)]
        spl_block+=[tmp]
    else:
        spl_block=[total_block]
    return spl_block

def refine_fits(old_list, length_cat, symbol_cat, 
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

def match_source(data):
    source,blank=get_best_sources(data)
    flux_catalogue = {'1256-0547' : 15.0,
                      '3C279'     : 15.0,
                      '1921-293'  : 10.0,
                      '3C273'     : 20.0,
                      'G0634-2335': 1.0,
                      '1004-500'  : 0.4}
    s = [i for i in flux_catalogue]
    try:
        flux = flux_catalogue[difflib.get_close_matches(source, s)[0]]
    except (IndexError,TypeError) as e:
        flux = 1.0 #float(raw_input('Unknown flux for '+source+'. Please input (Jy)...\n'))
    return flux

##############################################################################
##############################################################################

if __name__=='__main__':
    main()

##############################################################################
##############################################################################
