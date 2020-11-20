#!/usr/local/bin/parseltongue

#version 2.0 2019/01/29
#
#2.0 takes into account "blank" source input for all

import sys, os, pdb
import numpy as np
from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
import scipy.io as sio

###################################################
################### parameters ####################
###################################################

def main():
    dogeo  = 1
    docont = 1
    doline = 1

    AIPS.userno = 280+13+1
    EXPERIMENT  = "S001D"
    EXPCLASS    = "UVDATA"
    EXPSEQ      = 1
    EXPDISK     = 1
    file_path   = "/home/observer/analysis/s001d/files/"

    n = 2 #number of files
    #                     name         : spectral resolution fine=0.9766kHz, coarse=500kHz
    file_catalogue =   {'s001d-cont'  : 'coarse',
    	                's001d-line'  : 'fine'}

    data_names     =   {'s001d-cont'  : 's001d_lowres.fits',
                        's001d-line'  : 's001d_hires.fits'}
    #
    target_source       = [] #['G232.62-0.99']                  #masers/target             [] => ALL G*/J*
    target_calibrators  = [] #['J1725','J729','J1735','J1748']  #phasereference quasars    [] => ALL G*/J*
    correlation_key     = 'fine'   # correlation type to source from for spectra (ONLY RELEVANT FOR LINE DATA)
    geoblock_timer      = [0]      #[0, 8, 0, 0, 0, 20, 10, 0]

    file_names=[]
    for name in file_catalogue:
        file_names+=[name]

    if not len(file_names)==n:
        sys.exit('Incorrect number of files!')

    source_lists = {}
    for name in file_names:
        source_lists.update({name : generate_source_list(file_path, name)})
    #
    geoquasars=[]
    prefix=['F', 'G', 'J']

    if dogeo == 1:
        fringe_finders=[]
        print 'CREATING GEO FILE'
        for name in file_names:
            correlation_type=file_catalogue[name]
            if not correlation_type=='coarse':
                print name+' incorrect spectral type'
                continue
            else:
                sources = source_lists[name]
                for source in sources:
                    if prefix[0] in source:
                        fringe_finders+=[source]
                    if source[0] not in prefix:
                        geoquasars+=[source]
            total_geo = geoquasars + fringe_finders
            total_geo = remove_repeats(total_geo)
            for i in range(int(round(len(total_geo)/5.0))):
                if not 5*i+5==len(total_geo):
                    sublist=total_geo[5*i:5*i+5]
                    load(file_path+data_names[name], EXPERIMENT+'_G', 'UVDATA', 1, 1, sublist, 1, geoblock_timer)
        sort(EXPERIMENT+'_G', 'UVDATA', 1, 1)
        index(EXPERIMENT+'_G', 'UVDATA', 1, 1)
        print 'CREATED GEO FILE'
    else:
        print 'NOT CREATING GEO FILE'

    target_sources=[]
    target_calibratorss=[]
    all_cont_sources=[]
    if docont == 1:
        fringe_finders=[]
        print 'CREATING CONTINUUM FILE'
        for name in file_names:
            correlation_type=file_catalogue[name]
            if not correlation_type=='coarse':
                print name+' incorrect spectral type'
                continue
            else:
                sources = source_lists[name]
                if not target_source == []:
                    for target in target_source:
                        if not target in sources:
                            print target+' source not in FITS file '+name
                #
                if not target_calibrators == []:
                    for target in target_calibrators:
                        if not target in sources:
                            print target+' source not in FITS file '+name
                #
                for source in sources:
                    if prefix[0]==source[0]:
                        fringe_finders+=[source]
                    elif prefix[1]==source[0]:
                        target_sources+=[source]
                    elif prefix[2]==source[0]:
                        target_calibratorss+=[source]
            #
            if target_source==[]:
                target_source=target_sources
            if target_calibrators==[]:
                target_calibrators=target_calibratorss
            #
            total_sources=target_source+target_calibrators+fringe_finders
            total_sources=remove_repeats(total_sources)
            print total_sources
            total_sources=remove_repeats(total_sources)
            all_cont_sources+=total_sources
            if len(total_sources)>20:
                for i in range(int(round(len(total_sources)/5.0))):
                    if not 5*i+5==len(total_sources):
                        sublist=total_sources[20*i:5*i+20]
                        load(file_path+data_names[name], EXPERIMENT+'_C', 'UVDATA', 1, 1, sublist, 1, 0)
            else:
                load(file_path+data_names[name], EXPERIMENT+'_C', 'UVDATA', 1, 1, total_sources, 1, 0)
        sort(EXPERIMENT+'_C', 'UVDATA', 1, 1)
        index(EXPERIMENT+'_C', 'UVDATA', 1, 1)
        print 'CREATED CONT FILE'    
    else:
        print 'NOT CREATING CONTINUUM FILE'


    all_line_sources=[]
    if doline == 1:
        #pdb.set_trace()
        fringe_finders=[]
        print 'CREATING SPECTRAL LINE FILE'
        for name in file_names:
            error=0
            correlation_type=file_catalogue[name]
            if not correlation_type==correlation_key:
                    print name+' incorrect spectral type, requiring '+correlation_key
                    continue
            else:
                sources = source_lists[name]
                if not target_source == []:
                    for target in target_source:
                        if not target in sources:
                            print target+' source not in FITS file '+name
                            error=error+1
                #
                if not target_calibrators == []:
                    for target in target_calibrators:
                        if not target in sources:
                            print target+' source not in FITS file '+name
                            error=error+1
                #
                for source in sources:
                    if prefix[0]==source[0]:
                        fringe_finders+=[source]
                    elif prefix[1]==source[0]:
                        target_sources+=[source]
                    elif prefix[2]==source[0]:
                        target_calibratorss+=[source]
            # the next section deals with 'all' sources option
            if target_source==[]:
                target_source=target_sources
            if target_calibrators==[]:
                target_calibrators=target_calibratorss
            # 
            total_sources=target_source+target_calibrators+fringe_finders
            total_sources=remove_repeats(total_sources)
            print total_sources
            all_line_sources+=total_sources
            if len(total_sources)<=error:
                print 'Errors greater than number of required sources, not attempting load'
                continue
            else:
                if len(total_sources)>20:
                    for i in range(int(round(len(total_sources)/20.0))):
                        if not 20*i+20==len(total_sources):
                            sublist=total_sources[20*i:20*i+20]
                            print sublist
                            load(file_path+data_names[name], EXPERIMENT+'_L', 'UVDATA', 1, 1, sublist, 1, 0)
                else:
                    load(file_path+data_names[name], EXPERIMENT+'_L', 'UVDATA', 1, 1, total_sources, 1, 0)
        sort(EXPERIMENT+'_L', 'UVDATA', 1, 1)
        index(EXPERIMENT+'_L', 'UVDATA', 1, 1)
        print 'CREATED SPECTRAL LINE FILE'   
    else:
        print 'NOT CREATING SPECTRAL LINE FILE'

################################################### data manipulations

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

def make_vector(m, col):
    x=[]
    for i in m:
        x+=[i[col]]
    return x

def make_str(list):
    for i in range(len(list)):
        list[i]=str(list[i])
    return list

def make_float(list):
    for i in range(len(list)):
        try:
            list[i]=float(list[i])
        except ValueError:
            try:
                list[i]=list[i-1]
            except ValueError:
                list[i]=list[i+1]
    return list

def make_vector(m, col):
    x=[]
    for i in m:
        x+=[i[col]]
    return x

def make_matrix(vec1, vec2):
    Matrix=range(len(vec1))
    for i in range(len(vec1)):
        Matrix[i]=[vec1[i], vec2[i]]
    return Matrix

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

def refine_fits(old_list, length_cat, symbol_cat, column_search_number, column_retrival_number):
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
                    tmp_list_2+=[column]
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
    return tmp_list_3

def printall(list):
    #prints all items in a list (matlab style)
    for i in list:
        print i
    return

def remove_repeats(list):
    track_count=[]
    for entry in list:
        if entry in track_count:
            continue
        else:
            track_count+=[entry]
    return track_count

################################################ loading data from reference

def generate_source_list(path, filename):
    ###
    scanfiles = path+'pipeline/'+filename+'.SCAN'
    if not os.path.exists(scanfiles):
        sys.exit('Pipeline files not present, need to make a .SCAN file in LISTR and put them in '+scanfiles)
    l=splitt(get_file(scanfiles))
    l2 = refine_fits(l, 11, ":", None, 2)
    l3 = remove_repeats(l2)
    return l3

################################################ AIPS functions 

def load(path, outname, outclass, outseq, outdisk, sources, doconcat, timerange):
    fitld = AIPSTask('fitld')
    fitld.datain = path
    fitld.outname = outname
    fitld.outclass = outclass
    fitld.outseq = outseq
    fitld.outdisk = outdisk
    fitld.doconcat = doconcat #1 or -1
    fitld.sources = [None]+sources
    if timerange == 0:
        timerange = [0]
    else:
        ''
    fitld.timer = [None]+timerange
    fitld.go()

def sort(inname, inclass, inseq, indisk):
    uvsrt = AIPSTask('uvsrt')
    uvsrt.default
    uvsrt.inname = inname
    uvsrt.inclass = inclass
    uvsrt.inseq = inseq
    uvsrt.indisk = indisk
    uvsrt.outname = inname
    uvsrt.outclass = inclass
    uvsrt.outseq = inseq
    uvsrt.outdisk = indisk
    uvsrt.sort = 'TB'
    uvsrt.go()

def index(inname, inclass, inseq, indisk):
    indxr = AIPSTask('indxr')
    indxr.default
    indxr.inname = inname
    indxr.inclass = inclass
    indxr.inseq = inseq
    indxr.indisk = indisk
    indxr.go()

################################################

if __name__=='__main__':
    main()

