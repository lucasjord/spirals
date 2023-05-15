#!/usr/local/bin/python3

#######################################################################

'''Script for importing an SN and SU files output from AIPS task TBOUT
containing real and imaginary solutions from FRING (no gain amplitudes)
then fitting a plane over said phases in time. Fit equation is conside-
red to be:

\\phi = \\phi_0 + A*\\Delta\\alpha*cos(\\delta) + B*\\Delta\\delta

- \\phi                           :: FRING phase
- \\phi_0                         :: Interpolated phase to (\\Delta
\\alpha*cos(\\delta),\\Delta\\delta)=(0,0)
- A                               :: Phase slope in RA direction
- \\Delta\\alpha*cos(\\delta)     :: Angular distance from phase cent-
re in RA direction
- B                               :: Phase slope in DEC direction
- \\Delta\\delta                  :: Angular distance from phase cent-
re in DEC direction

Which is a simple plane over the angular space.

Weighted Least-Squares is performed after removing 'bad' points to
calculate the \\phi_0, A and B parameters, which are then interpolated 
in time and angular extend back to the objects and returned as another 
SN table to be re-enterted into AIPS via tbin.

'''

#######################################################################

import matplotlib.pyplot as plt
import numpy as np
import matplotlib

import sys, os, random, copy, pdb
from math import exp, log, pi, atan2, sin, cos
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import inv

from itertools import product

############################################################
####################   Preamble   ##########################
############################################################

snfile = sys.argv[1]
number_of_stations = 4

#######################################################################
####################      MAIN Function      ##########################
#######################################################################

def main():
    head,data,end = read_multiband_delays2(snfile)
    inver = int(head[10].split()[2])
    refant = determine_ref_ant(data)
    d = clean_INDE_and_D(data)
    D = np.array(d)
    T = D[:,1]
    box = np.array([1,3,4,12,13,16])
    D2 = D[:,box]
    antenna_data1 = extract_telescope_phase(D2,number_of_stations)

    sufile = snfile.replace('SN2','SU1')
    su    = read_source_positions(sufile)
    names = read_source_names(sufile)
    #sid   = {v: k for k, v in names.iteritems()}
    if not len(su)==len(names):
        sys.exit('Source and Names files different lengths')
    number_of_orbits = len(su)-1
    for name in names:
        if 'G' in name:
            centre_source=name
            print(centre_source)
    # so we will be fitting each baseline independantly, and clustering groups of <number of orbits> sources together and fitting those

    # creating data structure to store parmameters
    stat_parms = {}
    for station_number in antenna_data1:
        prms = []
        for j in range(3):
            prms.append([])
        stat_parms.update({station_number:prms})

    P   = {}
    A   = {}
    B   = {}
    Chi = {}
    T   = {}
    for station_number in range(len(antenna_data1)):
        reverse_index = []
        for i in range(len(antenna_data1[station_number])):
            if not antenna_data1[station_number][i][1]==int(names[centre_source]):
                reverse_index+=[i]
        reverse_index = np.array(reverse_index)
        relevant_data = np.array(antenna_data1[station_number])
        antenna_data  = relevant_data[reverse_index]
        t  = []
        p  = []
        a  = []
        b  = []
        c  = []
        Z  = []

        #matplotlib.pyplot.close("all")

        n = range(len(antenna_data)-number_of_orbits)
        for N in n:
            #cutting out a bit of data to fit
            DATA = antenna_data[N:(number_of_orbits+N)]              
            tm   = DATA[:,0].mean()
            dt   = DATA[:,0].max() - DATA[:,0].min()
            #if range of data time larger than 20 mins, skip iteration
            if dt>=0.015: continue
            sids = DATA[:,1] 
            real = DATA[:,3]
            imag = DATA[:,4]
            snr  = DATA[:,5]
            x = []; y = []; z = []; w = []; s = []
            #nan capturing
            for row in range(len(DATA)):
                if 9999.9 in [real[row], imag[row]] or snr[row]==0.0:
                    continue
                else:
                    pos = list(np.array(su[str(int(sids[row]))])-np.array(su[names[centre_source]]))
                    x+= [pos[0]]; y+= [pos[1]]
                    z+= phase2([real[row]],[imag[row]])
                    w+= [1/((57/snr[row])**2+100.0)]
                    s+= [sids[row]]
            sid = np.array(s)

            perm = list(product([True,False],repeat=len(w)))
            final_perm = []
            for i in range(len(perm)):
                if perm[i].count(False)<=2:
                    final_perm+=[list(perm[i])]
            #best-fit capturing: itiratively remove one point at a time
            C = np.zeros(shape=(len(final_perm)))
            #print('{0:10s}{1:10s}{2:10s}{3:10s}{4:10s}'.format('         k','     CHI^2','     PHASE','       A_t','       B_t'))
            #print('         -----------------------------------------')
            best='blank'
            wrst='knalb'
            for k in range(len(final_perm)):
                index = np.array(final_perm[k])
                #print(index)
                X = 57*np.asarray(x)[index] ; Y = 57*np.asarray(y)[index]*cos(su[names[centre_source]][1])
                F = 57*np.matrix(np.array(z)[index]).T
                W = np.matrix(np.diag(np.array(w)[index]))
                M = np.matrix([X,X,Y]).T
                M[:,0]= M[:,0]*0.0 + 1.0
                phase_t = []
                A_t     = []
                B_t     = []
                #parms   = inv(M.T*W*M)*M.T*W*(F)
                try: 
                    parms   = inv(M.T*W*M)*M.T*W*(F)
                except np.linalg.LinAlgError: 
                    continue
                phase_t = parms.tolist()[0]  #deg
                A_t     = parms.tolist()[1]  #deg/deg
                B_t     = parms.tolist()[2]  #deg/deg
                time    = np.array(range(len(phase_t)))*0.0+tm

                phi     = np.asarray(z)[index]*57
                #rho     = correlations(inv(M.T*W*M))
                phi_m   = phase_t[0] + (A_t)*X + (B_t)*Y
                res     = (phi_m - phi)*(np.array(w)[index]**0.5)
                if len(X)-3<=0.0:
                    chis = np.nan
                else:
                    chis    = sum(res**2)/(len(X)-3)
                C[k]=chis
                #print('{0:10.0f}{1:10.2f}{2:10.2f}{3:10.2f}{4:10.2f}'.format(k,chis,phase_t[0],A_t[0],B_t[0]))
            best = np.array(final_perm[np.nanargmin(C)])
            wrst = np.array(final_perm[np.nanargmax(C)])
            #re-solving with best values
            index = ~(~best*wrst)
            X = 57*np.asarray(x)[index] ; Y = 57*np.asarray(y)[index]*cos(su[names[centre_source]][1])
            F = 57*np.matrix(np.array(z)[index]).T
            W = np.matrix(np.diag(np.array(w)[index]))
            M = np.matrix([X,X,Y]).T
            M[:,0]= M[:,0]*0.0 + 1.0
            phase_t = []
            A_t     = []
            B_t     = []
            #parms   = inv(M.T*W*M)*M.T*W*(F)
            try: parms   = inv(M.T*W*M)*M.T*W*(F)
            except np.linalg.LinAlgError: 
                t+=time.tolist()
                p+=[0]
                a+=[0]
                b+=[0]
                c+=[0]
                continue
            phase_t = parms.tolist()[0]  #deg
            A_t     = parms.tolist()[1]  #deg/deg
            B_t     = parms.tolist()[2]  #deg/deg
            time    = np.array(range(len(phase_t)))*0.0+tm
            phi     = np.asarray(z)[index]*57
            #rho     = correlations(inv(M.T*W*M))
            phi_m   = phase_t[0] + (A_t)*X + (B_t)*Y
            res     = (phi_m - phi)*(np.array(w)[index]**0.5)
            chis    = sum(res**2)/(len(X)-3)
            t+=time.tolist()
            p+=phase_t
            a+=A_t
            b+=B_t
            c+=[chis]
        
        T.update({station_number:np.array(t)})
        P.update({station_number:p})
        A.update({station_number:a})
        B.update({station_number:b})
        Chi.update({station_number:c})
        
    ############################################################

    outfile = snfile[:-11]+'3'+snfile[-10:]
    count=0
    content=[]
    for j in P:
        for i in range(len(T[j])):
            for k in su:
                x, y = list(np.array(su[k])-np.array(su[names[centre_source]]))
                X = 57*np.asarray(x) ; Y = 57*np.asarray(y)*cos(su[names[centre_source]][1])
                phase = P[j][i] + A[j][i]*X + B[j][i]*Y
                count=count+1
                content+=[('{:8.0f}{:>20.15f}D-01{:>15s}{:>11d}{:>11d}{:>11d}{:>11d}{:>11f}E+00{:>11.0f}{:>11f}E+00{:>11f}E+00{:>11f}E+00{:>11f}E+00{:>11f}E+00{:>11f}E+00{:>11f}E+00{:>11f}E+01{:>11s}').format(
                         count,T[j][i]*10.0,'5.324036E-04',int(k),j+1,1,-1,0,0,0,0,0,cos(phase*pi/180.0),sin(phase*pi/180.0),0,0,1.0,refant)]
    head[4] = f'NAXIS2  ={count:>21.0f} / Number of entries in table'
    head[10]= f'EXTVER  ={inver+1:>21.0f} / Version Number of table'

    os.system('rm '+outfile)
    outf = open(outfile,'a+')
    for line in head:    print(line,file=outf)
    outf.close()
    outf = open(outfile,'a+')
    for line in content: print(line,file=outf)
    outf.close()
    outf = open(outfile,'a+')
    print(end,file=outf)
    outf.close()

#######################################################################
####################   Defining Functions    ##########################
#######################################################################

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

def remove_repeats(list):
    track_count=[]
    for entry in list:
        if entry in track_count:
            continue
        else:
            track_count+=[entry]
    return track_count

def printall(thing):
    for i in thing:
        print(i)

##################################################
### AIPS slash script specific functions
##################################################

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
    tmp_list_4 = []
    for k in range(len(tmp_list_3)):
        if 'localhos' not in tmp_list_3[k]:
            tmp_list_4+=[tmp_list_3[k]]
    return tmp_list_4

def read_multiband_delays(file):
    Mtable = refine_fits(splitt(get_file(file))[20:],7,None,None,None)
    return Mtable

def read_multiband_delays2(file):
    SN     = list(get_file(file))
    indx   = int(SN.index('***BEGIN*PASS***')+1)
    header = SN[:indx] #keep header intact for later
    data   = splitt(SN[indx:-1])
    ender  = SN[-1]
    return header,data,ender

def clean_INDE_and_D(list_list):
    for i in range(len(list_list)):
        for j in range(len(list_list[i])):
            if list_list[i][j]=="'INDE'":
                list_list[i][j]=9999.9
            else:
                list_list[i][j]=float(list_list[i][j].replace('D','E'))
    return list_list

def mask_antenna(array_to_mask,n_antennas):
    antenna_key = {}
    for i in range(n_antennas):
        antenna_key.update({i:np.zeros(shape=(array_to_mask.shape[0],))})
    for j in range(array_to_mask.shape[0]):
        antenna_key[int(array_to_mask[j,2]-1)][j]=1.0
    return antenna_key       

def mask_source(array_to_mask,su_dict):
    source_key = {}
    for i in su_dict:
        source_key.update({int(i):np.zeros(shape=(array_to_mask.shape[0],))})
    for j in range(array_to_mask.shape[0]):
        source_key[int(array_to_mask[j,1])][j]=1.0
    return source_key

def read_source_positions(file):
    deg_rad = pi/180.0
    SUtable = refine_fits(splitt(get_file(file))[17:],5,None,None,None)
    ra = make_vector(SUtable,3)
    dec= make_vector(SUtable,4)
    RA = replaceDE(ra)
    DEC= replaceDE(dec)
    source_position = {}
    for i in range(len(SUtable)):
        source_position.update({SUtable[i][1]:(float(RA[i])*deg_rad,float(DEC[i])*deg_rad)})
    ###
    SUtable = refine_fits(splitt(get_file(file))[17:],6,None,None,None)
    ra = make_vector(SUtable,4)
    dec= make_vector(SUtable,5)
    RA = replaceDE(ra)
    DEC= replaceDE(dec)
    for i in range(len(SUtable)):
        source_position.update({SUtable[i][1]:(float(RA[i])*deg_rad,float(DEC[i])*deg_rad)})
    return source_position

def read_source_names(file):
    SUtable = refine_fits(splitt(get_file(file))[17:],5,None,None,None)+refine_fits(splitt(get_file(file))[17:],6,None,None,None)
    source_names = {}
    for i in range(len(SUtable)):
        source_names.update({SUtable[i][2]:SUtable[i][1]})
    return source_names

def extract_telescope_phase(table,n_antennas):
    ptable={}
    for j in range(n_antennas):
        ptable.update({j:[]})
    for i in range(len(table)):
        ptable[int(table[i][2]-1)]+=[table[i]]
    return ptable

def replaceDE(vector):
    vector_2 = []
    for i in range(len(vector)):
        vector_2+=[vector[i].replace('D','E')]
    return vector_2

def replaceINDEv(vector):
    vector_2 = []
    for i in range(len(vector)):
        if vector[i]=='INDE':
            vector_2+=[None]
        else:
            vector_2+=[vector[i]]
    return vector_2

def length(vec1,vec2):
    L = []
    for i in range(len(vec1)):
        if vec2[i]==None or vec1[i]==None:
            L+=[None]
        else:
            L+=[(float(vec2[i])**2+float(vec1[i])**2)**0.5]
    return L

def phase(vec1,vec2):
    p = []
    for i in range(len(vec1)):
        if vec2[i]==None or vec1[i]==None:
            p+=[None]
        elif vec2[i]=='INDE' or vec1[i]=='INDE':
            p+=[np.inf]
        else:
            p+=[atan2(float(vec2[i]),float(vec1[i]))]
    return p

def phase2(vec1,vec2):
    p = []
    for i in range(len(vec1)):
        if vec2[i]==None or vec1[i]==None:
            p+=[None]
        elif vec2[i]=='INDE' or vec1[i]=='INDE':
            p+=[np.inf]
        else:
            p+=[2*np.arctan(float(vec2[i])/(((float(vec1[i]))**2+(float(vec2[i]))**2)**0.5+float(vec1[i])))]
    return p

def replaceINDE(element):
    if element=='INDE':
        new_element=None
    else:
        new_element=element
    return new_element

def determine_ref_ant(mdel_table):
    sample = make_vector(mdel_table,-1)
    choice = random.choice(sample)
    return choice

def extract_station_source_data(repository,antenna,sourceID):
    section = repository[antenna]
    part    = [] 
    for i in range(len(section)):
        if section[i][3]==sourceID:
            part+=[section[i]]
    return part

def area(input1,input2,input3):
    '''Determines the area of a triangle formed by three points'''
    cx = 1/3*(input1[0]+input2[0]+input3[0]) # determine centroid x
    cy = 1/3*(input1[1]+input2[1]+input3[1]) # determine centroid y
    ang1 = np.arctan2(input1[0]-cx,input1[1]-cy)
    if ang1<0.0: ang1=ang1+2*pi
    ang2 = np.arctan2(input2[0]-cx,input2[1]-cy)
    if ang2<0.0: ang2=ang2+2*pi
    ang3 = np.arctan2(input3[0]-cx,input3[1]-cy)
    if ang3<0.0: ang3=ang3+2*pi
    refdict = {0:input1,1:input2,2:input3}
    s = np.array([ang1, ang2, ang3]).argsort()
    p0 = refdict[s[0]]
    p1 = refdict[s[1]]
    p2 = refdict[s[2]]
    area = 0.5*(-p1[1]*p2[0]+p0[1]*(-p1[0]+p2[1])+p0[0]*(p1[1]-p2[1])+p1[0]*p2[1])
    #return ang1, ang2, ang3
    return area
    
def baryocentric4zero(input1,input2,input3):
    '''Returns Baryocentric coordinates for zero
    if not 0 <= any <= 1 then outside or on triangle'''
    cx = 1/3*(input1[0]+input2[0]+input3[0]) # determine centroid x
    cy = 1/3*(input1[1]+input2[1]+input3[1]) # determine centroid y
    ang1 = np.arctan2(input1[0]-cx,input1[1]-cy)
    if ang1<0.0: ang1=ang1+2*pi
    ang2 = np.arctan2(input2[0]-cx,input2[1]-cy)
    if ang2<0.0: ang2=ang2+2*pi
    ang3 = np.arctan2(input3[0]-cx,input3[1]-cy)
    if ang3<0.0: ang3=ang3+2*pi
    refdict = {0:input1,1:input2,2:input3}
    s = np.array([ang1, ang2, ang3]).argsort()
    p0 = refdict[s[0]]
    p1 = refdict[s[1]]
    p2 = refdict[s[2]]
    A = 0.5*(-p1[1]*p2[0]+p0[1]*(-p1[0]+p2[1])+p0[0]*(p1[1]-p2[1])+p1[0]*p2[1])
    b = range(3)
    b[0] = 1/(2*A)*(p0[1]*p2[0] - p0[0]*p2[1])
    b[1] = 1/(2*A)*(p0[0]*p1[1] - p0[1]*p1[0])
    b[2] = 1-b[0]-b[1]
    return b[s[0]],b[s[1]],b[s[2]]

def calculate_baryocentre(xs,ys,pointxy):
    '''Need array of x and y, length3, and a point'''
    if not len(xs)==3 and len(ys)==3:
        sys.exit('Need three points for triangle')
    M = np.matrix([xs,ys,[1,1,1]])
    p = np.matrix([pointxy[0],pointxy[1],1]).T
    l = inv(M)*p
    return l

##################################################
### Least squares functions 
##################################################

def model_plane(p, A, B, ra_i, dec_i):
    phi_i = p + A*ra_i + B*dec_i
    return phi_i

def Pmatrix_plane(data):
    #creates partial derivative matrix
    #specific to model plane (const + Ax + By)
    x = data[0]
    y = data[1]
    z = data[2]
    n = len(x)
    p = range(n)
    for i in range(n):
        p[i] = [1, x[i], y[i]]
    P = np.matrix(p)
    return P

def model_curve(p, A, B, C, D, E, ra_i, dec_i):
    phi_i = p + A*ra_i + B*dec_i + C*ra_i*dec_i + D*ra_i**2 + E*dec_i**2
    return phi_i

def Pmatrix_curve(data):
    #creates partial derivative matrix
    #specific to model curve (const + Ax + By + Cxy + Dx^2 + Ey^2)
    x = data[0]
    y = data[1]
    z = data[2]
    n = len(x)
    p = range(n)
    for i in range(n):
        p[i] = [1, x[i], y[i], x[i]*y[i], x[i]**2, y[i]**2]
    P = np.matrix(p)
    return P

def Design_matrix(matrix):
    #matrix must be numpy matrix object
    from numpy.linalg import inv
    D = inv((matrix.transpose())*matrix)
    return D

def res1(y_model, y, error):
    #calculates residuals
    residules=(y-y_model)/float(error)
    return residules

def correlations(Design_matrix):
    d   = np.matrix.tolist(Design_matrix)
    rho = copy.deepcopy(d)
    for j in range(len(d)):
        for k in range(len(d)):
            rho[j][k]=d[j][k]/(d[j][j]*d[k][k])**0.5
    return rho

def printMatrixE(a):
   #print "Matrix["+("%d" %a.shape[0])+"]["+("%d" %a.shape[1])+"]"
   rows = a.shape[0]
   cols = a.shape[1]
   for i in range(0,rows):
      for j in range(0,cols):
         print("%6.3f" %a[i,j]),
      print
   print

#######################################################################
####################        Execution        ##########################
#######################################################################

if __name__=='__main__':
    main()

#######################################################################
