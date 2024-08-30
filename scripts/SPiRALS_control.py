#!/home/observer/anaconda2/bin/ParselTongue

##############################################################################
# ParselTongue Script for Spirals calibration                          #
#                                                                            #
# Reads in data, runs indxr, checks antenna tables, writes station_file.inp, #
# control_file.inp and calibrator_file.inp (for continuum data), downloads   #
# TEC-maps and EOPs, and runs TECOR and EOP corrections, runs manual         #
# phase-cal, finge-fit on geo-sources, writes out files for fitting          #
#                                                                            #
# Version 3.1.0 (2021/10/30)                                                 #
#                                                                            #
# Changes:                                                                   #
#                                                                            #
# 1.0.1: Allow shifting of multiple sources                                  #
# 1.0.2: checkatmos added to look for zeros in zenith delay offsets          #
# 1.1.0: Slightly changed structure, full EVN capability                     #
# 2.0.0: Moved main script to definitions file                               #
# 2.0.1: Added DELZN, and imput model for fringe                             #
# 2.1.0: Added plot_tables, new options and experimental self-cal            #
# 2.2.0: Added AIPSVER, niter, get_key_flag, first epoch products            #
# 2.3.0: Added RDBE_check                                                    #
# 3.0.0: Added argparse, inverse multiview parms                             #
# 3.1.0: Added terminal arguments for predictable tasks                      #
#                                                                            #
##############################################################################

from AIPS import AIPS
import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument("epoch",
                    help="Experiment code e.g. s006i",
                    type=str)
parser.add_argument("aipsID",
                    help="AIPS ID, e.g 666",
                    type=int)
parser.add_argument("-i","--inter",
                    help="Interactive mode",
                    action="store_true",default=False)
parser.add_argument("-l","--log",
                    help="Logfile",
                    type=str,default='ParselTongue.log')
parser.add_argument("-L","--linedisk",
                    help="Disk to save/find line data",
                    type=int,default=1)
parser.add_argument("--load",
                    help="Load in data and make LISTR",
                    action='store_true',default=False)
parser.add_argument("--geo",
                    help="Disk to save/find line data",
                    action='store_true',default=False)
parser.add_argument("--prepr",
                    help="Calibrate line and cont data \
                    up to external amplitude cal",
                    action='store_true',default=False)
parser.add_argument("--mpcal",
                    help="Run manual phase cal",
                    action='store_true',default=False)
parser.add_argument("--cvel",
                    help="Run cvel",
                    action='store_true',default=False)
parser.add_argument("--mafring",
                    help="Fringe fit on maser",
                    action='store_true',default=False)
parser.add_argument("--splitimage",
                    help="Split and image calibrators \
                    and target fringe channel",
                    action='store_true',default=False)
parser.add_argument("--imulti",
                    help="Run inverse multiview and image",
                    action='store_true',default=False)

args = parser.parse_args()

def_file    = '/home/observer/spirals/definitions_spirals.py'
aipsver     = '31DEC19'

####################
# Input Parameters #
####################

logfile     = args.log 
AIPS.userno = args.aipsID
file_path   = '/home/observer/analysis/{}/files/'.format(args.epoch.lower())
eop_path    = '/home/observer/analysis/{}/eop/'.format(args.epoch.lower())
fit_path    = '/home/observer/analysis/fit/'
mv_path     = '/home/observer/analysis/multiview_4x/'

inter_flag  = int(args.inter) # interactive (1) or non-interactive (0) mode
n           = 3               # Number of UV-data files either to download, on disk
                              # or already in AIPS and to be analyzed
defdisk     = 1               # Default AIPS disk to use (can be changed later)

#############################################################################
###                  Do not change or move this part                     ####
[filename, outname, outclass]    = [range(n),range(n),range(n)]
[nfiles, ncount, doconcat]       = [range(n),range(n),range(n)]
[outdisk, flagfile, antabfile]   = [range(n),range(n),range(n)]
for i in range(n):
    [flagfile[i], antabfile[i], outdisk[i]] = ['','',defdisk]
    [nfiles[i], ncount[i], doconcat[i]]     = [0,1,-1]
#############################################################################

antname     = 'LBA'                   # Antenna order for FITLD    
refant      = 1                       # {1: 'CD', 2: 'HO', 3: 'KE', 4: 'YG'} if AUS, {1:'AT',2:'CD',3:'HO',4:'MP',5:'PA',6:'WA'} if LBA
calsource   = ['G232.62-0.99'][0]     # calibrator        '' => automatically
target      = ['J0725','J0735','J0748','F0730-1141'] # continuum sources '' => automatically
mp_source   = ['F0730-1141']          # fringe finder     '' => automatically
mp_timera   = [0, 6, 32, 2, 0, 6, 33, 1] #[0,0,0,0,0,0,0,0]  # constrain time range for fringe finder?
bandcal     = ['']                    # Bandpass calibrator

delzn_flag  =  0                  # Use DELZN (1) or fit_geblocks (0)?dd
dual_geo    =  0                  # Using dual frequency geodetic blocks?

geo_data_nr =  0                  # data file with geo data? (<0 for no geoblock)
cont        =  1                  # data file with continuum data?
line        =  2                  # data file with line data?
                                  # if you have only one dataset, use cont = line

cvelsource  = ['G232.62-0.99']    #  line sources '' => calsource
velchan     =  513                #  channel used for cvel, default 0=>Nchan/2+1
vel         = [23.08]             #  Velocity for cvel in band center

channel     = 513                 # channelqq used for fringe fit after cvel

pos_shift   = {'G232.62-0.99' : [ 0.016220, 0.245883],
               'J0725'        : [ 0.000452,-0.000098],
               'J0729'        : [ 0.000257,-0.000057],
               'J0735'        : [ 0.000508, 0.000271],
               'J0748'        : [-0.000178,-0.000474],
               'F0730-1141'   : [-0.000524, 0.000359]}
               					#arcsecond
                                # Which source should be shifted ''=calsource
                                # e.g. {'G211':[0.2,0.3], 'J0007':[0.1,0.3]}
                                # {'' : [0,0]} => no shift
#################
# Split Options #
#################

smooth       = [0,0,0]          # Smooth during split for line data
split_outcl  = 'SPLIT'          # outclass in SPLIT '' => 'SPLIT'

##################################
# Optional inputs for fringe fit #
##################################

[fr_n, fr_c, fr_d, fr_s] = ['','',1,1]
                                             # Input image to use in FRINGE:
                                             # ['NAME','CLASS',DISK,SEQUENCE]
smodel = [1,0]                 # SMODEL in FRING                                
solint = 0                     # SOLINT in FRING
nmaps  = 1                     # NMAPS in FRING

##############################################
# Optional imput for phased VLA observations #
##############################################

flux = {'' : [0,0,0,0]}         # Flux of cal source (for VLA)

######################
# Imaging Parameters #
######################

cellsize     = 0.00005            # CELLSIZE for IMAGR
imsize       = 512                # IMSIZE for IMAGR
bchan        = 501                # Start channel for data cube
echan        = 530                # End channel for data cube
niter        = 200                # niter in imagr
uvwtfn       = 'N'                # UVWTFN in IMAGR
robust       = 0                  # ROBUST in IMAGR
beam         = [0.001,0.001,0]    # Restoring beam: bmaj, bmin, PA
antennas     = [1,2,3,4]          # Also used for RPOSSM
ant_bls      = [1,4]              # only used for co_imagr_flag so far
imgr_timer   = [0,0,0,0,0,0,0,0]  # used for co_imagr and ma_imgr and imv_imagr

######################
#   SAD Parameters   #
######################

min_snr      = 20              # Fit only Peaks with SNR>min_snr
dyna         = 0.2             # Don't fit peaks less than 
                               # dyna*peak in channel
####################
# Self Calibration #
####################

phase_loop = [30./60,20./60.]      # Solints for Phase selfcal
amp_loop   = [1200.,360, 60.,30.]
dofit      = [[0], [0], [0],[0]]

#####################
# inverse Multiview #
#####################

mvwin = 15.0    # mins

###############
# Input Files #
###############

filename[0] = args.epoch+'-geo.fits'
outname[0]  = args.epoch.upper()+'_G'
outclass[0] = 'UVDATA'
outdisk[0]  = 1

filename[1] = args.epoch+'-cont.fits'
outname[1]  = args.epoch.upper()+'_C'
outclass[1] = 'UVDATA'
outdisk[1]  = 1

filename[2] = args.epoch+'-line.fits'
outname[2]  = args.epoch.upper()+'_L'
outclass[2] = 'UVDATA'
outdisk[2]  = args.linedisk

####################
# Data preparation #
####################
download_flag   = 0      # Download data from archive?
load_flag       = 0      # Load data from disk?
listr_flag      = 0      # Print out LISTR?
get_key_flag    = 0      # Download key-file from archive

#####################
# geoblock analysis #
#####################
RDBE_check      = 0      # Check Geoblock data for RDBE errors?
geo_prep_flag   = 0      # Run TECOR and EOP corrections?xs
geo_fringe_flag = 0      # Fringe fit the data?
doprt_flag      = 0      # Print out files for fitting?
dofit_flag      = 0      # Start geoblock fitting?
doplot_flag     = 0      # Plot fit_geoblock output?

#####################
# phaseref analysis #
#####################
tasav_flag      = 0      # Run tasav on original tables?
restore_su_flag = 0      # Restore original SU table from TASAV?
restore_fg_flag = 0      # Restore original FG table from TASAV?
#
pr_prep_flag    = 0      # Run TECOR, EOPs, ATMOS, PANG, and position shift?
apcal_flag      = 0      # Do amplitude calibration?  #cl6
pr_fringe_flag  = 0      # Do manual phase cal?
do_band_flag    = 0      # Do bandpass calibration?
cvel_flag       = 0      # Run CVEL on line data?
#
possm_flag      = 0      # Run POSSM to select channel?
ma_fringe_flag  = 0      # Fringe fit on maser and calibrate data?
co_fringe_flag  = 0      # Fringe fit on continuum and calibrate data?
#
snflg_flag      = 0      # Run snflg after fringe?
print_sn_flag   = 0      # Print out last SN table?
plot_tables     =-1      # Plot SN-Tables of which data set (-1 for no plots)
#
split_flag      = 0      # Split calibrated data?
fittp_flag      = 0      # Write calibrated data to disk?
co_imagr_flag   = 0      # Image continuum (target) sources?
#
ma_imagr_flag   = 0      # Image one channel of line (cvel-)sources?
cube_imagr_flag = 0      # Image data cube of line (cvel-)sources?
#
imv_prep_flag   = 1  	 # Run calib on calibrators and print out files for fitting?
imultiv_flag    = 0      # Run inverse MultiView on data? 
imv_app_flag    = 0      # apply imv solutions to target?
imv_imagr_flag  = 0      # image inverse multivew target (calsour) after imv?
#
rpossm_flag     = 0      # Produce first epoch spectrum?
ma_sad_flag     = 0      # Run SAD on cube?
plot_map        = 0      # Make spot map
#
phase_cal_flag  = 0      # Phase self-calibration on calibrator (cont)
amp_cal_flag    = 0      # Amplitude self-calibration on calibrator (cont)
apply_selfcal   = 0      # Apply self calibration to target sources?

####################
# Download options #
####################
# Only for archive download use arch_user = 'nopass' for e2earchive
if download_flag==1:
    file        = range(n)
    arch_user   = 'nopass'
    arch_pass   = 'nopass'
    file[0]=''
    file[1]=''
    file[2]=''

##############################################################################
# check if optional flags are set and overwrite parms if they are

if args.load       ==True:
    load_flag       = 1
    listr_flag      = 1
if args.geo        ==True:
    geo_prep_flag   = 1
    geo_fringe_flag = 1
    doprt_flag      = 1
    dofit_flag      = 1
    doplot_flag     = 1
if args.prepr      ==True:
    tasav_flag      = 1
    restore_su_flag = 1
    restore_fg_flag = 1
    pr_prep_flag    = 1
    apcal_flag      = 1
if args.mpcal      ==True:
    pr_fringe_flag  = 1
if args.cvel       ==True:
    cvel_flag       = 1
if args.mafring    ==True:
    ma_fringe_flag  = 1
if args.splitimage ==True:
    split_flag      = 1
    co_imagr_flag   = 1
    ma_imagr_flag   = 1
if args.imulti     ==True:
    imv_prep_flag   = 1
    imultiv_flag    = 1
    imv_app_flag    = 1
    imv_imagr_flag  = 1

##############################################################################
# Start main script
if __name__=='__main__':
    execfile(r''+def_file)

##############################################################################
