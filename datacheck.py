###############################################################################
# ParselTongue datacheck.py
'''
Check SPIRALS correlated data.
v0.1 2020-08-25
v0.2 2020-08-30
v0.3 2020-08-31
v0.4 2020-09-01
v0.5 2020-09-03
AIPS part is fine. Needs still:
- run from corr_scripts/.
- help file and comments in code (finished up to --load)
'''
###############################################################################

#poo

# Call in some python and parseltongue packages.
from AIPS import AIPS
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
import os, sys, subprocess, re, glob, argparse

# Help file and input arguments.
parser = argparse.ArgumentParser(add_help=True, description = '''
    A description of what the script does.
    ''')
parser.add_argument('exper', help='experiment code in lower case, e.g. s001a')
parser.add_argument('-l', '--load', help='fitld data into AIPS and list scans (fitld)', action='store_true')
parser.add_argument('-c', '--calibrate', help='fringe fit data and find maser peak (fring)', action='store_true')
parser.add_argument('-p', '--plot', help='make some plots of spectra (possm)', action='store_true')
parser.add_argument('-d', '--delete', help='delete data in AIPS catalog (zap)', action='store_true')
# Print help if no argument is given.
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

###############################################################################

# Get experiment code from positional argument.
experiment = args.exper

# Define correlation folder and data archive paths.
path_correlations = '/home/observer/correlations2/'
path_archive = '/mnt/parallax/spirals/'

# Set AIPS ID based on experiment code. s=28 (extended hex),source#,epoch#
if len(experiment) == 5:
    AIPS.userno=int(str(int(experiment[0],36))+experiment[2:4]+str(ord(experiment[4])-(ord('a')-1)))
else:
    print '\nExperiment code should be in s001a format. See --help. \n'
    sys.exit(0)

# Check if geo fits exists and define path.
if os.path.exists(path_correlations+experiment+'/geo/'+experiment+'-geo.fits'):
    path_geo = path_correlations+experiment+'/geo/'
elif os.path.exists(path_correlations+experiment+'/'+experiment+'-geo.fits'):
    path_geo = path_correlations+experiment+'/'
elif os.path.exists(path_archive+experiment+'/'+experiment+'-geo.fits'):
    path_geo = path_archive+experiment+'/'
else:
    print "\nCannot find, experiment"+"-geo.fits in correlation or archive (parallax) folders. See --help.\n"

# Check if cont fits exists and define path.
if os.path.exists(path_correlations+experiment+'/cont/'+experiment+'-cont.fits'):
    path_cont = path_correlations+experiment+'/cont/'
elif os.path.exists(path_correlations+experiment+'/'+experiment+'-cont.fits'):
    path_cont = path_correlations+experiment+'/'
elif os.path.exists(path_archive+experiment+'/'+experiment+'-cont.fits'):
    path_cont = path_archive+experiment+'/'
else:
    print "\nCannot find, experiment"+"-cont.fits in correlation or archive (parallax) folders. See --help.\n"

# Check if line fits exists and define path.
if os.path.exists(path_correlations+experiment+'/line/'+experiment+'-line.fits'):
    path_line = path_correlations+experiment+'/line/'
elif os.path.exists(path_correlations+experiment+'/'+experiment+'-line.fits'):
    path_line = path_correlations+experiment+'/'
elif os.path.exists(path_archive+experiment+'/'+experiment+'-line.fits'):
    path_line = path_archive+experiment+'/'
else:
    print "\nCannot find, experiment"+"-line.fits in correlation or archive (parallax) folders. See --help.\n"
    sys.exit(0)

# Define data tags for AIPS.
fits_geo = path_geo + experiment + '-geo.fits'
fits_cont = path_cont + experiment + '-cont.fits'
fits_line = path_line + experiment + '-line.fits'
name_geo = experiment.upper() + '.G'
name_cont = experiment.upper() + '.C'
name_line = experiment.upper() + '.L'
data_geo = AIPSUVData(name_geo,'UVDATA',1,1)
data_cont = AIPSUVData(name_cont,'UVDATA',1,1)
data_line = AIPSUVData(name_line,'UVDATA',1,1)

# Make output folder for datacheck files in experiment correlation folder.
if os.path.exists(path_correlations+experiment+'/datacheck/'):
    path_datacheck = path_correlations+experiment+'/datacheck/'
else:
    os.mkdir(path_correlations+experiment+'/datacheck/')
    path_datacheck = path_correlations+experiment+'/datacheck/'

# Remove existing auxiliary files (scans lists, source tables, spectra, AIPS plots).
for item in glob.glob(path_datacheck+'*.txt'):
    os.remove(item)
for item in glob.glob(path_datacheck+'*.ps'):
    os.remove(item)

###############################################################################

# (--load) Load data into AIPS and make some auxiliary files.
if args.load == True or (args.load == False and args.calibrate == False and args.plot == False and args.delete == False):

    # Delete esixting data in AIPS catalog.
    if data_geo.exists():
        data_geo.clrstat()
        data_geo.zap()
    if data_cont.exists():
        data_cont.clrstat()
        data_cont.zap()
    if data_line.exists():
        data_line.clrstat()
        data_line.zap()

    # Fitld geo, cont, line data.
    fitld = AIPSTask('FITLD')
    fitld.datain = fits_geo
    fitld.outname = name_geo
    fitld()
    fitld.datain = fits_cont
    fitld.outname = name_cont
    fitld()
    fitld.datain = fits_line
    fitld.outname = name_line
    fitld()

    # Make scan lists for geo, cont, line data.
    listr = AIPSTask('LISTR')
    listr.optype = 'SCAN'
    listr.indata = data_geo
    listr.outprint = path_datacheck + 'listr-geo.txt'
    listr()
    listr.indata = data_cont
    listr.outprint = path_datacheck + 'listr-cont.txt'
    listr()
    listr.indata = data_line
    listr.outprint = path_datacheck + 'listr-line.txt'
    listr()

    # Print SU source table of cont/line data.
    prtab = AIPSTask('PRTAB')
    prtab.indata = data_cont
    prtab.inext = 'SU'
    prtab.outprint = path_datacheck + 'SU.txt'
    prtab()

###############################################################################

# (--calibrate) Fring calibrators in geo, cont and line. Find and fring maser peak channel. Apply solutions.
if args.calibrate == True or (args.load == False and args.calibrate == False and args.plot == False and args.delete == False):

    # Find maser source name in SU table (G source).
    maser = subprocess.check_output(['awk', '/G[0-9]/ {print $3}', path_datacheck+'SU.txt']).split()
    
    # Find delay calibrator source names in SU table (F sources and initial fringe check source).
    fringe_finder = subprocess.check_output(['awk', '/F[0-9]/ {print $3}', path_datacheck+'SU.txt']).split()
    f1 = subprocess.Popen(['awk', '{print $2, $3}', path_datacheck+'SU.txt'], stdout=subprocess.PIPE)
    f2 = subprocess.Popen(['awk', '$1 == "1" {print $2}'], stdin=f1.stdout, stdout=subprocess.PIPE)
    f1.stdout.close()
    f3,err = f2.communicate()
    fringe_check = f3.split()
    calibrators = fringe_check + fringe_finder

    # Delete existing SN tables.
    data_geo.zap_table('SN',-1)
    data_cont.zap_table('SN',-1)
    data_line.zap_table('SN',-1)

    # Fringe fit all geo sources and delay calibrators (F sources + fringe check source) in cont/line.
    fring = AIPSTask('FRING')
    fring.aparm = AIPSList([3,0,0,0,0,2,3,0,0,0])
    fring.solint = -1
    fring.indata = data_geo
    fring()
    fring.indata = data_cont
    fring.calsour = AIPSList(calibrators)
    fring()
    fring.indata = data_line
    fring.calsour = AIPSList(calibrators)
    fring.dparm[8] = 5
    fring()

    # Delete existing CL tables, except CL 1 (default starting table).
    while data_geo.table_highver('AIPS CL') > 1:
        data_geo.zap_table('AIPS CL', 0)
    while data_cont.table_highver('AIPS CL') > 1:
        data_cont.zap_table('AIPS CL', 0)
    while data_line.table_highver('AIPS CL') > 1:
        data_line.zap_table('AIPS CL', 0)

    # Apply delay solutions (SN 1) to data. 
    clcal = AIPSTask('CLCAL')
    clcal.snver = 1
    clcal.gainver = 1
    clcal.gainuse = 2
    clcal.indata = data_geo
    clcal()
    clcal.indata = data_cont
    clcal()
    clcal.indata = data_line
    clcal()

    # Find and fring maser peak
    possm = AIPSTask('POSSM')
    possm.indata = data_line
    possm.solint = 0
    possm.stokes = 'I'
    possm.sources = AIPSList(maser)
    possm.nplots = 0
    possm.docal = -1
    possm.aparm[1] = -1
    possm.outtext = path_datacheck+'maser-spectrum.txt'
    possm()

    maser_peak = subprocess.check_output(f"sort -k6 -n {path_datacheck}maser-spectrum.txt | tail -n 1", shell=True).split()
    maser_channel = int(maser_peak[0])
    print("Maser peak channel: ", maser_channel)

    fring = AIPSTask('FRING')
    fring.indata = data_line
    fring.calsour = AIPSList(maser)
    fring.aparm = AIPSList([3,0,3,0,0,2,3,0,0,0])
    fring.dparm = AIPSList([3,-1,0])
    fring.solint = -1
    fring.docal = 1
    fring.gainuse = 2
    fring.bchan = maser_channel
    fring.echan = maser_channel
    fring()

    clcal = AIPSTask('CLCAL')
    clcal.indata = data_line
    clcal.snver = 2
    clcal.gainver = 2
    clcal.gainuse = 3
    clcal.sources = AIPSList(maser)
    clcal()


# Plotting some stuff
if args.plot == True or (args.load == False and args.calibrate == False and args.plot == False and args.delete == False):
    cowsplot

    data_geo.zap_table('PL',-1)
    data_cont.zap_table('PL',-1)
    data_line.zap_table('PL',-1)

# SN tables from fring
    snplt = AIPSTask('SNPLT')
    snplt.inext = 'SN'
    snplt.inver = 1
    snplt.opcode = 'ALSI'
    snplt.do3co = 1
    snplt.nplots = 4
    snplt.indata = data_geo
    snplt.optype = 'DELA'
    snplt.pixrange = AIPSList([-200e-9,200e-9])
    snplt()
    snplt.indata = data_cont
    snplt.optype = 'DELA'
    snplt.pixrange = AIPSList([-200e-9,200e-9])
    snplt()
    snplt.indata = data_line
    snplt.inver = 2
    snplt.optype = 'PHAS'
    snplt.pixrange = AIPSList([-180,180])
    snplt()

    lwpla = AIPSTask('LWPLA')
    lwpla.plver = 1
    lwpla.lpen = 1
    lwpla.dparm = AIPSList([0,1,0,0,0,4,41,8,0])
    lwpla.indata = data_geo
    lwpla.inver = data_geo.table_highver('PL')
    lwpla.outfile = path_datacheck+'/fringSN-geo.ps'
    lwpla()
    lwpla.indata = data_cont
    lwpla.inver = data_cont.table_highver('PL')
    lwpla.outfile = path_datacheck+'/fringSN-cont.ps'
    lwpla()
    lwpla.indata = data_line
    lwpla.inver = data_line.table_highver('PL')
    lwpla.outfile = path_datacheck+'/fringSN-line.ps'
    lwpla()

    data_geo.zap_table('PL',-1)
    data_cont.zap_table('PL',-1)
    data_line.zap_table('PL',-1)

# Raw spectra for each baseline
    possm = AIPSTask('POSSM')
    possm.docal = 0
    possm.solint = -1
    possm.nplots = 6
    possm.aparm = AIPSList([0,1,0,0,-180,180,0,0,3,0])
    possm.indata = data_geo
    possm()
    possm.indata = data_cont
    possm.sources = AIPSList(calibrators)
    possm()
    possm.indata = data_line
    possm.sources = AIPSList(maser)
    possm.codetype = 'AMP'
    possm.aparm[1] = -1
    possm.aparm[8] = 0
    possm.solint = 0
    possm()
    possm.nplots = 4
    possm.aparm[8] = 1
    possm.stokes = 'HALF'
    possm()

    lwpla = AIPSTask('LWPLA')
    lwpla.plver = 1
    lwpla.lpen = 1
    lwpla.dparm = AIPSList([0,1,0,0,0,4,41,8,0])
    lwpla.indata = data_geo
    lwpla.inver = data_geo.table_highver('PL')
    lwpla.outfile = path_datacheck+'/spectrumRAW-geo.ps'
    lwpla()
    lwpla.indata = data_cont
    lwpla.inver = data_cont.table_highver('PL')
    lwpla.outfile = path_datacheck+'/spectrumRAW-cont.ps'
    lwpla()
    lwpla.indata = data_line
    lwpla.inver = data_line.table_highver('PL')
    lwpla.outfile = path_datacheck+'/spectrumRAW-line.ps'
    lwpla()

    data_geo.zap_table('PL',-1)
    data_cont.zap_table('PL',-1)
    data_line.zap_table('PL',-1)

# Calibrated spectra and superarray
    possm = AIPSTask('POSSM')
    possm.docal = 1
    possm.solint = -1
    possm.nplots = 6
    possm.aparm = AIPSList([0,1,0,0,-180,180,0,0,3,0])
    possm.stokes = 'HALF'
    possm.indata = data_geo
    possm.gainuse = 2
    possm()
    possm.indata = data_cont
    possm.gainuse = 2
    possm.sources = AIPSList(calibrators)
    possm()
    possm.indata = data_line
    possm.sources = AIPSList(maser)
    possm.aparm[1] = 0
    possm.aparm[8] = 0
    possm.stokes = 'HALF'
    possm.solint = 0
    possm.gainuse = 3
    possm.bchan = maser_channel-100
    possm.echan = maser_channel+100
    possm()
    possm.stokes = 'I'
    possm.aparm[1] = 0
    possm.nplots = 0
    possm()
    possm.aparm[1] = -1
    possm()

    lwpla = AIPSTask('LWPLA')
    lwpla.plver = 1
    lwpla.lpen = 1
    lwpla.dparm = AIPSList([0,1,0,0,0,4,41,8,0])
    lwpla.indata = data_geo
    lwpla.inver = data_geo.table_highver('PL')
    lwpla.outfile = path_datacheck+'/spectrumCAL-geo.ps'
    lwpla()
    lwpla.indata = data_cont
    lwpla.inver = data_cont.table_highver('PL')
    lwpla.outfile = path_datacheck+'/spectrumCAL-cont.ps'
    lwpla()
    lwpla.indata = data_line
    lwpla.inver = data_line.table_highver('PL')
    lwpla.outfile = path_datacheck+'/spectrumCAL-line.ps'
    lwpla()

    data_geo.zap_table('PL',-1)
    data_cont.zap_table('PL',-1)
    data_line.zap_table('PL',-1)

# VPLOT on maser channel
    vplot = AIPSTask('VPLOT')
    vplot.indata = data_line
    vplot.sources = AIPSList(maser)
    vplot.bchan = maser_channel
    vplot.echan = maser_channel
    vplot.bparm = AIPSList([0,2,1,0,0,-180,180,-180,180,0])
    vplot.nplots = 6
    vplot.docal = 0
    vplot()
    vplot.docal = 1
    vplot.gainuse = 3
    vplot()

    lwpla = AIPSTask('LWPLA')
    lwpla.indata = data_line
    lwpla.plver = 1
    lwpla.lpen = 1
    lwpla.dparm = AIPSList([0,1,0,0,0,4,41,8,0])
    lwpla.inver = data_line.table_highver('PL')
    lwpla.outfile = path_datacheck+'/maserpeakVPLOT.ps'
    lwpla()

    data_line.zap_table('PL',-1)

    # Clean up AIPS catalog
if args.delete == True:

    if data_geo.exists():
        data_geo.clrstat()
        data_geo.zap()
    if data_cont.exists():
        data_cont.clrstat()
        data_cont.zap()
    if data_line.exists():
        data_line.clrstat()
        data_line.zap()

    print "\nData of", experiment, "deleted from AIPS user number:", AIPS.userno

print("\n######### \nFinished. See --help for details. \n######### \n")
