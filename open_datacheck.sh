#!/bin/bash
# Open output plots from ParselTongue datacheck.py, found in
# /home/observer/correlations2/<exper>/datacheck/
# $1: <exper>, e.g. s001a

cd /home/observer/correlations2/$1/datacheck/
gv fringSN-geo.ps
gv fringSN-cont.ps
gv fringSN-line.ps
gv spectrumRAW-geo.ps
gv spectrumRAW-cont.ps
gv spectrumRAW-line.ps
gv spectrumCAL-geo.ps
gv spectrumCAL-cont.ps
gv spectrumCAL-line.ps
gv maserpeakVPLOT.ps
