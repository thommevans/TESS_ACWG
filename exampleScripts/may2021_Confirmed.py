import os, sys, pdb
from TESS_ACWG import downloadTargetLists, processTargetLists
from TESS_ACWG import surveyGrids
from TESS_ACWG import massRadiusFigure
import surveySetup

ADIR = os.getcwd()

########################################################################
# Step 1: Download the target list from NASA Exoplanet Archive.
# TODO: Have routines that download the target lists without
# having to do it manually.
csvFpath = downloadTargetLists.targetsWithPublishedConfirmation() # dummy call for now
#CSV_FNAME = 'PS_2021.05.03_06.59.24.csv'
CSV_FNAME = 'PS_2021.05.29_07.25.27.csv'
CSV_FPATH = os.path.join( ADIR, CSV_FNAME )
csvFpath = CSV_FPATH # for now use file that was downloaded manually

########################################################################
# Step 2: Process the csv file downloaded from NASA Exoplanet Archive.
if 1:
    confirmedPickle = processTargetLists.Confirmed( csvIpath=csvFpath, pklOdir=ADIR )
else:
    confirmedPickle = 'confirmedProperties.pkl'
    
########################################################################
# Step 3: Make the figures using the processed pickle file as input.
#figFpaths = surveyGrids.Confirmed( ipath=confirmedPickle, surveyGrid='ACWG', \
#                                   obsSample='PublishedMassesOnly', thresholdTSM='ACWG' )
survey = { 'surveyName':'ACWG', 'obsSample':'PublishedMassesOnly', 'framework':'ACWG', \
           'gridEdges':surveySetup.gridEdges, 'thresholdTSM':surveySetup.thresholdTSM, \
           'preCuts':surveySetup.preCutsConfirmed }
figFpaths = surveyGrids.Confirmed( ipath=confirmedPickle, survey=survey )

########################################################################
# Step 4: Make a plot comparing the assumed mass-radius relation to
# the measured masses and radii of the confirmed planets.
massRadiusFpath = massRadiusFigure.Confirmed( ipath=confirmedPickle )

