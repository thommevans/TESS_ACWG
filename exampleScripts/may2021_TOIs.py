import os, sys, pdb
from TESS_ACWG import downloadTargetLists, processTargetLists
from TESS_ACWG import surveyGrids
import surveySetup

ADIR = os.getcwd()

########################################################################
# Step 1: Download the target list from NASA Exoplanet Archive.
# TODO: Have routines that download the target lists without
# having to do it manually.
toiFpath = downloadTargetLists.targetsUnpublishedTOIs() # dummy call for now
#TOI_FNAME = 'TOI_2021.05.03_04.50.08.csv'
#TOI_FNAME = 'TOI_2021.05.29_07.32.30.csv'
#TOI_FPATH = os.path.join( ADIR, TOI_FNAME )
#toiFpath = TOI_FPATH # for now use file that was downloaded manually

########################################################################
# Step 2: Process the csv file downloaded from NASA Exoplanet Archive.

toiPickle = processTargetLists.TOIs( csvIpath=toiFpath, pklOdir=ADIR )

    
########################################################################
# Step 3: Make the figures using the processed pickle file as input.
survey = { 'surveyName':'ACWG', 'framework':'ACWG', \
           'gridEdges':surveySetup.gridEdges, 'preCuts':surveySetup.preCutsTOIs, \
           'thresholdTSM':surveySetup.thresholdTSM, 'thresholdESM':surveySetup.thresholdESM }
RARanges = 'completeSet'
figFpaths = surveyGrids.TOIs( ipath=toiPickle, survey=survey, RARanges=RARanges, SMFlag = 'TSM' )
