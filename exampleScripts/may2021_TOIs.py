import os, sys, pdb
from TESS_ACWG import downloadTargetLists, processTargetLists
from TESS_ACWG import surveyGrids
import surveySetup

ADIR = os.getcwd()

########################################################################
# Step 1: Download the target list from NASA Exoplanet Archive.

toiFpath = downloadTargetLists.targetsUnpublishedTOIs()
TOI_FNAME = toiFpath
TOI_FPATH = os.path.join( ADIR, TOI_FNAME )

########################################################################
# Step 2: Process the csv file downloaded from NASA Exoplanet Archive.
if 1:
    toiPickle = processTargetLists.TOIs( csvIpath=toiFpath, pklOdir=ADIR )
else:
    toiPickle = 'toiProperties.pkl'
    
########################################################################
# Step 3: Make the figures using the processed pickle file as input.
survey = { 'surveyName':'ACWG', 'framework':'ACWG', \
            'gridEdges':surveySetup.gridEdges, 'preCuts':surveySetup.preCutsTOIs, \
            'thresholdTSM':surveySetup.thresholdTSM }
RARanges = 'completeSet'
figFpaths = surveyGrids.TOIs( ipath=toiPickle, survey=survey, RARanges=RARanges )

