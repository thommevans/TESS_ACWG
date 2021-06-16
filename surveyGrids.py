import pdb, sys, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle
from . import Utils

FIGDIR = os.path.join( os.getcwd(), 'Figures' )
    

def quickCycle1():
    d = np.loadtxt( 'cycle1.csv', delimiter=',', dtype=str )
    pl = d[:,0]
    RpRE = np.array(d[:,2],dtype=float)
    TeqK = np.array(d[:,3],dtype=float)
    titleStr = 'JWST Cycle 1 targets that will be observed GTO + GO'
    titleStr += '\nIncludes all modes and observation types '
    titleStr += '(transits, eclipses, phase curves)'
    fig, ax, ax2 = generateBasicAxis( wideFormat=True, showLegend=False, \
                                      titleStr=titleStr )
    Tgrid, Rgrid = drawGrid( ax, surveyGrid='ACWG' )
    ax.plot( TeqK, RpRE, 'ok' )
    pdb.set_trace()


def Confirmed( ipath='confirmedProperties.pkl', survey={} ):
    """
    """
    wideFormat = True
    showNeptuneRadius = False
    showJupiterRadius = False
    addSignature = False
    figPaths = transmissionGridConfirmed( ipath=ipath, wideFormat=wideFormat, \
                                          showNeptuneRadius=showNeptuneRadius, \
                                          showJupiterRadius=showJupiterRadius, \
                                          survey=survey, addSignature=addSignature )
    for f in figPaths: # PDFs and PNGs
        for k in list( f.keys() ):
            fnew = f[k].replace( 'Confirmed_', 'Confirmed_{0}_'\
                                 .format( survey['obsSample'] ) )
            os.rename( f[k], fnew )
            plt.close( 'all' )
    return None


def TOIs( ipath='toiProperties.pkl', survey={}, months='all' ):
    """
    """
    wideFormat = True
    addSignature = False
    DecRestrictions = [ ['DecAll',None,None], ['DecNth',-20,None], ['DecSth',None,20] ]
    if months=='completeSet':
        months = [ 'RAall', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', \
                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' ]
    opaths = {}
    for i in DecRestrictions:
        opaths[i[0]] = {}
        for m in months:
            if m=='RAall':
                RA = [ None, None ]
            else:
                RA = Utils.getRARange( m )
            figPaths = transmissionGridTOIs( ipath=ipath, wideFormat=wideFormat, \
                                             addSignature=addSignature, survey=survey, \
                                             RAMin_hr=RA[0], RAMax_hr=RA[1], \
                                             DecMin_deg=i[1], DecMax_deg=i[2] )
            opaths[i[0]][m] = []
            for f in figPaths: # PDFs and PNGs
                for k in list( f.keys() ):
                    opath = f[k]
                    if f[k].find( '.pdf' )>0:
                        fnew = f[k].replace( '.pdf', '_{0}_{1}.pdf'.format( i[0], m ) )
                    elif opath.find( '.png' )>0:
                        fnew = f[k].replace( '.png', '_{0}_{1}.png'.format( i[0], m ) )
                    os.rename( f[k], fnew )
                    opaths[i[0]][m] += [ fnew ]
            plt.close( 'all' )
    return opaths


#############################################################################
# Transmission spectroscopy survey:


def transmissionGridTOIs( ipath='toiProperties.pkl', wideFormat=True, \
                          addSignature=True, survey={}, \
                          RAMin_hr=None, RAMax_hr=None, \
                          DecMin_deg=None, DecMax_deg=None ):
    """
    TOIs that have not been confirmed.
    """
    showGrid = True
    z, dateStr = readTOIProperties( ipath=ipath )
    ostr = 'TOIs'

    n0 = len( z['planetName'] )
    ixs0, cutStr, titleStr = survey['preCuts']( z )
    # Exclude targets outside the RA limits:
    RAStr, RAMin_hr, RAMax_hr = Utils.processRARestriction( RAMin_hr, RAMax_hr )
    ixsRA = ( z['RA_hr'][ixs0]>=RAMin_hr )*( z['RA_hr'][ixs0]<=RAMax_hr )
    # Exclude targets outside the Dec limits:
    DecStr, DecMin_deg, DecMax_deg = Utils.processDecRestriction( DecMin_deg, DecMax_deg )
    ixsDec = ( z['Dec_deg'][ixs0]>=DecMin_deg )*( z['Dec_deg'][ixs0]<=DecMax_deg )
    RADecStr = '{0}\n{1}'.format( RAStr, DecStr )
    ixs = np.arange( n0 )[ixs0][ixsRA*ixsDec]
    
    pl = z['planetName'][ixs]
    TSM = z['TSM'][ixs]
    Teq = z['TeqK'][ixs]
    Ts = z['TstarK'][ixs]
    MpVal = z['MpValME'][ixs]
    RpVal = z['RpValRE'][ixs]
    
    onames = {}

    # Radius-temperature grid plot listing the top-ranked planets in each cell:
    fig2, ax2 = plotTeqRpGrid( Teq, RpVal, Ts, TSM, pl, titleStr=titleStr, \
                               dateStr=dateStr, survey=survey, RADecStr=RADecStr )
    onames['2'] = '{0}_gridTopTSMs.pdf'.format( ostr )

    toiNote = 'TOIs with "PC" TFOPWG Disposition shown in darker font\n'
    toiNote += 'Masses estimated from empirical relation (adapted from Chen & Kipping 2017)'
    fig2.text( 0.08, 0.91-0.10, toiNote, \
               c='black', fontsize=14, horizontalalignment='left', \
               verticalalignment='bottom' )
    
    if addSignature==True:
        for ax in [ax2]:
            addSignatureToAxis( ax )
    
    figs = { '2':fig2 }
    print( '\nSaved:' )
    if wideFormat==True:
        odirExt = 'survey{0}/wideFormat'.format( survey['surveyName'] )
    else:
        odirExt = 'survey{0}/narrowFormat'.format( survey['surveyName'] )
    odir = os.path.join( FIGDIR, odirExt )
    if os.path.isdir( odir )==False:
        os.makedirs( odir )
    opathsPDF = {}
    opathsPNG = {}
    sourceStr = 'Source: NASA Exoplanet Archive ({0})'.format( dateStr )
    for k in ['2']:
        figs[k].text( 0.97, 0.01, sourceStr, fontsize=10, \
                      horizontalalignment='right', verticalalignment='bottom' )
        if addSignature==True:
            onames[k] = onames[k].replace( '.pdf', '_wSignature.pdf' )
        opathk = os.path.join( odir, onames[k] )
        figs[k].savefig( opathk )
        opathk_png = opathk.replace( '.pdf', '.png' )
        figs[k].savefig( opathk_png )
        opathsPDF[k] = opathk
        opathsPNG[k] = opathk_png
        print( '{0}\n{1}'.format( opathk, opathk_png ) )
        print( 'RADecStr = {0}'.format( RADecStr ) )
    return opathsPDF, opathsPNG
    
    
def transmissionGridTESS( publishedMasses=True, wideFormat=True, addSignature=True ):
    """
    Confirmed TESS planets without published mass.
    """
    showGrid = True
    z = readConfirmedTESSProperties( publishedMasses=publishedMasses )
    if publishedMasses==True:
        ostr = 'ConfirmedWithMassTESS'
        titleStr = 'Confirmed TESS planets with peer-reviewed published masses'
    else:
        ostr = 'ConfirmedNoMassTESS'
        titleStr = 'Confirmed TESS planets without peer-reviewed published masses'
    
    n0 = len( z['planetName'] )
    # Exclude very big and small stars:
    ixs0 = ( z['RsRS']>=0.05 )*( z['RsRS']<10 )
    # TODO = Add option to apply a bright limit?
    ixs = np.arange( n0 )[ixs0]#[ixs2][ixs3]
    print( '{0:.0f} planets have >5-sigma mass measurements.'.format( len( ixs ) ) )
    pl = z['planetName'][ixs]
    TSM = z['TSM'][ixs]
    Teq = z['TeqK'][ixs]
    Ts = z['TstarK'][ixs]
    MpVal = z['MpValME'][ixs]
    RpVal = z['RpValRE'][ixs]
    RpLsig = z['RpLsigRE'][ixs]
    RpUsig = z['RpUsigRE'][ixs]

    onames = {}

    # Radius-temperature grid plot listing the top-ranked planets in each cell:
    fig2, ax2 = plotTeqRpGrid( Teq, RpVal, Ts, TSM, pl, titleStr=titleStr )
    onames['2'] = '{0}_gridTopTSMs.pdf'.format( ostr )

    if publishedMasses==False:
        fig2.text( 0.10, 0.905-0.025, 'Masses estimated from empirical relation', \
                   c='black', fontsize=14, horizontalalignment='left', \
                   verticalalignment='bottom' )
    
    if addSignature==True:
        for ax in [ax2]:
            addSignatureToAxis( ax )
    
    figs = { '2':fig2 }
    print( '\nSaved:' )
    if wideFormat==True:
        odirExt = 'survey{0}/wideFormat'.format( surveyName )
    else:
        odirExt = 'survey{0}/narrowFormat'.format( surveyName )
    odir = os.path.join( FIGDIR, odirExt )
    if os.path.isdir( odir )==False:
        os.makedirs( odir )
    for k in ['2']:
        #if wideFormat==True:
        #    onames[k] = onames[k].replace( '.pdf', '_Wide.pdf' )
        if addSignature==True:
            onames[k] = onames[k].replace( '.pdf', '_wSignature.pdf' )
        opath = os.path.join( odir, onames[k] )
        figs[k].savefig( opath )
        print( opath )

    return None
    
    
def transmissionGridConfirmed( ipath='confirmedProperties.pkl', wideFormat=True, \
                               survey={}, addSignature=False, showGrid=True, \
                               showNeptuneRadius=False, showJupiterRadius=False ):
    """
    
    """
    z, dateStr = readConfirmedProperties( ipath=ipath )
    ostr = 'Confirmed'
    #titleStr = 'Confirmed planets (empirical mass-radius relation used where necessary)'

    # Not applying Dec restrictions to Confirmed planets for now:
    #DecStr, DecMin_deg, DecMax_deg = processDecRestriction( None, None )
    RADecStr = ''
    
    ixs, cutStr, titleStr = survey['preCuts']( z, survey['obsSample'] )
    
    print( '{0:.0f} planets have mass measurements or estimates'.format( len( ixs ) ) )
    print( 'and orbit stars with radii 0.05-10 R_Sun' )
    pl = z['planetName'][ixs]
    TSM = z['TSM'][ixs]
    Teq = z['TeqK'][ixs]
    Ts = z['TstarK'][ixs]
    MpVal = z['MpValME'][ixs]
    MpLsig = z['MpLsigME'][ixs]
    MpUsig = z['MpUsigME'][ixs]
    RpVal = z['RpValRE'][ixs]
    RpLsig = z['RpLsigRE'][ixs]
    RpUsig = z['RpUsigRE'][ixs]
    TESS = np.array( z['TESS'][ixs], dtype=int )

    onames = {}
    
    # Radius-temperature plot for all planets with well-measured mass:
    fig1a, ax1a = plotTeqRpScatter( pl, Teq, RpVal, Ts, TSM, TESS, applyTSMcuts=False, \
                                    wideFormat=wideFormat, survey=survey, \
                                    showGrid=showGrid, titleStr=titleStr, \
                                    indicateTESS=False, dateStr=dateStr, \
                                    showNeptuneRadius=showNeptuneRadius, \
                                    showJupiterRadius=showJupiterRadius )
    onames['1a'] = '{0}_allPlanets.pdf'.format( ostr )
    
    # Radius-temperature plot for all planets with well-measured mass
    # and TSM cuts applied:
    fig1b, ax1b = plotTeqRpScatter( pl, Teq, RpVal, Ts, TSM, TESS, applyTSMcuts=True, \
                                    wideFormat=wideFormat, survey=survey, \
                                    showGrid=showGrid, titleStr=titleStr, \
                                    indicateTESS=False, dateStr=dateStr, \
                                    showNeptuneRadius=showNeptuneRadius, \
                                    showJupiterRadius=showJupiterRadius )
    #TSMStr = Utils.getTSMStr( thresholdTSM )
    #fig1b.text( 0.5, 0.93, TSMStr, c='green', fontsize=14, \
    #            horizontalalignment='center', verticalalignment='bottom' )
    onames['1b'] = '{0}_TSMcutsApplied.pdf'.format( ostr )
                                    

    # Radius-temperature grid plot listing the top-ranked planets in each cell:
    fig2, ax2 = plotTeqRpGrid( Teq, RpVal, Ts, TSM, pl, titleStr=titleStr, \
                               dateStr=dateStr, survey=survey, RADecStr=RADecStr  )
    fig2.text( 0.10, 0.995, cutStr, c='black', fontsize=12, \
               horizontalalignment='left', verticalalignment='top' )
    onames['2'] = '{0}_gridTopTSMs.pdf'.format( ostr )

    # Scatter plots without the grid:
    fig3a, ax3a = plotTeqRpScatter( pl, Teq, RpVal, Ts, TSM, TESS, applyTSMcuts=False, \
                                    wideFormat=wideFormat, survey=survey, \
                                    showGrid=showGrid, titleStr=titleStr, \
                                    indicateTESS=True, dateStr=dateStr, \
                                    showNeptuneRadius=showNeptuneRadius, \
                                    showJupiterRadius=showJupiterRadius )
    onames['3a'] = '{0}_allPlanets_showsTESS.pdf'.format( ostr )
    fig3b, ax3b = plotTeqRpScatter( pl, Teq, RpVal, Ts, TSM, TESS, applyTSMcuts=True, \
                                    wideFormat=wideFormat, survey=survey, \
                                    showGrid=showGrid, titleStr=titleStr,
                                    indicateTESS=True, dateStr=dateStr, \
                                    showNeptuneRadius=showNeptuneRadius, \
                                    showJupiterRadius=showJupiterRadius )
    #TSMStr = Utils.getTSMStr( thresholdTSM )
    #fig3b.text( 0.5, 0.935, TSMStr, c='green', fontsize=14, \
    #            horizontalalignment='center', verticalalignment='bottom' )
    onames['3b'] = '{0}_TSMcutsApplied_showsTESS.pdf'.format( ostr )

    if addSignature==True:
        for ax in [ax1a,ax1b,ax2,ax3a,ax3b]:
            addSignatureToAxis( ax )
    
    figs = { '1a':fig1a, '1b':fig1b, '2':fig2, '3a':fig3a, '3b':fig3b }
    print( '\nSaved:' )
    if wideFormat==True:
        odirExt = 'survey{0}/wideFormat'.format( survey['surveyName'] )
    else:
        odirExt = 'survey{0}/narrowFormat'.format( survey['surveyName'] )
    odir = os.path.join( FIGDIR, odirExt )
    if os.path.isdir( odir )==False:
        os.makedirs( odir )
    opathsPDF = {}
    opathsPNG = {}
    sourceStr = 'Source: NASA Exoplanet Archive ({0})'.format( dateStr )
    for k in ['1a','1b','2','3a','3b']:
        figs[k].text( 0.97, 0.01, sourceStr, fontsize=10, \
                      horizontalalignment='right', verticalalignment='bottom' )
        if addSignature==True:
            onames[k] = onames[k].replace( '.pdf', '_wSignature.pdf' )
        opathk = os.path.join( odir, onames[k] )
        figs[k].savefig( opathk )
        opathk_png = opathk.replace( '.pdf', '.png' )
        figs[k].savefig( opathk_png )
        opathsPDF[k] = opathk
        opathsPNG[k] = opathk_png
        print( '{0}\n{1}'.format( opathk, opathk_png ) )
    return opathsPDF, opathsPNG


def transmissionPredictedTESS( showSolarSystem=True, wideFormat=False, \
                               surveyModule='ACWG', showGrid=True, \
                               showStellarTrack=True, showNeptuneRadius=True, \
                               showJupiterRadius=True, crossHair=None, \
                               addSignature=False ):
    """
    Focusing on TESS predicted yield, not worrying about survey grid.
    This routine is currently out of date.
    """
    z = readConfirmedProperties()
    ostr = 'predictedTESS'
    
    n0 = len( z['planetName'] )
    ixs1 = np.isfinite( z['MpValME'] )*np.isfinite( z['MpLsigME'] )
    ixs2 = ( z['MpValME'][ixs1]/z['MpLsigME'][ixs1]>=5 )
    # Exclude very big and small stars:
    ixs3 = ( z['RsRS'][ixs1][ixs2]>=0.05 )*( z['RsRS'][ixs1][ixs2]<10 )
    ixs = np.arange( n0 )[ixs1][ixs2][ixs3]
    print( '{0:.0f} planets have >5-sigma mass measurements.'.format( len( ixs ) ) )
    pl = z['planetName'][ixs]
    TSM = z['TSM'][ixs]
    Teq = z['TeqK'][ixs]
    Ts = z['TstarK'][ixs]
    MpVal = z['MpValME'][ixs]
    MpLsig = z['MpLsigME'][ixs]
    MpUsig = z['MpUsigME'][ixs]
    RpVal = z['RpValRE'][ixs]
    RpLsig = z['RpLsigRE'][ixs]
    RpUsig = z['RpUsigRE'][ixs]
    TESS = np.array( z['TESS'][ixs], dtype=int )
    
    # Since we're illustrating the *predicted* contribution by TESS,
    # remove the planets that have been detected by TESS already:
    ixs = ( TESS==0 )

    
    # Radius-temperature plot for all planets with well-measured mass:
    titleStr = 'Top-ranked transmission spectroscopy targets '
    if wideFormat==True:
        titleStr += 'with published ($>5\\sigma$) masses'
    else:
        titleStr += '\nwith published ($>5\\sigma$) masses'
        
    fig0a, ax0a = plotTeqRpScatter( pl, Teq[ixs], RpVal[ixs], Ts[ixs], TSM[ixs], \
                                    TESS[ixs], applyTSMcuts=False, ms=6, alpha=1, \
                                    starColors=True, showSolarSystem=showSolarSystem, \
                                    showStellarTrack=showStellarTrack, \
                                    wideFormat=wideFormat, titleStr=titleStr, \
                                    surveyModule=surveyModule, showGrid=showGrid, \
                                    indicateTESS=False, dateStr=dateStr, \
                                    showNeptuneRadius=showNeptuneRadius, \
                                    showJupiterRadius=showJupiterRadius )
    if crossHair is not None:
        ix = ( pl==crossHair )
        ax0a.axvline( Teq[ix], c='HotPink' )
        ax0a.axhline( RpVal[ix], c='HotPink' )
        pdb.set_trace()
    if wideFormat==True:
        odirExt = 'survey{0}/wideFormat'.format( surveyModule )
    else:
        odirExt = 'survey{0}/narrowFormat'.format( surveyModule )
    odir = os.path.join( FIGDIR, odirExt )
    if os.path.isdir( odir )==False:
        os.makedirs( odir )
    onames = {}
    for k in ['0a','0b','1a','1b','2a','2b','3']:
        onames[k] = '{0}_{1}.pdf'.format( ostr, k )
    if showSolarSystem==True:
        for k in ['0a','0b','1a','1b','2a','2b','3']:
            onames[k] = onames[k].replace( '.pdf', '_wSS.pdf' )
    if showStellarTrack==True:
        for k in ['0a','0b','1a','1b','2a','2b','3']:
            onames[k] = onames[k].replace( '.pdf', '_wST.pdf' )
    opaths = {}
    for k in ['0a','0b','1a','1b','2a','2b','3']:
        opaths[k] = os.path.join( odir, onames[k] )
        #if wideFormat==True:
        #    opaths[k] = opaths[k].replace( '_RpTeq_', '_RpTeq_Wide_' )
    fig0a.savefig( opaths['0a'] )
    
    # Radius-temperature plot for all planets with well-measured mass:
    if wideFormat==True:
        titleStr = 'Planets with published masses plus predicted TESS planets '
    else:
        titleStr = 'Planets with published masses\nplus predicted TESS planets '
    fig0b, ax0b = plotTeqRpScatter( pl[ixs], Teq[ixs], RpVal[ixs], Ts[ixs], TSM[ixs], \
                                    TESS[ixs], applyTSMcuts=False, ms=6, alpha=1, \
                                    starColors=True, showSolarSystem=showSolarSystem, \
                                    showStellarTrack=showStellarTrack, \
                                    wideFormat=wideFormat, titleStr=titleStr, \
                                    surveyModule=surveyModule, showGrid=showGrid, \
                                    indicateTESS=False, dateStr=dateStr, \
                                    showNeptuneRadius=showNeptuneRadius, \
                                    showJupiterRadius=showJupiterRadius )
    plotTeqRpTESS( ax0b, thresholdTSM=thresholdTSM )
    fig0b.savefig( opaths['0b'] )

    # Same as previous but with low alpha value for known planets:
    titleStr = 'Planets with published masses (faint)' 
    if wideFormat==True:
        titleStr += ' plus predicted TESS planets (solid)'
    else:
        titleStr += '\nplus predicted TESS planets (solid)'
    fig1, ax1 = plotTeqRpScatter( pl[ixs], Teq[ixs], RpVal[ixs], Ts[ixs], TSM[ixs], \
                                  TESS[ixs], applyTSMcuts=False, ms=6, alpha=0.3, \
                                  starColors=True, showSolarSystem=showSolarSystem, \
                                  showStellarTrack=showStellarTrack, \
                                  wideFormat=wideFormat, titleStr=titleStr, \
                                  surveyModule=surveyModule, showGrid=showGrid, \
                                  indicateTESS=False, dateStr=dateStr, \
                                  showNeptuneRadius=showNeptuneRadius, \
                                  showJupiterRadius=showJupiterRadius )
    #opath1a = os.path.join( FIGDIR, onames['1a'] )
    fig1.savefig( opaths['1a'] )

    # Add the bright predicted TESS planets:
    titleStr = 'Planets with published masses (grey)' 
    if wideFormat==True:
        titleStr += ' plus predicted TESS planets'
    else:
        titleStr += '\nplus predicted TESS planets'
    #titleStr = 'Top-ranked transmission spectroscopy targets '
    #titleStr += 'with published ($>5\\sigma$) masses + TESS predicted'
    plotTeqRpTESS( ax1, showSolarSystem=False, showNeptuneRadius=False, \
                   showJupiterRadius=False )
    #opath1b = os.path.join( FIGDIR, onames['1b'] )
    fig1.savefig( opaths['1b'] )
    
    # Radius-temperature plot for all planets with well-measured mass:
    fig2, ax2 = plotTeqRpScatter( pl[ixs], Teq[ixs], RpVal[ixs], Ts[ixs], TSM[ixs], \
                                  TESS[ixs], applyTSMcuts=False, ms=3, alpha=1, \
                                  starColors=False, showSolarSystem=showSolarSystem, \
                                  showStellarTrack=showStellarTrack, \
                                  wideFormat=wideFormat, titleStr=titleStr, \
                                  surveyModule=surveyModule, showGrid=showGrid, \
                                  indicateTESS=False, dateStr=dateStr, \
                                  showNeptuneRadius=showNeptuneRadius, \
                                  showJupiterRadius=showJupiterRadius )
    if wideFormat==True:
        odirExt = 'survey{0}/wideFormat'.format( surveyModule )
    else:
        odirExt = 'survey{0}/narrowFormat'.format( surveyModule )
    odir = os.path.join( FIGDIR, odirExt )
    if os.path.isdir( odir )==False:
        os.makedirs( odir )
    #opath2a = os.path.join( FIGDIR, onames['2a'] )
    fig2.savefig( opaths['2a'] )

    # Add the bright predicted TESS planets:
    plotTeqRpTESS( ax2, showSolarSystem=False, showNeptuneRadius=False, \
                   showJupiterRadius=False )
    #opath2b = os.path.join( FIGDIR, onames['2b'] )
    fig2.savefig( opaths['2b'] )

    # Make a plot with the TESS planets only:
    titleStr = 'TESS predicted planets (Barclay et al., 2018)'
    fig3, ax3, ax3Legend = generateBasicAxis( wideFormat=wideFormat, titleStr=titleStr )
    plotTeqRpTESS( ax3, titleStr=titleStr, showSolarSystem=showSolarSystem, \
                   showNeptuneRadius=showNeptuneRadius, \
                   showJupiterRadius=showJupiterRadius )
    if showStellarTrack==True:
        ax3 = addStellarTrack( ax3 )
    #opath3 = os.path.join( FIGDIR, onames['3'] )
    fig3.savefig( opaths['3'] )

    print( '\nSaved:' )
    for k in ['0a','0b','1a','1b','2a','2b','3']:
        print( onames[k] )
    #print( '\nSaved:\n{0}\n{1}\n{2}\n{3}\n{4}\n{5}\n{6}\n'\
    #       .format( opath0a, opath0b, opath1a, opath1b, opath2a, opath2b, opath3 ) )
    return None


#############################################################################
# Emission spectroscopy survey:


#############################################################################
# Utility routines:

def printTopPredictedSubNeptunes( z, onlySubNeptunes=True ):
    ixs = ( z['cad2min']==1 )*( z['RpValRE']>1.5 )*( z['RpValRE']<4 )
    brightLim = 6
    ixs *= ( z['Vmag']>brightLim )*( z['Jmag']>brightLim )*( z['Kmag']>brightLim )
    TSM = z['TSM'][ixs]
    RpRE = z['RpValRE'][ixs]
    MpME = z['MpValME'][ixs]
    aRs = z['aRs'][ixs]
    ideg = z['ideg'][ixs]
    T14hr = z['T14hr'][ixs]
    Pday = z['Pday'][ixs]
    Teq = z['TeqK'][ixs]
    RsRS = z['RsRS'][ixs]
    Ts = z['TstarK'][ixs]
    V = z['Vmag'][ixs]
    J = z['Jmag'][ixs]
    K = z['Kmag'][ixs]
    Tedges = np.arange( 200, 1100, 100 )
    nbins = len( Tedges )-1
    ntop = 5
    if onlySubNeptunes==True:
        print( '\nTop predicted **sub-Neptunes** cooler than 1000K in 100K bins:' )
    else:
        print( '\nTop predicted planets cooler than 1000K in 100K bins:' )
    for i in range( nbins ):
        Tl = Tedges[i]
        Tu = Tedges[i+1]
        ixs1 = ( Teq>=Tl )*( Teq<Tu )
        print( '\n>>>> Between {0:.0f}-{1:.0f}K:'.format( Tl, Tu ) )
        hstr = '         TSM  Mp(ME)  Rp(RE)  Rs(RS)   Ts(K)   Tp(K)     '
        hstr += 'P(d)     aRs   i(deg)  T14(h)      V       J       K'
        print( hstr )
        print( '{0}'.format( 110*'-' ) )
        m = int( ixs1.sum() )
        ixs2 = np.arange( m )[np.argsort( TSM[ixs1] )][::-1]
        n = int( min( [ len( ixs2 ), ntop ] ) )
        for j in range( ntop ):
            ostr = '{0:.0f}. '.format( j+1 ).rjust( 5 )
            ostr += '{0:.1f} '.format( TSM[ixs1][ixs2][j] ).rjust( 8 )
            ostr += '{0:.2f} '.format( MpME[ixs1][ixs2][j] ).rjust( 8 )
            ostr += '{0:.2f} '.format( RpRE[ixs1][ixs2][j] ).rjust( 8 )
            ostr += '{0:.2f} '.format( RsRS[ixs1][ixs2][j] ).rjust( 8 )
            ostr += '{0:.0f} '.format( Ts[ixs1][ixs2][j] ).rjust( 8 )
            ostr += '{0:.0f} '.format( Teq[ixs1][ixs2][j] ).rjust( 8 )
            ostr += '{0:.2f} '.format( Pday[ixs1][ixs2][j] ).rjust( 9 )
            ostr += '{0:.2f} '.format( aRs[ixs1][ixs2][j] ).rjust( 8 )
            ostr += '{0:.1f} '.format( ideg[ixs1][ixs2][j] ).rjust( 8 )
            ostr += '{0:.1f} '.format( T14hr[ixs1][ixs2][j] ).rjust( 8 )
            ostr += '{0:.1f} '.format( V[ixs1][ixs2][j] ).rjust( 8 )
            ostr += '{0:.1f} '.format( J[ixs1][ixs2][j] ).rjust( 8 )
            ostr += '{0:.1f} '.format( K[ixs1][ixs2][j] ).rjust( 8 )
            print( ostr )

    return None


def plotTeqRpTESS( ax, showSolarSystem=True, showNeptuneRadius=True, dateStr='', \
                   showJupiterRadius=True, titleStr='' ):
    t = readPredictedProperties()
    m = 10
    ixs = ( t['cad2min']==1 )*( ( t['Vmag']<m )+( t['Jmag']<m )+( t['Kmag']<m ) )
    TSM = t['TSM'][ixs]
    Teq = t['TeqK'][ixs]
    RpVal = t['RpValRE'][ixs]
    printTopPredictedSubNeptunes( t )
    Ts = t['TstarK'][ixs]
    n = len( Teq )
    alpha = 1
    ms = 6
    z0 = 1000
    applyTSMcuts = False
    for i in range( n ):
        c = Utils.getStarColor( Ts[i] )
        if thresholdTSM=='ACWG':
            thresholdTSM = surveyModule.getThresholdTSM( RpVal[i], framework=thresholdTSM )
            pdb.set_trace() # TEMPORARY
        if applyTSMcuts==False: # plotting everything regardless of TSM
            ax.plot( [Teq[i]], [RpVal[i]], 'o', ms=ms, alpha=alpha, \
                     mfc=c, mec=c, zorder=z0+i )
        elif TSM[i]>thresholdTSM: # if TSM cuts applied, this one is high enough
            ax.plot( [Teq[i]], [RpVal[i]], 'o', ms=ms, alpha=alpha, \
                     mfc=c, mec=c, zorder=z0+i )
        else: # otherwise plot as a smaller background point
            ax.plot( [Teq[i]], [RpVal[i]], 'o', ms=0.5*ms, alpha=alpha, \
                     mfc=c0, mec=c0, zorder=0 )
    ax = addSolarSystem( ax, showSolarSystem=showSolarSystem, \
                         showNeptuneRadius=showNeptuneRadius, \
                         showJupiterRadius=showJupiterRadius )
    return ax
    
    
def drawGrid( ax, cgrid=None, zorder=0, survey={} ):
    if cgrid is None:
        cgrid = np.array( [ 201, 148, 199 ] )/256.
    
    Tgrid, Rgrid = survey['gridEdges']( survey['surveyName'] )

    # Number of cells along each axis:
    nT = len( Tgrid )
    nR = len( Rgrid )
    
    for i in range( nT ):
        ax.plot( [Tgrid[i],Tgrid[i]], [Rgrid.min(),Rgrid.max()], '-', \
                 c=cgrid, zorder=zorder )
    for i in range( nR ):
        ax.plot( [Tgrid.min(),Tgrid.max()], [Rgrid[i],Rgrid[i]], '-', \
                 c=cgrid, zorder=zorder )
    return Tgrid, Rgrid


def addSignatureToAxis( ax ):
    c = 0.2*np.ones( 3 )
    ax.text( 0.02, 1, 'Figure by T. Mikal-Evans', fontsize=12, color=c, \
             verticalalignment='top', transform=ax.transAxes, zorder=0 )
    return None


def plotTeqRpGrid( TeqK, RpRE, TstarK, TSM, pl, cgrid=None, titleStr='', \
                   RADecStr='', dateStr='', wideFormat=True, survey={} ):
    if cgrid is None:
        cgrid = np.array( [ 201, 148, 199 ] )/256.
    fig, ax, ax2 = generateAxisGrid( wideFormat=wideFormat, titleStr=titleStr, \
                                     RADecStr=RADecStr )
    
    Tgrid, Rgrid = survey['gridEdges']( survey['surveyName'] )
    nT = len( Tgrid )
    nR = len( Rgrid )
    xLines = np.arange( 0.5, nT+0.5 )
    yLines = np.arange( 0.5, nR+0.5 )
    for i in range( nT ):
        ax.plot( [xLines[i],xLines[i]], [yLines.min(),yLines.max()], '-', \
                 c=cgrid, zorder=1 )
    for i in range( nR ):
        ax.plot( [xLines.min(),xLines.max()], [yLines[i],yLines[i]], '-', \
                 c=cgrid, zorder=1 )
    ax, TSMstr = addTopTSMs( ax, pl, TSM, TeqK, RpRE, TstarK, Tgrid, Rgrid, \
                             xLines, yLines, survey=survey )
    formatAxisTicks( ax )
    ax.xaxis.set_ticks( xLines, minor=False )
    ax.yaxis.set_ticks( yLines, minor=False )
    ax.set_xticklabels( Tgrid )
    ax.set_yticklabels( Rgrid )
    if wideFormat==False:
        subtitleY = 0.94
        dySubTitle = 0.01
    else:
        subtitleY = 0.925
        dySubTitle = 0.015
    fig.text( 0.08, subtitleY, TSMstr, c='green', fontsize=14, \
              horizontalalignment='left', verticalalignment='bottom' )
    otherNotes = 'Grey points to not meet TSM thresholds'
    fig.text( 0.08, subtitleY-dySubTitle, otherNotes, c='black', \
              fontsize=14, horizontalalignment='left', verticalalignment='top' )
    otherNotes = 'No bright limits have been applied\n'
    otherNotes += 'Numbers in square brackets are TSM values'    
    fig.text( 0.08, subtitleY-2.7*dySubTitle, otherNotes, c='black', \
              fontsize=14, horizontalalignment='left', verticalalignment='top' )

    dx = 0.02*( xLines.max()-xLines.min() )
    dy = 0.03*( yLines.max()-yLines.min() )
    ax.set_xlim( [ xLines.min()-dx, xLines.max()+dx ] )
    ax.set_ylim( [ yLines.min()-dy, yLines.max()+dy ] )

    return fig, ax


def formatAxisTicks( ax ):
    tick_fs = 14
    tl = 10
    tlm = 5
    tw = 2
    ax.spines['bottom'].set_linewidth( tw )
    ax.spines['left'].set_linewidth( tw )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params( labelsize=tick_fs )
    ax.xaxis.set_tick_params( length=tl, width=tw, which='major' )
    ax.xaxis.set_tick_params( length=tlm, width=tw, which='minor' )
    ax.yaxis.set_tick_params( length=tl, width=tw, which='major' )
    ax.yaxis.set_tick_params( length=tlm, width=tw, which='minor' )
    return ax

def addTopTSMs( ax, pl, TSM, TeqK, RpRE, TstarK, Tgrid, Rgrid, \
                xLines, yLines, survey={} ):
        
    framework = survey['framework']
    nx = len( xLines )-1
    ny = len( yLines )-1
    n = len( pl )
    ixs0 = np.arange( n )
    text_fs = 12
    nList = 5
    ms = 8
    for i in range( nx ): # loop over temperature columns
        ixsi = ( TeqK>=Tgrid[i] )*( TeqK<Tgrid[i+1] )
        #xtxt = xLines[i] + 0.05*( xLines[i+1]-xLines[i] )
        #xsymb = xLines[i] + 0.90*( xLines[i+1]-xLines[i] )
        xsymb = xLines[i] + 0.06*( xLines[i+1]-xLines[i] )
        xtxt = xLines[i] + 0.12*( xLines[i+1]-xLines[i] )
        #print( '\n\nNew cell:' )
        for j in range( ny ): # loop over radius rows
            # Indices inside box above the TSM threshold:
            RpREj = 0.5*( Rgrid[j]+Rgrid[j+1] )
            TSMj, TSMstr = survey['thresholdTSM']( RpREj, framework=framework )
            #if ( thresholdTSM=='ACWG' )+( thresholdTSM=='TOIs' ):
            #    TSMj = Utils.getThresholdTSM( RpREj, framework=thresholdTSM )
            #else:
            #    TSMj = thresholdTSM
            ixsj = ( RpRE>=Rgrid[j] )*( RpRE<Rgrid[j+1] )
            #ixsij = ixs0[ixsi*ixsj*( TSM>0)]
            ixsij = ixs0[ixsi*ixsj*( TSM>TSMj)]
            nij = len( ixsij )
            if nij>0:
                # Order by decreasing TSM:
                ixso = np.argsort( TSM[ixsij] )
                ixs = ixsij[ixso][::-1]
                nwrite = min( [ nij, nList ] )
                dy = ( yLines[j+1]-yLines[j] )/float(nList+0.5)
                y0 = yLines[j]+4.8*dy
                for k in range( nwrite ):
                    ytxt = y0-k*dy
                    plStr = pl[ixs][k].replace( ' ', '' )
                    plStr = '{0} [{1:.0f}]'.format( plStr, TSM[ixs][k] )
                    if ( plStr.find( '(APC' )>0 )+( plStr.find( '(CP)' )>0 ):
                        c = 'Silver'
                        wt = 'normal'
                    else:
                        c = 'Black'
                        wt = 'normal'
                    ax.text( xtxt, ytxt, plStr, fontsize=text_fs, weight=wt, color=c, \
                             horizontalalignment='left', verticalalignment='center' )
                    ck = Utils.getStarColor( TstarK[ixs][k] )
                    ax.plot( [xsymb], [ytxt], 'o', ms=ms, mec=ck, mfc=ck )
                #pdb.set_trace()
    return ax, TSMstr


def generateAxisScatter( xlim=[0,3100], ylim=[0,26], wideFormat=False, \
                         whichType='RpTeq', titleStr='', DecStr='', \
                         showLegend=True ):
    fig, ax, axLegend = generateAxes( wideFormat=wideFormat, whichType=whichType, \
                                      showLegend=showLegend )
    title_fs = 18
    toplineY = 0.98
    fig.text( 0.02, toplineY-0.02, titleStr, fontsize=title_fs, weight='heavy', \
              rotation=0, horizontalalignment='left', verticalalignment='bottom' )
    return fig, ax, axLegend


def generateAxisGrid( xlim=[0,3100], ylim=[0,26], wideFormat=False, whichType='RpTeq', \
                      RADecStr='', titleStr='', showLegend=True ):
    fig, ax, axLegend = generateAxes( wideFormat=wideFormat, whichType=whichType, \
                                      showLegend=showLegend )
    title_fs = 18
    toplineY = 0.98
    fig.text( 0.02, toplineY-0.02, titleStr, fontsize=title_fs, weight='heavy', \
              rotation=0, horizontalalignment='left', verticalalignment='bottom' )
    subtitle_fs = 14
    fig.text( 0.98, toplineY, RADecStr, fontsize=subtitle_fs, weight='normal', \
              rotation=0, horizontalalignment='right', verticalalignment='top' )
    return fig, ax, axLegend


def generateAxes( wideFormat=True, whichType='RpTeq', showLegend=True ):
    if wideFormat==False:
        fig = plt.figure( figsize=[11,9] )
        xlow = 0.09
        ylow = 0.085
        axh = 0.8
        axw = 0.90
        dxl = 0.06
        xlow2 = xlow+0.5*axw
        ylow2 = ylow+axh+0.005
        axw2 = 0.5*axw
        subtitleY = 0.94
        dyNewLine = 0.01
    else:
        fig = plt.figure( figsize=[16,9] )
        xlow = 0.064
        ylow = 0.085
        axh = 0.715
        axw = 0.93
        dxl = 0.044
        xlow2 = xlow+0.7*axw
        ylow2 = ylow+axh+0.02
        axw2 = 0.25*axw
        subtitleY = 0.925
        dySubTitle = 0.015
    ax = fig.add_axes( [ xlow, ylow, axw, axh ] )
    axLegend = fig.add_axes( [ xlow2, ylow2, axw2, 0.09*axh ] )
    addStellarSpectralTypeLegend( axLegend, ms=8, text_fs=10 )
    ax = formatAxes( ax, whichType=whichType )
    label_fs = 16
    if whichType=='RpTeq':
        fig.text( xlow-dxl, ylow+0.5*axh, '$R_p$ ($R_E$)', fontsize=label_fs, \
                  rotation=90, horizontalalignment='right', verticalalignment='center' )
        fig.text( xlow+0.5*axw, 0.001, '$T_{\\rm{eq}}$ (K)', fontsize=label_fs, \
                  rotation=0, horizontalalignment='center', verticalalignment='bottom' )
    elif ( whichType=='RpInsolLog' )+( whichType=='RpInsol' ):
        fig.text( xlow-dxl, ylow+0.5*axh, '$R_p$ ($R_E$)', fontsize=label_fs, \
                  rotation=90, horizontalalignment='right', verticalalignment='center' )
        fig.text( xlow+0.5*axw, 0.001, 'Insolation relative to Earth', rotation=0,
                  fontsize=label_fs, horizontalalignment='center', \
                  verticalalignment='bottom' )
    if showLegend==True:
        subtitle_fs = 14
        subtitleStr = 'Circles indicate host star spectral type'
        if wideFormat==True:
            sptY = 0.87
        else:
            sptY = 0.90
        fig.text( xlow2+0.5*axw2, sptY, subtitleStr, fontsize=subtitle_fs, \
                  horizontalalignment='center', verticalalignment='bottom', \
                  weight='normal', rotation=0 )
    return fig, ax, axLegend


def formatAxes( ax, whichType='RpTeq', xlim='default', ylim='default', \
                xticksMajor='default', yticksMajor='default', \
                xticksMinor='default', yticksMinor='default' ):
                
    tick_fs = 14
    tl = 10
    tlm = 5
    tw = 2
    if whichType=='RpTeq':
        if xlim is 'default':
            xlim = [ 0, 3100 ]
        if ylim is 'default':
            ylim = [ 0, 26 ]            
        if xticksMinor is 'default':
            xticksMinor = np.arange( 0, 4000, 100 )
        if xticksMajor is 'default':
            xticksMajor = np.arange( 0, 4000, 500 )
        if yticksMinor is 'default':
            yticksMinor = np.arange( 0, 30, 1 )
        if yticksMajor is 'default':
            yticksMajor = np.arange( 0, 30, 2 )
    elif whichType=='RpInsol':
        if xlim is 'default':
            xlim = [ 0, 11e3 ]
        if ylim is 'default':
            ylim = [ 0, 26 ]            
        if xticksMinor is 'default':
            xticksMinor = np.arange( 0, 10500, 100 )
        if xticksMajor is 'default':
            xticksMajor = np.arange( 0, 10500, 1000 )            
        if yticksMinor is 'default':
            yticksMinor = np.arange( 0, 30, 1 )
        if yticksMajor is 'default':
            yticksMajor = np.arange( 0, 30, 2 )
    elif whichType=='RpInsolLog':
        ax.set_xscale( 'log' )
        ax.xaxis.set_major_formatter( matplotlib.ticker.FuncFormatter( tickLogFormat ) )
        if xlim is 'default':
            xlim = [ 0.1, 11e3 ]
        if ylim is 'default':
            ylim = [ 0, 26 ]            
        if xticksMinor is 'default':
            #xticksMinor = np.arange( 100, 10500, 100 )
            xticksMinor = np.arange( 0.1, 1, 0.1 )
            xticksMinor = np.concatenate( [ xticksMinor, np.arange( 1, 10, 1 ) ] )
            xticksMinor = np.concatenate( [ xticksMinor, np.arange( 10, 100, 10 ) ] )
            xticksMinor = np.concatenate( [ xticksMinor, np.arange( 100, 1000, 100 ) ] )
            xticksMinor = np.concatenate( [ xticksMinor, np.arange( 1000, 10000, 1000 ) ] )
        if xticksMajor is 'default':
            #xticksMajor = np.arange( 100, 10500, 1000 )
            xticksMajor = np.logspace( -1, 4, 6 )
        if yticksMinor is 'default':
            yticksMinor = np.arange( 0, 30, 1 )
        if yticksMajor is 'default':
            yticksMajor = np.arange( 0, 30, 2 )
    else:
        pdb.set_trace()

    ax.xaxis.set_ticks( xticksMinor, minor=True )
    ax.xaxis.set_ticks( xticksMajor, minor=False )
    ax.yaxis.set_ticks( yticksMinor, minor=True )
    ax.yaxis.set_ticks( yticksMajor, minor=False )
        
    ax.tick_params( labelsize=tick_fs )
    ax.xaxis.set_tick_params( length=tl, width=tw, which='major' )
    ax.xaxis.set_tick_params( length=tlm, width=tw, which='minor' )
    ax.yaxis.set_tick_params( length=tl, width=tw, which='major' )
    ax.yaxis.set_tick_params( length=tlm, width=tw, which='minor' )

    ax.spines['bottom'].set_linewidth( tw )
    ax.spines['left'].set_linewidth( tw )
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim( ylim )
    ax.set_xlim( xlim )
    return ax

def tickLogFormat( y, pos ):
    # Find the number of decimal places required
    decimalplaces = int(np.ceil(np.maximum(-np.log10(y),0))) # =0 for numbers >=1
    # Insert that number into a format string
    formatstring = '{{:.{:1d}f}}'.format(decimalplaces)
    # Return the formatted tick label
    return formatstring.format(y)


    

def plotTeqRpScatter( planetNames, Teq, RpVal, Ts, TSM, TESS, ms=8, alpha=1, \
                      starColors=True, applyTSMcuts=False, survey={}, \
                      showGrid=True, indicateTESS=False, showSolarSystem=False, \
                      showStellarTrack=False, showNeptuneRadius=True, \
                      showJupiterRadius=True, titleStr='', wideFormat=False, \
                      dateStr='' ):
    """
    """
    n = len( Teq )
    nTESS = np.sum( TESS )
    cTESS = np.array( [ 213, 128, 255 ] )/256.
    zTESS = 2*n
    z0 = 200
    msTESS = 2*ms
    
    #fig, ax, ax2 = generateBasicAxis( wideFormat=wideFormat, titleStr=titleStr )
    fig, ax, ax2 = generateAxisScatter( wideFormat=wideFormat, titleStr=titleStr )
    if showGrid==True:
        Tgrid, Rgrid = drawGrid( ax, survey=survey, zorder=1 )
    framework = survey['framework']
        
    c0 = 0.9*np.ones( 3 )
    for i in range( n ): # loop over each exoplanet
        if starColors==True:
            c = Utils.getStarColor( Ts[i] )
        else:
            c = c0
        TSMi, TSMstr = survey['thresholdTSM']( RpVal[i], framework=framework )
        if applyTSMcuts==False: # plotting everything regardless of TSM
            if ( indicateTESS==True )*( TESS[i]==1 ):
                ax.plot( [Teq[i]], [RpVal[i]], 'o', ms=ms, alpha=alpha, \
                         mfc=c, mec=c, zorder=zTESS+nTESS+i )
                ax.plot( [Teq[i]], [RpVal[i]], 'o', ms=msTESS, alpha=alpha, \
                         zorder=zTESS, mfc=cTESS, mec=cTESS )
            else:
                ax.plot( [Teq[i]], [RpVal[i]], 'o', ms=ms, alpha=alpha, \
                         mfc=c, mec=c, zorder=z0+i )
                
        elif TSM[i]>TSMi: # if TSM cuts applied, this one is high enough
            if ( indicateTESS==True )*( TESS[i]==1 ):
                ax.plot( [Teq[i]], [RpVal[i]], 'o', ms=ms, alpha=alpha, \
                         mfc=c, mec=c, zorder=zTESS+nTESS+i )
                ax.plot( [Teq[i]], [RpVal[i]], 'o', ms=msTESS, alpha=alpha, \
                         zorder=zTESS, mfc=cTESS, mec=cTESS )
            else:
                ax.plot( [Teq[i]], [RpVal[i]], 'o', ms=ms, alpha=alpha, \
                         mfc=c, mec=c, zorder=z0+i )
        else: # otherwise plot as a smaller background point
            ax.plot( [Teq[i]], [RpVal[i]], 'o', ms=0.5*ms, alpha=alpha, \
                     mfc=c0, mec=c0, zorder=0 )
    if showStellarTrack==True:
        ax = addStellarTrack( ax )
    ax = addSolarSystem( ax, showSolarSystem=showSolarSystem, \
                         showNeptuneRadius=showNeptuneRadius, \
                         showJupiterRadius=showJupiterRadius )
    if wideFormat==False:
        subtitleY = 0.94
        dyNewLine = 0.01
    else:
        subtitleY = 0.925
        dySubTitle = 0.015
    if applyTSMcuts==True:
        fig.text( 0.08, subtitleY, TSMstr, c='green', fontsize=14, \
                  horizontalalignment='left', verticalalignment='bottom' )
        otherNotes = 'Grey points do not meet TSM thresholds'
        fig.text( 0.08, subtitleY-dySubTitle, otherNotes, c='black', \
                  fontsize=14, horizontalalignment='left', verticalalignment='top' )
        
    return fig, ax

def addStellarTrack( ax ):
    Rs0, Ts0 = Utils.readStellarTrack()
    Ts = np.linspace( Ts0.min(), Ts0.max(), 1000 )
    Rs = np.interp( Ts, Ts0, Rs0 )
    ixs1 = ( Ts>2200 )
    ixs2 = ( Ts<=2200 )*( Ts>1300 )
    ixs3 = ( Ts<=1300 )
    lw = 5
    c1 = np.array( [178,24,43] )/256. #Cyan
    c2 = np.array( [191,129,45] )/256. #'HotPink'
    c3 = np.array( [140,81,10] )/256. #'GreenYellow'
    ax.plot( Ts[ixs1], Rs[ixs1], '-', c=c1, lw=lw, zorder=1e5 )
    ax.plot( Ts[ixs2], Rs[ixs2], '-', c=c2, lw=lw, zorder=1e5 )
    ax.plot( Ts[ixs3], Rs[ixs3], '-', c=c3, lw=lw, zorder=1e5 )
    return ax

def addSolarSystem( ax, showSolarSystem=True, showNeptuneRadius=True, \
                    showJupiterRadius=True ):
    ss = Utils.solarSystem()
    if showSolarSystem==True:
        #ssc = 'Black'
        ssc = 'HotPink'
        ssc = np.array( [0,204,0] )/256.
        ssPlanets = list( ss['TeqK'].keys() )
        ssms = 8
        for k in ssPlanets:
            ax.plot( [ss['TeqK'][k]], [ss['RpRE'][k]], 'o', ms=ssms, \
                     zorder=1e6, mfc=ssc, mec=ssc, alpha=0.5 )
            ax.plot( [ss['TeqK'][k]], [ss['RpRE'][k]], 'o', ms=ssms, \
                     zorder=1e6, mfc='none', mec=ssc, mew=2, alpha=1 )
    cgrey = 0.6*np.ones( 3 )
    if showNeptuneRadius==True:        
        ax.axhline( ss['RpRE']['Neptune'], ls='--', c=cgrey, lw=2, zorder=0 )
    if showJupiterRadius==True:
        ax.axhline( ss['RpRE']['Jupiter'], ls='--', c=cgrey, lw=2, zorder=0 )
    return ax

def addStellarSpectralTypeLegend( ax, ms=10, text_fs=12 ):
    c, SpTs = Utils.getAllStarColors()
    for i in ['top','bottom','right','left']:
        ax.spines[i].set_visible(False)
    plt.setp( ax.xaxis.get_ticklabels(), visible=False )
    plt.setp( ax.yaxis.get_ticklabels(), visible=False )
    n = len( SpTs )
    for i in range( n ):
        k = SpTs[i]
        ax.plot( i+1, 0.5, 'o', ms=ms, mfc=c[k], mec=c[k] )
        ax.text( i+1, 0.4, k, rotation=45, fontsize=text_fs, \
                 horizontalalignment='right', verticalalignment='top' )
    ax.set_ylim( [ 0, 1 ] )
    ax.tick_params(axis = 'x', which = 'both', bottom = False, top = False)
    ax.tick_params(axis = 'y', which = 'both', left = False, right = False)
    return None
        


def readConfirmedProperties( ipath='confirmedProperties.pkl' ):
    """
    Returns properties for all confirmed planets. For planets with
    available radius but no mass measurement, an empirical relation
    is used to estimate the mass. This has been done for the purpose
    of identifying targets for a UV host star survey, which is why
    it's OK if the planet masses haven't been published.
    """
    ifile = open( ipath, 'rb' )
    z0 = pickle.load( ifile )
    ifile.close()
    z = z0['allVals']
    planetName = z['planetName']
    RsRS = z['RsRS']
    aAU = z['aAU']
    TeqK = z['TeqK']
    Insol = z['Insol']
    TstarK = z['TstarK']
    T14hr = z['T14hr']
    TSM = z['TSM']
    b = z['b']
    RpRs = z['RpRs']
    Vmag = z['Vmag']
    Jmag = z['Jmag']
    Hmag = z['Hmag']
    Kmag = z['Kmag']
    RpValRE = z['RpValRE']
    RpLsigRE = z['RpLowErrRE']
    RpUsigRE = z['RpUppErrRE']
    MpValME = z['MpValME']
    MpLsigME = z['MpLowErrME']
    MpUsigME = z['MpUppErrME']
    TESS = z['discoveredByTESS']
    n0 = len( planetName )
    ixs = np.arange( n0 )[np.isnan( MpValME )*np.isfinite( RpValRE )]
    n = len( ixs )
    for i in range( n ):
        MpValME[ixs[i]] = Utils.planetMassFromRadius( RpValRE[ixs[i]] )
        MpUncFill = 1e9 #min( [ 0.1, MpValME[ixs[i]]/5.01 ] )
        MpLsigME[ixs[i]] = MpUncFill
        MpUsigME[ixs[i]] = MpUncFill
    TSM[ixs] = Utils.computeTSM( RpValRE[ixs], MpValME[ixs], \
                                 RsRS[ixs], TeqK[ixs], Jmag[ixs] )
    ixs = np.isfinite( TeqK )*np.isfinite( TSM )*np.isfinite( RpValRE )
    print( '\nReading in {0:.0f} planets total.'.format( n0 ) )
    print( 'Returning {0:.0f} planets with radii, TSM, and Teq values.'\
           .format( ixs.sum() ) )
    outp = { 'planetName':planetName[ixs], 'TESS':TESS[ixs], \
             'TSM':TSM[ixs], 'T14hr':T14hr[ixs], 'Vmag':Vmag[ixs], \
             'Jmag':Jmag[ixs], 'Hmag':Hmag[ixs], 'Kmag':Kmag[ixs], \
             'b':b[ixs], 'RpRs':RpRs[ixs], 'TeqK':TeqK[ixs], 'Insol':Insol[ixs], \
             'TstarK':TstarK[ixs], 'RsRS':RsRS[ixs], 'aAU':aAU[ixs], \
             'RpValRE':RpValRE[ixs], 'RpLsigRE':RpLsigRE[ixs], 'RpUsigRE':RpUsigRE[ixs], \
             'MpValME':MpValME[ixs], 'MpLsigME':MpLsigME[ixs], 'MpUsigME':MpUsigME[ixs] }
    return outp, z0['dateStr']


def readTOIProperties( ipath='toiProperties.pkl' ):
    """
    """
    ifile = open( ipath, 'rb' )
    z0 = pickle.load( ifile )
    ifile.close()
    z = z0['allVals']
    planetName = z['planetName']
    #print( planetName )
    #pdb.set_trace()
    RA = z['RA_deg']
    RAhr = RA*(24/360.)
    Dec = z['Dec_deg']
    n0 = len( planetName )
    RsRS = z['RsRS']
    TeqK = z['TeqK']
    Jmag = z['Jmag']
    TstarK = z['TstarK']
    TSM = z['TSM']
    RpValRE = z['RpValRE']
    MpValME = z['MpValME']
    RpRs = z['RpRs']
    Jmag = z['Jmag']
    #planetName = []
    #n = len( toiNumbers )
    #for i in range( n ):
    #    planetName += [ 'TOI-{0}'.format( toiNumbers[i] ) ]
    #planetName = np.array( planetName, dtype=str )
    ixs = np.isfinite( TeqK )*np.isfinite( TSM )*np.isfinite( RpValRE )
    print( '\nReading in {0:.0f} TOIs total.'.format( n0 ) )
    print( 'Returning {0:.0f} TOIs with radii, TSM, and Teq values.'\
           .format( ixs.sum() ) )
    outp = { 'planetName':planetName[ixs], 'TSM':TSM[ixs], 'RpRs':RpRs[ixs], \
             'RA_deg':RA[ixs], 'RA_hr':RAhr[ixs], 'Dec_deg':Dec[ixs], \
             'TeqK':TeqK[ixs], 'TstarK':TstarK[ixs], 'RsRS':RsRS[ixs], 'Jmag':Jmag[ixs], \
             'RpValRE':RpValRE[ixs], 'MpValME':MpValME[ixs] }
    return outp, z0['dateStr']

def readNoMassTESSProperties():
    """
    Returns properties of confirmed TESS planets lacking published masses.
    """
    ifile = open( 'confirmedProperties.pkl', 'rb' )
    z0 = pickle.load( ifile )
    ifile.close()
    z = z0['allVals']
    planetName = z['planetName']
    RsRS = z['RsRS']
    TeqK = z['TeqK']
    Jmag = z['Jmag']
    TstarK = z['TstarK']
    TSM = z['TSM']
    RpValRE = z['RpValRE']
    RpLsigRE = z['RpLowErrRE']
    RpUsigRE = z['RpUppErrRE']
    MpValME = z['MpValME']
    MpLsigME = z['MpLowErrME']
    MpUsigME = z['MpUppErrME']
    TESS = z['discoveredByTESS']
    print( '\nReading confirmed TESS planets lacking peer-reviewed published masses:' )
    ixs = np.isfinite( TeqK )*np.isfinite( RpValRE )*np.isnan( MpValME )*( TESS==1 )
    n = ixs.sum()
    for i in range( n ):
        print( '{0:.0f}. {1}'.format( i+1, planetName[ixs][i] ) )
    #print( '\nReading in {0:.0f} planets total.'.format( len( TSM ) ) )
    print( 'Returning {0:.0f} planets with measured radii and Teq values.'\
           .format( ixs.sum() ) )
    # Add in the estimated masses and then use these to compute estimated TSM values
    # (ESM values should already be included as it doesn't require the mass):
    MpValME = np.zeros( n )
    for i in range( n ):
        MpValME[i] = Utils.planetMassFromRadius( RpValRE[ixs][i] )
    TSM = Utils.computeTSM( RpValRE[ixs], MpValME, RsRS[ixs], TeqK[ixs], Jmag[ixs] )
    print( 'Masses taken from empirical relation and used for TSM calculation.\n' )
    outp = { 'planetName':planetName[ixs], 'TESS':TESS[ixs], \
             'TeqK':TeqK[ixs], 'TstarK':TstarK[ixs], 'RsRS':RsRS[ixs], \
             'RpValRE':RpValRE[ixs], 'RpLsigRE':RpLsigRE[ixs], 'RpUsigRE':RpUsigRE[ixs], \
             'MpValME':MpValME, 'TSM':TSM }
    return outp


def readConfirmedTESSProperties( publishedMasses=True ):
    """
    Returns properties of confirmed TESS planets.
    """
    ifile = open( 'confirmedProperties.pkl', 'rb' )
    z0 = pickle.load( ifile )
    ifile.close()
    z = z0['allVals']
    planetName = z['planetName']
    RsRS = z['RsRS']
    TeqK = z['TeqK']
    Jmag = z['Jmag']
    TstarK = z['TstarK']
    TSM = z['TSM']
    RpValRE = z['RpValRE']
    RpLsigRE = z['RpLowErrRE']
    RpUsigRE = z['RpUppErrRE']
    MpValME = z['MpValME']
    MpLsigME = z['MpLowErrME']
    MpUsigME = z['MpUppErrME']
    TESS = z['discoveredByTESS']
    if publishedMasses==True:
        print( '\nReading confirmed TESS planets with peer-reviewed published masses:' )
        ixs = np.isfinite( TeqK )*np.isfinite( RpValRE )*np.isfinite( MpValME )*( TESS==1 )
        n = ixs.sum()
        for i in range( n ):
            print( '{0:.0f}. {1}'.format( i+1, planetName[ixs][i] ) )
        #print( '\nReading in {0:.0f} planets total.'.format( len( TSM ) ) )
        print( 'Returning {0:.0f} planets with measured radii, TSM, and Teq values.'\
               .format( ixs.sum() ) )

        outp = { 'planetName':planetName[ixs], 'TESS':TESS[ixs], \
                 'TeqK':TeqK[ixs], 'TstarK':TstarK[ixs], 'RsRS':RsRS[ixs], \
                 'RpValRE':RpValRE[ixs], 'RpLsigRE':RpLsigRE[ixs], \
                 'RpUsigRE':RpUsigRE[ixs], 'MpValME':MpValME[ixs], 'TSM':TSM[ixs] }        
    else:
        print( '\nReading confirmed TESS planets lacking peer-reviewed published masses:' )
        ixs = np.isfinite( TeqK )*np.isfinite( RpValRE )*np.isnan( MpValME )*( TESS==1 )
        n = ixs.sum()
        for i in range( n ):
            print( '{0:.0f}. {1}'.format( i+1, planetName[ixs][i] ) )
        #print( '\nReading in {0:.0f} planets total.'.format( len( TSM ) ) )
        print( 'Returning {0:.0f} planets with measured radii and Teq values.'\
               .format( ixs.sum() ) )
        # Add in the estimated masses and then use these to compute estimated TSM values
        # (ESM values should already be included as it doesn't require the mass):
        MpValME = np.zeros( n )
        for i in range( n ):
            MpValME[i] = Utils.planetMassFromRadius( RpValRE[ixs][i] )
        TSM = Utils.computeTSM( RpValRE[ixs], MpValME, RsRS[ixs], \
                                TeqK[ixs], Jmag[ixs] )
        print( 'Masses taken from empirical relation and used for TSM calculation.\n' )
        outp = { 'planetName':planetName[ixs], 'TESS':TESS[ixs], \
                 'TeqK':TeqK[ixs], 'TstarK':TstarK[ixs], 'RsRS':RsRS[ixs], \
                 'RpValRE':RpValRE[ixs], 'RpLsigRE':RpLsigRE[ixs], \
                 'RpUsigRE':RpUsigRE[ixs], 'MpValME':MpValME, 'TSM':TSM }
    return outp


def readWithMassTESSProperties():
    """
    Returns properties of confirmed TESS planets lacking published masses.
    """
    ifile = open( 'confirmedProperties.pkl', 'rb' )
    z0 = pickle.load( ifile )
    ifile.close()
    z = z0['allVals']
    planetName = z['planetName']
    RsRS = z['RsRS']
    TeqK = z['TeqK']
    Jmag = z['Jmag']
    TstarK = z['TstarK']
    TSM = z['TSM']
    RpValRE = z['RpValRE']
    RpLsigRE = z['RpLowErrRE']
    RpUsigRE = z['RpUppErrRE']
    MpValME = z['MpValME']
    MpLsigME = z['MpLowErrME']
    MpUsigME = z['MpUppErrME']
    TESS = z['discoveredByTESS']
    print( '\nReading confirmed TESS planets lacking peer-reviewed published masses:' )
    ixs = np.isfinite( TeqK )*np.isfinite( RpValRE )*np.isfinite( MpValME )*( TESS==1 )
    n = ixs.sum()
    for i in range( n ):
        print( '{0:.0f}. {1}'.format( i+1, planetName[ixs][i] ) )
    #print( '\nReading in {0:.0f} planets total.'.format( len( TSM ) ) )
    print( 'Returning {0:.0f} planets with measured radii, TSM, and Teq values.'\
           .format( ixs.sum() ) )
    
    outp = { 'planetName':planetName[ixs], 'TESS':TESS[ixs], \
             'TeqK':TeqK[ixs], 'TstarK':TstarK[ixs], 'RsRS':RsRS[ixs], \
             'RpValRE':RpValRE[ixs], 'RpLsigRE':RpLsigRE[ixs], 'RpUsigRE':RpUsigRE[ixs], \
             'MpValME':MpValME[ixs], 'TSM':TSM[ixs] }
    return outp


def readPredictedProperties( version=2 ):
    iname = 'predictedProperties_v{0:.0f}.pkl'.format( version )
    ifile = open( iname, 'rb' )
    z = pickle.load( ifile )
    ifile.close()
    RsRS = z['RsRS']
    TeqK = z['TeqK']
    Insol = z['Insol']
    TstarK = z['TstarK']
    TSM = z['TSM']
    aRs = z['aRs']
    b = z['b']
    ideg = np.rad2deg( np.arccos( b/aRs ) )
    T14hr = z['T14hr']
    Pday = z['Pday']
    MsMS = z['MsMS']
    RpValRE = z['RpValRE']
    MpValME = z['MpValME']
    cad2min = z['cad2min']
    Vmag = z['Vmag']
    Jmag = z['Jmag']
    Kmag = z['Kmag']
    ixs = np.isfinite( TeqK )*np.isfinite( TSM )*np.isfinite( RpValRE )
    print( '\nReading in {0:.0f} planets total.'.format( len( TSM ) ) )
    print( 'Returning {0:.0f} planets with radii, TSM, and Teq values.'\
           .format( ixs.sum() ) )
    outp = { 'TSM':TSM[ixs], 'cad2min':cad2min, 'TeqK':TeqK[ixs], \
             'aRs':aRs[ixs], 'Pday':Pday[ixs], 'Insol':Insol[ixs], \
             'MsValMS':MsMS[ixs], 'RsRS':RsRS[ixs], 'TstarK':TstarK[ixs], \
             'RpValRE':RpValRE[ixs], 'MpValME':MpValME[ixs], \
             'ideg':ideg[ixs], 'b':b[ixs], 'T14hr':T14hr[ixs], \
             'Vmag':Vmag[ixs], 'Jmag':Jmag[ixs], 'Kmag':Kmag[ixs] }
    return outp



