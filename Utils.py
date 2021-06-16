import pdb, sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import pysynphot
import math

"""
readStellarTrack()
planetMassFromRadius()
solarSystem()
computeTSM()
"""

RSUN_SI = 6.955e8
MSUN_SI = 1.9889e30
MJUP_SI = 1.8986e27
MEARTH_SI = 5.9723e24
RJUP_SI = 7.149e7
REARTH_SI = 6.371e6
AU_SI = 1.496e11
GRAV_SI = 6.67428e-11 # gravitational constant in m^3 kg^-1 s^-2


def photBands():
    idir = os.path.dirname( __file__ )
    tessPath = os.path.join( idir, 'tess-response-function-v2.0.csv' )
    Tband = np.loadtxt( tessPath, delimiter=',' )
    Tband[:,0] /= 1000. # convert from nm to micron.
    Jband = np.loadtxt( os.path.join( idir, '2MASS_J.dat' ) )
    Jband[:,0] /= 1000. # convert from nm to micron.
    Hband = np.loadtxt( os.path.join( idir, '2MASS_H.dat' ) )
    Hband[:,0] /= 1000. # convert from nm to micron.
    Kband = np.loadtxt( os.path.join( idir, '2MASS_Ks.dat' ) )
    Kband[:,0] /= 1000. # convert from nm to micron.
    return Tband, Jband, Hband, Kband


def modelStellarSpectrum( TeffK, loggCGS, FeH=0 ):
    if TeffK<3500:
        print( '  WARNING: set Teff={0:.0f}K --> Teff=3500K'.format( TeffK ) )
        TeffK = 3500
    if loggCGS>5:
        print( '  WARNING: set logg={0:.3f} --> logg=5'.format( loggCGS ) )
        loggCGS = 5
    sp = pysynphot.Icat( 'k93models', TeffK, FeH, loggCGS )
    star = [ sp.wave/1e4, sp.flux ]
    return star


def convertTmag( Tmag, TeffK, loggCGS, outputMag='J', vega=None, star=None ):
    """
    Routine to convert Tmag to JHK mag.

    For example:

    Kmag = Tmag + 2.5*log10( [ vega_K/vega_T ]*[ star_T/star_K ] )

    where vega_K/vega_T is flux of Vega in the K band relative to T band
    and star_K/star_T is flux of Vegathe star of interest in the K band 
    relative to T band.
    """

    # Read in spectra for Vega and the star of interest:
    # NOTE: This is the most time-consuming step.    
    t1 = time.time()
    if vega is None:
        vega = spectrumVega( makePlot=False )
    t2 = time.time()
    if star is None:
        star = modelStellarSpectrum( TeffK, loggCGS, FeH=0 )
    t3 = time.time()
    # Read in the photometric passbands:
    T, J, H, Ks = photBands()
    if outputMag=='J':
        M = J
    elif outputMag=='H':
        M = H
    elif outputMag=='Ks':
        M = Ks

    # Interpolate the Vega spectrum onto the photometric passbands
    # and then sum to get the relative fluxes in each passband:
    t4 = time.time()
    Fvega_T = np.sum( T[:,1]*np.interp( T[:,0], vega[0], vega[1] ) )
    Fvega_M = np.sum( M[:,1]*np.interp( M[:,0], vega[0], vega[1] ) )
    
    # Do the same for the star of interest:
    t5 = time.time()
    Fstar_T = np.sum( T[:,1]*np.interp( T[:,0], star[0], star[1] ) )
    Fstar_M = np.sum( M[:,1]*np.interp( M[:,0], star[0], star[1] ) )

    t6 = time.time()
    # Use the relative brightnesses to convert from Tmag to the output magnitude:
    Mmag = Tmag + 2.5*np.log10( ( Fvega_M/Fvega_T )*( Fstar_T/Fstar_M ) )
    t7 = time.time()

    if 0:
        print( t2-t1 )
        print( t3-t2 )
        print( t4-t3 )
        print( t5-t4 )
        print( t6-t5 )
        print( t7-t6 )
        print( t7-t1 )
    return Mmag


def spectrumVega( makePlot=False ):
    sp = pysynphot.Icat( 'k93models', 9600, 0, 4.1 )
    wavMicr = sp.wave/1e4
    Flam = sp.flux
    if makePlot:
        # Compare to observed/model Vega from HST calibration:
        ipath = os.path.join( os.environ['PYSYN_CDBS'], 'calspec', \
                              'alpha_lyr_stis_010.fits' )
        #ipath = os.path.join( os.environ['PYSYN_CDBS'], 'calspec', \
        #'alpha_lyr_mod_004.fits' )
        hst = pysynphot.FileSpectrum( ipath )
        plt.figure()
        plt.plot( wavMicr, Flam/Flam.max(), '-k', label='Kurucz model' )
        plt.plot( hst.wave/1e4, hst.flux/hst.flux.max(), '-r', label='HST cal' )
        plt.xlim( [ 0, 2 ] )
        plt.title( 'Vega' )
    return wavMicr, Flam
    
    

def densityContours():
    idir = os.path.dirname( __file__ )
    h2 = np.loadtxt( os.path.join( idir, 'contours_h2.txt' ) )
    h2o = np.loadtxt( os.path.join( idir, 'contours_h2o.txt' ) )
    mgsio3 = np.loadtxt( os.path.join( idir, 'contours_mgsio3.txt' ), skiprows=1 )
    fe = np.loadtxt( os.path.join( idir, 'contours_fe.txt' ), skiprows=1 )
    uM = 1 # ( MEARTH_SI/MJUP_SI )
    uR = 1 # ( REARTH_SI/RJUP_SI )
    h2[:,0] = h2[:,0]*uM
    h2[:,1] = h2[:,1]*uR
    h2o[:,0] = h2o[:,0]*uM
    h2o[:,1] = h2o[:,1]*uR
    mgsio3[:,0] = mgsio3[:,0]*uM
    mgsio3[:,1] = mgsio3[:,1]*uR
    fe[:,0] = fe[:,0]*uM
    fe[:,1] = fe[:,1]*uR
    return h2, h2o, mgsio3, fe

def readStellarTrack():
    d = np.loadtxt( 'ames_dusty_5Gyr.txt' )
    MsSI = d[:,0]*MSUN_SI
    TeffK = d[:,1]
    loggCGS = d[:,3]
    gCGS = 10**loggCGS
    gSI = gCGS/100.
    RsSI = np.sqrt( GRAV_SI*MsSI/gSI )
    RsRE = RsSI/REARTH_SI
    return RsRE, TeffK
    

def massRadiusChenKipping2017( RpRE_in ):
    """
    Evaluates the mean of the Chen & Kipping (2017) distribution.

    NOTE:
    The S3 index has been adjusted to be slightly positive. This
    is done purely for convenience, because otherwise it is not
    possible to quickly obtain a deterministic mass for a given
    radius, i.e. when there's a combination of negative and 
    positive indices there will be mass degeneracies for certain
    input radii.
    """
    
    # Power law indices:
    S1 = 0.2790
    S2 = 0.589
    #S3 = -0.044 # value quoted in Chen & Kipping (2017)
    S3 = 0.01 # mild tweak done purely for convenience
    S4 = 0.881
    # Other transition points from Table 1
    T12ME = np.log10( 2.04 )
    T23ME = np.log10( 0.414*( MJUP_SI/MEARTH_SI ) )
    T34ME = np.log10( 0.080*( MSUN_SI/MEARTH_SI ) )
    # Terran power law constant from Table 1:
    C1curl = np.log10( 1.008 )
    # Iteratively derive other power law constants:
    C2curl = C1curl + ( S1-S2 )*T12ME
    C3curl = C2curl + ( S2-S3 )*T23ME
    C4curl = C3curl + ( S3-S4 )*T34ME

    log10MpME = np.linspace( -3, 5, 1000 )

    log10M12 = np.log10( 2.04 )
    log10M23 = np.log10( 0.414*( MJUP_SI/MEARTH_SI ) )
    log10M34 = np.log10( 0.080*( MSUN_SI/MEARTH_SI ) )
    ixs1 = ( log10MpME<=log10M12 )
    ixs2 = ( log10MpME>log10M12 )*( log10MpME<=log10M23 )
    ixs3 = ( log10MpME>log10M23 )*( log10MpME<=log10M34 )
    ixs4 = ( log10MpME>log10M34 )


    log10RpRE = np.ones_like( log10MpME )
    log10RpRE[ixs1] = C1curl + ( log10MpME[ixs1]*S1 )
    log10RpRE[ixs2] = C2curl + ( log10MpME[ixs2]*S2 )
    log10RpRE[ixs3] = C3curl + ( log10MpME[ixs3]*S3 )
    log10RpRE[ixs4] = C4curl + ( log10MpME[ixs4]*S4 )

    log10MpME_out = np.interp( np.log10( RpRE_in ), log10RpRE, log10MpME )
    MpME_out = 10**log10MpME_out

    if 0:
        plt.figure()
        plt.plot( log10MpME[ixs1], log10RpRE[ixs1], '-r' )
        plt.plot( log10MpME[ixs2], log10RpRE[ixs2], '-k' )
        plt.plot( log10MpME[ixs3], log10RpRE[ixs3], '-g' )
        plt.plot( log10MpME[ixs4], log10RpRE[ixs4], '-c' )
        #plt.plot( [log10MpME_out], np.log10( RpRE_in ), 'ok' )
        print( RpRE_in, MpME_out )
        pdb.set_trace()
        
    return MpME_out
    

def planetMassFromRadius( RpRE, whichRelation='Chen&Kipping2017' ):
    """
    Taken from Eq 2 of Kempton et al (2017).
    """
    if whichRelation=='Chen&Kipping2017':
        MpME = massRadiusChenKipping2017( RpRE )
    elif whichRelation=='Kempton+2018':
        if np.ndim( RpRE )==0:
            if RpRE<1.23:
                MpME = 0.9718*( RpRE**3.58 )
            elif ( RpRE>=1.23 )*( RpRE<30 ):
                MpME = 1.436*( RpRE**1.70 )
            else:
                print( '\n\nPlanet radius too large... {0:.1f}RE'.format( RpRE ) )
                MpME = np.nan
        else:
            MpME = np.zeros_like( RpRE )
            ixs1 = ( RpRE<1.23 )
            MpME[ixs1] = 0.9718*( RpRE[ixs1]**3.58 )
            ixs2 = ( RpRE>=1.23 )*( RpRE<14.26 )
            MpME[ixs2] = 1.436*( RpRE[ixs2]**1.70 )
            ixs3 = ( RpRE>=14.26 )
            MpME[ixs3] = np.nan
    return MpME


def solarSystem():
    z = {}
    z['TstarK'] = 5800.
    z['RsSI'] = RSUN_SI
    z['aAU'] = {}
    z['aAU']['Mercury'] = 0.3870993
    z['aAU']['Venus'] = 0.723336
    z['aAU']['Earth'] = 1.000003
    z['aAU']['Mars'] = 1.52371
    z['aAU']['Jupiter'] = 5.2029
    z['aAU']['Saturn'] = 9.537
    z['aAU']['Titan'] = z['aAU']['Saturn']
    z['aAU']['Uranus'] = 19.189
    z['aAU']['Neptune'] = 30.0699
    planets = list( z['aAU'].keys() )
    z['aSI'] = {}
    for k in planets:
        z['aSI'][k] = AU_SI*z['aAU'][k]
    z['aRs'] = {}
    for k in planets:
        z['aRs'][k] = z['aSI'][k]/z['RsSI']
    z['TeqK'] = {}
    for k in planets:
        z['TeqK'][k] = calcTeqK( z['TstarK'], z['aRs'][k] )
    z['RpSI'] = {}
    z['RpSI']['Mercury'] = ( 4879e3 )/2.
    z['RpSI']['Venus'] = ( 12104e3 )/2.
    z['RpSI']['Earth'] = ( 12756e3 )/2.
    z['RpSI']['Mars'] = ( 6792e3 )/2.
    z['RpSI']['Jupiter'] = ( 142984e3 )/2.
    z['RpSI']['Saturn'] = ( 120536e3 )/2.
    z['RpSI']['Titan'] = 2575e3
    z['RpSI']['Uranus'] = ( 51118e3 )/2.
    z['RpSI']['Neptune'] = ( 49528e3 )/2.
    z['RpRE'] = {}
    for k in planets:
        z['RpRE'][k] = z['RpSI'][k]/REARTH_SI
    #for k in planets:
    #    print( '{0}: Rp={1:.1f}, T={2:.0f}'.format( k, z['RpRE'][k], z['TeqK'][k] ) )
    return z


def computeTSM( RpValRE, MpValME, RsRS, TeqK, Jmag ):
    nAll = len( RpValRE )
    # Indices for different radii of each scale factor:
    ixsA = np.arange( nAll )[np.isfinite( RpValRE )]
    ixsB = np.arange( nAll )[np.isfinite( RpValRE )==False]
    ixs1 = ixsA[( RpValRE[ixsA]<1.5 )]
    ixs2 = ixsA[( RpValRE[ixsA]>=1.5 )*( RpValRE[ixsA]<2.75 )]
    ixs3 = ixsA[( RpValRE[ixsA]>=2.75 )*( RpValRE[ixsA]<4.0 )]
    ixs4 = ixsA[( RpValRE[ixsA]>=4.0 )*( RpValRE[ixsA]<10. )]
    ixs5 = ixsA[( RpValRE[ixsA]>=10 )]
    # Scale factors provided in Table 1 of Kempton et al (2018):
    c1 = 0.190
    c2 = 1.26
    c3 = 1.28
    c4 = 1.15
    c5 = 1.0
    # TSM before applying scale factor:
    Rp3 = RpValRE**3.
    MpRs2 = MpValME*( RsRS**2. )
    y = ( Rp3*TeqK/MpRs2 )*( 10**( -Jmag/5. ) )
    TSM = np.zeros( nAll )
    TSM[ixs1] = c1*y[ixs1]
    TSM[ixs2] = c2*y[ixs2]
    TSM[ixs3] = c3*y[ixs3]
    TSM[ixs4] = c4*y[ixs4]
    TSM[ixs5] = c5*y[ixs5]
    TSM[ixsB] = np.nan
    return TSM


def computeESM( TeqK, RpRs, TstarK, Kmag ):
    nAll = len( TeqK )
    wavRefSI = 7.5e-6
    TdayK = 1.10*TeqK # from Section 3.2 of Kempton+2018
    Bday =  PlanckFuncSI( wavRefSI, TdayK )
    Bstar =  PlanckFuncSI( wavRefSI, TstarK )
    ESM = 4.29*(1e6)*( Bday/Bstar )*( RpRs**2. )*( 10**( -Kmag/5. ) )
    return ESM


def PlanckFuncSI( wavSI, T ):
    """
    Returns the Planck spectrum in cgs units given 
    a wavelength range and temperature.
    """
    hSI = np.longdouble( 6.62607015e-34 ) # Planck constant (J*s)
    cSI = np.longdouble( 2.9979245800e8 ) # speed of light (m/s)
    kSI = np.longdouble( 1.380649e-23 ) # Boltzman constant (J/K)
    c0 =  2.*hSI*( cSI**2. )/( wavSI**5. )
    c1 = hSI*cSI/kSI/T
    irrSI = c0/( np.exp( c1/wavSI ) - 1. ) # radiance
    fluxSI = np.pi*irrSI*wavSI # surface flux
    return fluxSI


def calcTeqK( TstarK, aRs ):
    TeqK = ( TstarK/np.sqrt( aRs ) )*( 0.25**0.25 )
    return TeqK


def getThresholdTSM_REDUNDANT( RpRE, framework='ACWG' ):
    """
    Thresholds from Figure 5 of Kempton et al. (2018).
    """
    if framework=='ACWG':
        if RpRE<1.50: # 1. Terrestrials
            return 10
        elif ( RpRE>=1.50 )*( RpRE<2.75 ): # 2. Small sub-Neptunes
            return 90
        elif ( RpRE>=2.75 )*( RpRE<4.00 ): # 3. Large sub-Neptunes
            return 90
        elif ( RpRE>=4.00 )*( RpRE<10.00 ): # 4. Sub-Jovians
            return 90
        elif ( RpRE>=10.00 ): # 5. Gas giants
            return 100
    elif framework=='TOIs':
        if RpRE<1.50: # 1. Terrestrials
            return 10
        else:
            return 50 # 2. Everything else

    
def getTSMStr_REDUNDANT( thresholdTSM ):
    if thresholdTSM=='ACWG':
        TSMStr = '* Kempton et al. (2018) TSM cuts applied'
        #pdb.set_trace()
    elif thresholdTSM=='TOIs':
        TSMStr = 'Only targets with TSM>50 (Rp>1.5RE) or TSM>10 (Rp<1.5RE) shown'
    else:
        TSMStr = 'Only targets with TSM>{0:.0f} shown'.format( thresholdTSM )
    return TSMStr


def getRARanges():
    m = 4
    RAedges = np.arange( 0, 24+m, m )
    n = len( RAedges )-1
    RARanges = []
    for i in range( n ):
        RARanges += [ [ RAedges[i], RAedges[i+1] ] ]
    return RARanges
    
def getRARange( month ):
    # RAmid is approximately the sidereal angle (i.e. overhead RA) at
    # midnight on the 20th day of the month.
    if month=='Jan':
        RAmid = 8
    elif month=='Feb':
        RAmid = 10
    elif month=='Mar':
        RAmid = 12
    elif month=='Apr':
        RAmid = 14
    elif month=='May':
        RAmid = 16
    elif month=='Jun':
        RAmid = 18
    elif month=='Jul':
        RAmid = 20
    elif month=='Aug':
        RAmid = 22
    elif month=='Sep':
        RAmid = 0
    elif month=='Oct':
        RAmid = 2
    elif month=='Nov':
        RAmid = 4
    elif month=='Dec':
        RAmid = 6
    dRA = 6 # +/- RA hr from overhead at midnight.
    RARange = [ RAmid-dRA, RAmid+dRA ]
    return RARange


def processRARestriction( RAMin_hr, RAMax_hr ):
    if ( RAMin_hr is not None )*( RAMax_hr is not None ):
        RAStr = '{0:.0f}<RA(hr)<{1:.0f}'.format( RAMin_hr, RAMax_hr )
    elif ( RAMin_hr is not None )+( RAMax_hr is not None ):
        if RAMin_hr is None:
            RAMin_hr = -1e9 
            RAStr = 'Only RA(hr)<{0:.1f} targets'.format( RAMax_hr )
        else:
            RAMax_hr = 1e9
            RAStr = 'Only RA(hr)>{0:.0f} targets'.format( RAMin_hr )
    else:
        RAMin_hr = -1e9 
        RAMax_hr = 1e9
        RAStr = 'No RA restrictions applied'
    return RAStr, RAMin_hr, RAMax_hr

def processDecRestriction( DecMin_deg, DecMax_deg ):
    if ( DecMin_deg is not None )*( DecMax_deg is not None ):
        DecStr = '{0:.0f}<Dec(deg)<{1:.0f}'.format( DecMin_deg, DecMax_deg )
    elif ( DecMin_deg is not None )+( DecMax_deg is not None ):
        if DecMin_deg is None:
            DecMin_deg = -1e9 
            DecStr = 'Only Dec(deg)<{0:.0f} targets'.format( DecMax_deg )
        else:
            DecMax_deg = 1e9
            DecStr = 'Only Dec(deg)>{0:.0f} targets'.format( DecMin_deg )
    else:
        DecMin_deg = -1e9 
        DecMax_deg = 1e9
        DecStr = 'No Dec restrictions applied'
    return DecStr, DecMin_deg, DecMax_deg


def getStarColor( T ):
    if ( T<3400 ): # late-M
        #c = np.array( [227,26,28] )/256.
        c = np.array( [178,24,43] )/256.
    elif ( T>=3400 )*( T<3800 ): # early-M
        c = np.array( [252,78,42] )/256.
    elif ( T>=3800 )*( T<4600 ): # late-K
        c = np.array( [253,141,60] )/256.
    elif ( T>=4600 )*( T<5200 ): # early-K
        c = np.array( [254,178,76] )/256.
    elif ( T>=5200 )*( T<5700 ): # late-G
        c = np.array( [254,217,118] )/256.
    elif ( T>=5700 )*( T<6000 ): # early-G
        c = np.array( [255,237,160] )/256.
    elif ( T>=6000 )*( T<6700 ): # late-F
        c = np.array( [158,202,225] )/256.
    elif ( T>=6700 )*( T<7400 ): # early-F
        c = np.array( [107,174,214] )/256.
    else: # OBA
        c = np.array( [8,81,156] )/256.
    return c
        
def getAllStarColors():
    c = {}
    c['late-M'] = getStarColor( 3200 )
    c['early-M'] = getStarColor( 3600 )
    c['late-K'] = getStarColor( 4000 )
    c['early-K'] = getStarColor( 4800 )
    c['late-G'] = getStarColor( 5500 )
    c['early-G'] = getStarColor( 5900 )
    c['late-F'] = getStarColor( 6500 )
    c['early-F'] = getStarColor( 7200 )
    c['OBA'] = getStarColor( 7500 )
    #SpTs = [ 'OBA', 'early-F', 'late-F', 'early-G', 'late-G', \
    #         'early-K', 'late-K', 'early-M', 'late-M' ]
    SpTs = [ 'late-M', 'early-M', 'late-K', 'early-K', \
             'late-G', 'early-G', 'late-F', 'early-F', 'OBA' ]
    return c, SpTs

def computeStellarMass( RsRS, loggstarCGS ):

    gStarSI = 10**loggstarCGS / 100 # converts log(g)[cm/s^2] to g [m/s^2]
    RsRS = RsRS * RSUN_SI # converts stellar radius from solar unit to SI 
    MsMS = gStarSI * RsRS**2 / GRAV_SI
    MsMS /= MSUN_SI # converts stellar mass to solar unit
    
    return MsMS

def computeRVSemiAmp( Pday, MpME, MsMS ):
    """
    Equation from: https://exoplanetarchive.ipac.caltech.edu/docs/poet_calculations.html
    """
    MpMJ = MpME * MEARTH_SI / MJUP_SI # converts mass from Earth unit to Jupiter unit
    i = math.pi/2
    e = 0
    
    K = 203 * (Pday)**(-1/3) * \
        ((MpMJ * math.sin(i)) / (MsMS + 9.458e-4 * (MpMJ)**(2/3))) * \
        (1 / (1 - e**2)**(1/2))

    return K


