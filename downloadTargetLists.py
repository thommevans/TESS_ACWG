"""
TODO: 

Add routines that are simple to call and download the needed information
from the NASA Exoplanet Archive. The goal is to avoid the cumbersome process
of visiting the NASA Exoplanet Archive, manually selecting the columns of 
interest, then downloading the csv files each time. All of this should be 
automated with a routine that is called once each time the csv target lists
need to be updated, i.e. at least once per month.

There need to be 2x routines:

1. targetsPublishedConfirmation()

Source = 'Confirmed planets' on NASA Exoplanet Archive.

These targets pass a stringent confirmation condition, which I believe is 
a peer-reviewed, published confirmation/validation paper. A confirmation
paper posted to arXiv is not sufficient, until it has been published in
a peer-reviewed journal.

2. targetsUnpublishedTOIs()

Source = 'TESS Project Candidates'

These targets have been identified as TESS 'Objects of Interest' (TOIs). 
However, they retain a 'candidate' status until a confirmation paper has
been published in a peer-reviewed journal. As such, NASA Exoplanet Archive
will not list a mass for these candidates, because if a published mass was
available, a given target would no longer be listed as a candidate. 
Therefore, an empirical mass-radius relation is used to estimate what the 
mass of each candidate is likely to be, allowing subsequent calculations 
of transmission spectroscopy metrics to be performed (with the big caveat
that empirical .


"""

def targetsWithPublishedConfirmation():
    """
    Confirmed Planets from NASA Exoplanet Archive.    
    1. Remove condition 'Default parameter set = 1'.
    2. 'Select columns' > 'Clear all'
    3. Then select these columns:
        1. Planet name
        2. Default parameter set
        3. Discovery facility
        4. Detected by transits
        5. Orbital period (day)
        6. Orbit semi-major axis (AU)
        7. Radius value (Earth radius)
        8. Radius upper uncertainty (Earth radius)
        9. Radius lower uncertainty (Earth radius)
       10. Mass value (Earth mass)
       11. Mass upper uncertainty (Earth mass)
       12. Mass lower uncertainty (Earth mass)
       13. Mass/Mass*sin(i) Provenance
       14. Eccentricity
       15. Insolation flux (Earth flux)
       16. Inclination (deg)
       17. Impact parameter
       18. Transit duration (hour)
       19. Ratio of semi-major axis to stellar radius
       20. Ratio of planet to stellar radius
       21. Stellar effective temperature (Kelvin)
       22. Stellar radius (solar radius)
       23. Stellar mass (solar mass)
       24. Stellar log10 gravity (CGS units, i.e. cm^2/s)
       25. RA (deg)
       26. Dec (dec)
       27. Distance (parsec)
       28. V (mag)
       29. J (mag)
       30. H (mag)
       31. Ks (mag)
       Click 'Update'. Close 'Select columns' popup.
    4. Set condition 'Detected by transits = 1'.
    5. Then 'Download Table' > 'CSV format'.
    """
    confirmedFpath = None
    return confirmedFpath

def targetsUnpublishedTOIs():
    """
    TESS Project Candidates from NASA Exoplanet Archive.
    1. 'Select columns' > 'Clear all'
    2. Then select these columns:
        1. TESS Object of Interest
        2. TFOPWG Disposition
        3. RA (deg) [ click '+' next to 'RA (sexagesimal)' ]
        4. Dec (dec) [ click '+' next to 'Dec (sexagesimal)' ]
        5. Orbital Period (day)
        6. Transit duration (hour)
        7. Radius value (Earth radii)
        8. Radius upper uncertainty (Earth radius)
        9. Radius lower uncertainty (Earth radius)
       10. Insolation flux (Earth flux)
       11. Equilibrium Temperature (Kelvin)
       12. TESS magnitude
       13. Stellar effective temperature (Kelvin)
       14. Stellar log10 gravity (CGS units, i.e. cm^2/s)
       15. Stellar radius (solar radius)
       Click 'Update'. Close 'Select columns' popup.
    3. Then 'Download Table' > 'CSV format'.
    """
    toiFpath = None
    return toiFpath
