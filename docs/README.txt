legacy_coord_search takes 5 arguments:
- catalog: An astropy SkyCoord object containing the target coordinates
- output: The name of the output file FIT catalog, must end in .FIT, .fits, etc...
- tile: Name of the tile MOC to use
- extent: Name of the extent MOC to use
- catalog_location: the path to the unzipped Legacy Peak Catalog.

The code will use the MOC files provided to determine whether or not the given
coordinates were observed, whether they lie on 5-sigma emmission, and the distance to,
and properties of the nearest source of peaked emmission (Designated TARG_).

The output consists of a FIT catalog with the following column names:
- RA: The user supplied RA
- DEC: The user supplied Declination
- IN_TILE: (True/False) Are the coordinates observed by the JCMT
- IN_EXT: (True/False) Are the coordinates within the extent of 5-sigma emmission
- TILEN: HEALpix tile number (order 6) of the coordinates (supplied even if target not observed by JCMT)
- TARG_ID: JCMT ID of the nearest neighbor source
- TARG_EXT: Parent extent of the nearest neighbor source
- TARG_RA: RA of the nearest neighbor source
- TARG_DEC: Declination of the nearest neighbor source
- TARG_SEP: Separation (in arcseconds) between the supplied coordinates and the nearest neighbor
- TARG_PEAKFLUX: Peak flux density (mJy/arcsec.sq) of the nearest neighbor source
- TARG_TILEN: Healpix tile number (order 6) of the nearest neigbor source
