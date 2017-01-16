# ----------------------------------------------------------------------------
#
# TITLE - LRCoords
# AUTHOR - James Lane
# PROJECT - jcmtlrsearch
# CONTENTS:
#   1. legacy_peakcat_search
#	2. _get_moc_info_
#   3. _get_peakcat_info_
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
''''''
__author__ = "James Lane"


# Imports
import numpy as np
from healpy.pixelfunc import ang2pix
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from pymoc.moc import MOC

# Functions
def legacy_peakcat_search(catalog, output, tile='jcmts850um_lr1_tile.fits',
                        extent='jcmts850um_lr1_extent.fits',
                        catalog_location='./'):
    '''
    legacy_coord_search:

    Take a catalog of points and cross-reference with the JCMT Legacy
    Release Catalogs. Outputs FITS catalog of information for each pair
    of coordinates.

    Args:
        coords (Astropy SkyCoord) - an astropy coordinate object containing
                                    the RA and Dec of the target positions.
        tile_moc (string) - The name of the Tile MOC, default is JCMT LR1
        extent_moc (string) - The name of the Extent MOC, default is JCMT LR1
        catalog_location (string) - Location of the Legacy Catalog relative to
                                    the current working directory, default is
                                    the current working directory.

    Returns:
        None
    '''

    # Get information from the MOCs.
    cells, contains = _get_moc_info_(catalog)

    # Get information from the Legacy Peak Catalogs.
    ids, coords, peaks, tiles = _get_peakcat_info_(catalog, cells[0], contains[0],
                                                   catalog_location)

    # Create the output columns for the table.
    cols = fits.ColDefs([
    fits.Column(name='RA', format='E', array=catalog.ra.degree, unit='deg'),
    fits.Column(name='DEC', format='E', array=catalog.ra.degree, unit='deg'),
    fits.Column(name='IN_TILE', format='L', array=contains[0]),
    fits.Column(name='IN_EXT', format='L', array=contains[1]),
    fits.Column(name='TILEN', format='J', array=cells[0]),
    fits.Column(name='TARG_ID', format='23A', array=ids[0]),
    fits.Column(name='TARG_EXT', format='23A', array=ids[1]),
    fits.Column(name='TARG_RA', format='E', array=coords[0], unit='deg'),
    fits.Column(name='TARG_DEC', format='E', array=coords[1], unit='deg'),
    fits.Column(name='TARG_SEP', format='E', array=coords[2], unit='arcsec'),
    fits.Column(name='TARG_PEAKFLUX', format='E', array=peaks, unit='mJy/arcsec**2'),
    fits.Column(name='TARG_TILEN', format='J', array=tiles)
    ])

    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.writeto(output)
#def

def _get_moc_info_(coords, tile='jcmts850um_lr1_tile.fits', extent='jcmts850um_lr1_extent.fits'):

    '''
    get_moc_info:

    Get information about the user supplied coordinates using two MOCs: a tile
    MOC that determines whether or not the coordinates were observed and an
    extent MOC that determines whether or not the coordinates are within a
    region of 5-sigma emmission.

    Args:
        coords (Astropy SkyCoord) - an astropy coordinate object containing
                                    the RA and Dec of the target positions.
        tile (string) - The name of the Tile MOC, default is JCMT LR1
        extent (string) - The name of the Extent MOC, default is JCMT LR1

    Returns:
    '''

    # Definitions suited for JCMT-LR1 *** Subject to parameterization ***
    tile_order = 6
    tile_nside = 2 ** tile_order
    pix_order = 16
    pix_nside = 2 ** pix_order

    # Load the MOCs
    tile_moc = MOC(filename=tile)
    ext_moc = MOC(filename=extent)

    # Declare coordinates in radians for performance
    ncoords = len(coords)
    coords_phi = coords.ra.radian
    coords_theta = (np.pi/2) - coords.dec.radian

    # Convert the coordinates into HEALpix cells for both Tile and Extent orders
    # Be aware that the tile order is used for identifying legacy release tiles,
    # the tile MOC still uses the pixel order!
    tile_cells = ang2pix(   nside=tile_nside,
                            theta=coords_theta,
                            phi=coords_phi,
                            nest=True)

    # Be aware that both tile and extent MOCs use the pixel order!
    pix_cells = ang2pix(  nside=pix_nside,
                            theta=coords_theta,
                            phi=coords_phi,
                            nest=True)

    # Check if the MOCs contain the HEALpix cells
    in_moc = np.empty((2,ncoords), dtype=bool)
    for i in range(ncoords):
        in_moc[0,i] = tile_moc.contains(order=pix_order, cell=pix_cells[i])
        in_moc[1,i] = ext_moc.contains(order=pix_order, cell=pix_cells[i])

    return  np.array([tile_cells,pix_cells], dtype=np.int64), in_moc
#def

def _get_peakcat_info_(coords, tile_cells, in_moc, cat_location='./'):
    '''
    get_peakcat_info:

    Use the Legacy Release peak catalogs to get information about peaked
    emmission near a series of user supplied coordinates. Requires get_moc_info
    output to provide catalog selection information.

    Args:
        coords (Astropy SkyCoord) - an astropy coordinate object containing
                                    the RA and Dec of the target positions.
        tile_cells (array) - HEALpix cell numbers corresponding to the Legacy
                             Release tile numbers. One HEALpix cell for each
                             coordinate position.
        in_moc (array) - Boolean array corresponding to whether or not the
                         position is a part of the Legacy Release.
        cat_location (string) - Location of the Legacy Catalog relative to the
                                current working directory, default is the
                                current working directory.

    Returns:

    '''

    # Defintions suited for JCMT LR1 *** Subject to parameterization ***

    # Find the number of unique tiles that need to be opened. Only choose tiles
    # that were observed in the survey.
    unique_tiles, unique_inv = np.unique(tile_cells[in_moc], return_inverse=True)

    # Loop over the unique tiles and read all information into single arrays.
    cat_id = np.array([],dtype='S23')
    cat_ra = np.array([])
    cat_dec = np.array([])
    cat_peak = np.array([])
    cat_ext = np.array([],dtype='S23')
    cat_tile = np.array([],dtype=int)

    for tile in unique_tiles:

        # Open the corresponding catalog
        cat_num = str(tile).zfill(6)
        cat_name = cat_location+'jcmts850um_peak-cat'+cat_num+'_pub_000.fits'
        cat_data = fits.open(cat_name)[1].data

        # Fill arrays
        cat_id = np.append(cat_id,cat_data['ID'])
        cat_ra = np.append(cat_ra,cat_data['RA'])
        cat_dec = np.append(cat_dec,cat_data['DEC'])
        cat_peak = np.append(cat_peak,cat_data['PEAK_FLUX'])
        cat_ext = np.append(cat_ext,cat_data['PARENT_EXTENT'])
        cat_tile = np.append(cat_tile,np.full(len(cat_data),tile,dtype=int))

    # Create a catalog object for cross-referencing
    cat_coords = SkyCoord(ra=cat_ra, dec=cat_dec, unit=(u.deg, u.deg))

    # Match user-supplied coordinates to the catalog to find nearest neighbors.
    idx, sep, dist = coords.match_to_catalog_sky(cat_coords)

    targ_coords = np.array([cat_ra[idx],cat_dec[idx],sep.arcsecond])
    targ_ids = np.array([cat_id[idx],cat_ext[idx]])
    targ_peak = cat_peak[idx]
    targ_tile = cat_tile[idx]

    return targ_ids, targ_coords, targ_peak, targ_tile
#def
