# ----------------------------------------------------------------------------
#
# TITLE - LRSearch_test
# AUTHOR - James Lane
# PROJECT - jcmtlrsearch
#
# ----------------------------------------------------------------------------
#
# Docstrings and metadata:
'''Testing the LRSearch module'''
__author__ = "James Lane"

#Imports
from jcmtlrsearch import LRSearch
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

# Define the test coordinates to use: 4 points in Orion A South
# The first two sets of coordinates are in HH1/2 and L1641S respectively,
# the third set of coordinates is in Orion AS but not ontop of significant
# emmission, the fourth set of coordinates is about half a degree west of the
# observation boundary.
coord_ra = np.array(['05:36:23.042','05:40:48.733','05:41:16.235','05:33:23.836'])
coord_dec = np.array(['-06:46:10.04','-08:06:44.65','-08:15:56.26','-06:55:52.32'])

coords = SkyCoord(ra=coord_ra, dec=coord_dec, frame='icrs', unit=(u.hourangle, u.deg))

LRSearch.legacy_peakcat_search(coords, 'outcat.FIT', catalog_location='LR1_peak-cat/')
