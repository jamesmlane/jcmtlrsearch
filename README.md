# jcmtlrsearch

## Introduction

JCMTLRSearch is a module for searching through the JCMT Legacy Release Catalogs.

## Requirements

The module requires the following Python packages:
- numpy
- healpy
- pymoc
- astropy

Both extent and tile MOC files are required, they can be found [here](http://www.eaobservatory.org/jcmt/science/archive/lr1/):

A compressed version of the Legacy Release Peak Catalogs is
[included](catalogs/) for convenience.

## Using the Code

The currently consists of one function: legacy_peakcat_search

It can be used as follows:

```
from jcmtlrsearch import LRSearch
LRSearch.legacy_peakcat_search( *args )
```

This function will take in a catalog of coordinates and cross-reference with the
Legacy Release Peak Catalogs. The coordinates are required to be input as an
Astropy SkyCoord object, which allows for flexibility on the part of the user.

At this point there are some files that are required to be stored locally,
including the 2 MOC files and the Legacy Peak Catalog.

See the [documentation](docs/README.txt) for more information on arguments and
output.

Also see the test script for an example of how to use the code as well as the
output (Note you will need to download the MOC files in order to run the test
script).
