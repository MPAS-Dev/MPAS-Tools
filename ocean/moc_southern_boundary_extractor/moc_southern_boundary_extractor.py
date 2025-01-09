#!/usr/bin/env python

'''
This script takes a mesh file (-m flag) and a file with MOC regions masks
(-f flag) produce by the MPAS mask creator.  The script produces a copy of
the contents of the MOC mask file, adding transects that mark the southern
boundary of each region in a file indicated with the -o flag.  The transect
is applied only to vertices and edges, not cells, because the need for southern
boundary transect data on cells is not foreseen.
'''

from mpas_tools.ocean.moc import (
    moc_southern_boundary_extractor
)

if __name__ == '__main__':
    moc_southern_boundary_extractor()
