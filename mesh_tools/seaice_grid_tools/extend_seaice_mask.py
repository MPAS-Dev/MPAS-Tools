#!/usr/bin/env python
import argparse

from mpas_tools.seaice.mask import extend_seaice_mask


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Extend the ice presence variable')

    parser.add_argument('-m', '--inputmesh',  dest="filenameMesh",     required=True,  help='MPAS mesh file for source regridding mesh')
    parser.add_argument('-p', '--presence',   dest="filenamePresence", required=True,  help='File with ice presence')
    parser.add_argument('-e', '--extenddist', dest="extendDistance",   required=True,  help='distance (km) to extend ice present region', type=float)
    parser.add_argument('-u', '--unitsphere', dest="unitSphere",       required=False, help='Is the mesh file a unit sphere', action='store_false')

    args = parser.parse_args()

    extend_seaice_mask(args.filenameMesh,args.filenamePresence,args.extendDistance,args.unitSphere)
