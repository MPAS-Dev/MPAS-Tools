#!/usr/bin/env python
import sys

from seaice.mesh import make_mpas_scripfile_on_cells

if __name__ == "__main__":

    if (len(sys.argv) != 4):
        print("Usage: seaice_mesh.py meshFilename, scripFilename, title")
        sys.exit()

    meshFilename  = sys.argv[1]
    scripFilename = sys.argv[2]
    title         = sys.argv[3]

    make_mpas_scripfile_on_cells(meshFilename, scripFilename, title)
