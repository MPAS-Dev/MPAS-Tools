#!/usr/bin/env python
# Phillip J. Wolfram
# 12/22/2015
"""
multires_scaled_hex.py 

Takes an existing MPAS mesh and coarse or refines the mesh outside a radial
region.  Mesh input should be a periodic_hex mesh.  Topology of mesh is
retained, which ensures that there are no non-hexagons within the mesh.

Example usage is 

    python multires_scaled_hex.py -i 250m_mesh.nc -o scaled_250m_mesh.nc

Phillip J. Wolfram
12/22/2015
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree as KDTree
import netCDF4
import shutil


def multires_scaled_hex(infname, outfname, xc=25000/2.0, yc=50000/2.0, radius=5000., alpha=1.05, ntimes=15, nlayers=0, plot=False):
    """
    Scales a hexagonal mesh from a radial region (xc,yc) with radius via an alpha scaling factor over ntimes
    layers of cells.  Stretching can be alternated via nlayers (experimental), which may prove useful via
    CVT smoothing codes.

    Phillip J. Wolfram
    12/22/2015
    """
    shutil.copyfile(infname, outfname)
    ds = netCDF4.Dataset(outfname,'r+')

    tree = KDTree(zip(ds.variables['xCell'], ds.variables['yCell']))
    _, center = tree.query([xc,yc])
    centerall = tree.query_ball_point([xc,yc], radius)

    xcenter = ds.variables['xCell'][center]
    ycenter = ds.variables['yCell'][center]
    x = ds.variables['xVertex'] - xcenter
    y = ds.variables['yVertex'] - ycenter
    xc = ds.variables['xCell'] - xcenter
    yc = ds.variables['yCell'] - ycenter
    

    # form rings aound the center
    dx = np.median(ds.variables['dvEdge'])
    dv = np.median(ds.variables['dcEdge'])/2.0
    cells = np.array(centerall).tolist()
    vertices = []
    for acell in centerall: 
        vertices.append(ds.variables['verticesOnCell'][acell,:ds.variables['nEdgesOnCell'][acell]]-1)
    vertices = np.unique(vertices).tolist()

    dalpha = alpha - 1.0
    alpha = 1.0
    for i in 1+np.arange(ntimes):
        for acell in cells[:]:
            for cellneighs in ds.variables['cellsOnCell'][acell]-1:
                cells.append(cellneighs)
                for avertex in ds.variables['verticesOnCell'][cellneighs,:ds.variables['nEdgesOnCell'][cellneighs]]-1:
                    vertices.append(avertex)
        cells = np.unique(cells).tolist()
        vertices = np.unique(vertices).tolist()
        # now have list of vertices and cells to NOT scale
    
        # alpha factor scaling
        if i <= ntimes/2.0:
            alpha += dalpha/ntimes
        else:
            alpha -= dalpha/ntimes

        # number of layers to scale
        if nlayers == 0 or not np.mod(i,nlayers):
            x *= alpha
            y *= alpha
            xc *= alpha
            yc *= alpha

            x[vertices] /= alpha
            y[vertices] /= alpha
            xc[cells] /= alpha
            yc[cells] /= alpha


        # plot incremental changes
        if plot:
            #plt.plot(x,y,'b.')
            plt.plot(xc,yc,'bo')
            plt.plot(0.0,0.0,'rx')
            plt.axis('equal')
            plt.show()

    # now write output
    ds.variables['xCell'][:] = xc + xcenter
    ds.variables['yCell'][:]= yc + ycenter
    ds.variables['xVertex'][:] = x + xcenter
    ds.variables['yVertex'][:] = y + ycenter

    ds.close()


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-i", "--infile", dest="inputfilename", help="grid files to be opened of form 'input*.nc'", metavar="FILE")
    parser.add_option("-o", "--outfile", dest="outputfilename", help="grid files to be output of form 'output*.nc'", metavar="FILE")
    parser.add_option("-x", "--xc", dest="xc", help="X location of centered region", metavar="FLOAT")
    parser.add_option("-y", "--yc", dest="yc", help="Y location of centered region", metavar="FLOAT")
    parser.add_option("-r", "--radius", dest="radius", help="Radius of centered region", metavar="FLOAT")
    parser.add_option("-a", "--alpha", dest="alpha", help="Scale factor for multiresolution", metavar="FLOAT")
    parser.add_option("-n", "--ntimes", dest="ntimes", help="Number of cell layers to scale for multiresolution", metavar="INT")


    options, args = parser.parse_args()
    if not options.inputfilename:
        parser.error("Input filename or expression ('-i') is a required input... e.g., -f 'input*.npz'")
    if not options.outputfilename:
        parser.error("Output filename or expression ('-o') is a required input... e.g., -f 'output*.npz'")
    if not options.xc:
        parser.error("X location of centered region is a required input... e.g., -x '12500'")
    else:
        options.xc = float(options.xc)
    if not options.yc:
        parser.error("Y location of centered region is a required input... e.g., -y '25000'")
    else:
        options.yc = float(options.yc)
    if not options.radius:
        parser.error("Radius of centered region is a required input... e.g., -r '5000.'")
    else:
        options.radius = float(options.radius)
    if not options.alpha:
        parser.error("Scale factor ... e.g., -a '1.05'")
    else:
        options.alpha = float(options.alpha)
    if not options.ntimes:
        parser.error("Number of cell layers to scale ... e.g., -n '15'")
    else:
        options.ntimes = int(options.ntimes)


    multires_scaled_hex(options.inputfilename, options.outputfilename, options.xc, options.yc, options.radius, \
            options.alpha, options.ntimes)
