from netCDF4 import Dataset
import math
import argparse

parser = argparse.ArgumentParser(description='Set mesh radius')

parser.add_argument('-m', '--meshFilename', dest="meshFilename", required=True, help='MPAS mesh file for resizing')
parser.add_argument('-r', '--sphereRadius', dest="sphereRadius", help='New sphere radius (m)')

args = parser.parse_args()


    
filein = Dataset(args.meshFilename,"a")

oldSphereRadius = filein.sphere_radius

if (args.sphereRadius == None):
    sphereRadius = 6371229.0
else:
    sphereRadius = args.sphereRadius

rescale = sphereRadius / oldSphereRadius
rescale2 = math.pow(rescale,2)

filein.sphere_radius = sphereRadius

nCells = len(filein.dimensions["nCells"])

xCell = filein.variables["xCell"];
yCell = filein.variables["yCell"];
zCell = filein.variables["zCell"];
xCell[:] *= rescale
yCell[:] *= rescale
zCell[:] *= rescale

xEdge = filein.variables["xEdge"];
yEdge = filein.variables["yEdge"];
zEdge = filein.variables["zEdge"];
xEdge[:] *= rescale
yEdge[:] *= rescale
zEdge[:] *= rescale

xVertex = filein.variables["xVertex"];
yVertex = filein.variables["yVertex"];
zVertex = filein.variables["zVertex"];
xVertex[:] *= rescale
yVertex[:] *= rescale
zVertex[:] *= rescale

dcEdge = filein.variables["dcEdge"];
dvEdge = filein.variables["dvEdge"];
dcEdge[:] *= rescale
dvEdge[:] *= rescale

areaCell = filein.variables["areaCell"];
areaTriangle = filein.variables["areaTriangle"];
kiteAreasOnVertex = filein.variables["kiteAreasOnVertex"];
areaCell[:] *= rescale
areaTriangle[:] *= rescale
kiteAreasOnVertex[:] *= rescale

filein.close()
