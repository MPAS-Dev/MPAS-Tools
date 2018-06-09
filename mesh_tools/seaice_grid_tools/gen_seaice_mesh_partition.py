#!/usr/bin/env python

# tool for creating better than naive grid partitions for MPAS-Seaice

# requires python/anaconda-2.7 and that the MPAS-Tools mesh_conversion_tools are installed

from netCDF4 import Dataset
import os, math, string, sys
import numpy as np
import argparse
import subprocess

#-------------------------------------------------------------------

def degree_to_radian(degree):

    return (degree * math.pi) / 180.0

#-------------------------------------------------------------------

def add_cell_cull_array(filename, cullCell):

    mesh = Dataset(filename,"a")

    cullCellVariable = mesh.createVariable("cullCell","i4",("nCells"))

    cullCellVariable[:] = cullCell

    mesh.close()

#-------------------------------------------------------------------

def cull_mesh(meshToolsDir, filenameIn, filenameOut, cullCell):

    add_cell_cull_array(filenameIn, cullCell)

    os.system("%s/MpasCellCuller.x %s %s -c" %(meshToolsDir,filenameIn,filenameOut))

#-------------------------------------------------------------------

def load_partition(graphFilename):

    graphFile = open(graphFilename,"r")
    lines = graphFile.readlines()
    graphFile.close()

    partition = []
    for line in lines:
        partition.append(string.atoi(line))

    return partition

#-------------------------------------------------------------------

def get_cell_ids(culledFilename, originalFilename):

    culledFile = Dataset(culledFilename,"r")
    nCellsCulled = len(culledFile.dimensions["nCells"])

    originalFile = Dataset(originalFilename,"r")
    nCellsOriginal = len(originalFile.dimensions["nCells"])

    cellid = np.zeros(nCellsCulled, dtype=np.int)

    cellMapFile = open("cellMapForward.txt","r")
    cellMapLines = cellMapFile.readlines()
    cellMapFile.close()

    iCellOriginal = 0
    for cellMapLine in cellMapLines:

        if (iCellOriginal % 1000 == 0):
            print iCellOriginal, " of ", nCellsOriginal

        cellMap = string.atoi(cellMapLine)

        if (cellMap != -1):

            cellid[cellMap] = iCellOriginal

        iCellOriginal = iCellOriginal + 1

    return cellid

#-------------------------------------------------------------------

def get_cell_ids_orig(culledFilename, originalFilename):

    culledFile = Dataset(culledFilename,"r")
    latCellCulled = culledFile.variables["latCell"][:]
    lonCellCulled = culledFile.variables["lonCell"][:]
    nCellsCulled = len(culledFile.dimensions["nCells"])

    originalFile = Dataset(originalFilename,"r")
    latCellOriginal = originalFile.variables["latCell"][:]
    lonCellOriginal = originalFile.variables["lonCell"][:]
    nCellsOriginal = len(originalFile.dimensions["nCells"])

    cellid = np.zeros(nCellsCulled, dtype=np.int)

    for iCellCulled in range(0,nCellsCulled):

        if (iCellCulled % 1000 == 0):
            print "iCellCulled: ", iCellCulled, "of ", nCellsCulled

        for iCellOriginal in range(0,nCellsOriginal):

            if (latCellCulled[iCellCulled] == latCellOriginal[iCellOriginal] and \
                lonCellCulled[iCellCulled] == lonCellOriginal[iCellOriginal]):

                cellid[iCellCulled] = iCellOriginal
                break

    return cellid

#-------------------------------------------------------------------

# parsing
parser = argparse.ArgumentParser(description='Create sea ice grid partition')

parser.add_argument('-m', '--mesh',        dest="meshFilename",         required=True,  help='MPAS mesh file')
parser.add_argument('-r', '--regions',     dest="regionFilename",       required=True,  help='region file')
parser.add_argument('-n', '--nprocs',      dest="nProcs",               required=True,  help='number of processors', type=int)
parser.add_argument('-c', '--culler',      dest="mpasCullerLocation",   required=False, help='location of cell culler')
parser.add_argument('-o', '--outprefix',   dest="outputPrefix",         required=False, help='output graph file prefic', default="graph.info")
parser.add_argument('-p', '--plotting',    dest="plotting",             required=False, help='create diagnostic plotting file of partitions', action='store_true')
parser.add_argument('-g', '--metis',       dest="metis",                required=False, help='name of metis utility', default="gpmetis")
parser.add_argument('-e', '--equatorcull', dest="cullEquatorialRegion", required=False, help='create diagnostic plotting file of partitions', action='store_true')

args = parser.parse_args()

# required arguments
meshFilename = args.meshFilename
regionFilename = args.regionFilename
nProcsArray = [args.nProcs]

# optional arguments
if (args.mpasCullerLocation == None):
    meshToolsDir = os.path.dirname(os.path.realpath(__file__)) + "/../mesh_conversion_tools/"
else:
    meshToolsDir = args.mpasCullerLocation
decompName = args.outputPrefix
plotting = args.plotting
cullEquatorialRegion = args.cullEquatorialRegion
metis = args.metis

#./gen_seaice_mesh_partition.py -m /Users/l235697/Work/MPAS-CICE/Standalone_sim_configs/domains/domain_QU120km/seaice_QU_120km.nc -r /Users/l235697/Work/MPAS-CICE/Performance/sea_ice_coverage_partition/regions.nc -c /Users/l235697/Work/MPAS-CICE/MPAS-Tools/grid_gen/mesh_conversion_tools/ -n 32

#./gen_seaice_mesh_partition.py --mesh /Users/l235697/Work/MPAS-CICE/Standalone_sim_configs/domains/domain_QU120km/seaice_QU_120km.nc --regions region.nc --culler /Users/l235697/Work/MPAS-Tools/Feature_branches/generalize_partition_script/MPAS-Tools/grid_gen/mesh_conversion_tools --nprocs 32 --plotting --outprefix wibble2

plotFilename = "partition_diag.nc"

# get regions
regionFile = Dataset(regionFilename,"r")
nRegions = regionFile.nRegions
region = regionFile.variables["region"][:]
regionFile.close()

# diagnostics
if (plotting):
    os.system("cp %s %s" %(meshFilename,plotFilename))

# load mesh file
mesh = Dataset(meshFilename,"r")
nCells = len(mesh.dimensions["nCells"])
latCell = mesh.variables["latCell"][:]
mesh.close()

for nProcs in nProcsArray:

    if (cullEquatorialRegion):
        nBlocks = nRegions * nProcs
    else:
        nBlocks = nProcs

    combinedGraph = np.zeros(nCells)

    for iRegion in range(0,nRegions):

        # tmp file basename
        tmp = "%s_%2.2i_tmp" %(meshFilename,iRegion)

        # create precull file
        tmpFilenamesPrecull = tmp+"_precull.nc"
        os.system("cp %s %s" %(meshFilename,tmpFilenamesPrecull))

        # make cullCell variable
        cullCell = np.ones(nCells)

        for iCell in range(0,nCells):
            if (region[iCell] == iRegion):
                cullCell[iCell] = 0

        #for iCell in range(0,nCells):
        #    if (latCell[iCell] >= degree_to_radian(minLatitudeLimits[iRegion]) and \
        #        latCell[iCell] <  degree_to_radian(maxLatitudeLimits[iRegion])):
        #        cullCell[iCell] = 0

        # cull the mesh
        tmpFilenamesPostcull = tmp+"_postcull.nc"
        cull_mesh(meshToolsDir, tmpFilenamesPrecull, tmpFilenamesPostcull, cullCell)

        # preserve the initial graph file
        os.system("mv culled_graph.info culled_graph_%i_tmp.info" %(iRegion))

        # partition the culled grid
        try:
            graphFilename = "culled_graph_%i_tmp.info" %(iRegion)
            subprocess.call([metis, graphFilename, str(nProcs)])
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                print "metis program %s not found" %(metis)
                sys.exit()
            else:
                print "metis error"
                raise

        # get the cell IDs for this partition
        cellid = get_cell_ids(tmpFilenamesPostcull, meshFilename)

        # load this partition
        graph = load_partition("culled_graph_%i_tmp.info.part.%i" %(iRegion,nProcs))

        # add this partition to the combined partition
        for iCellPartition in range(0,len(graph)):
            if (cullEquatorialRegion):
                combinedGraph[cellid[iCellPartition]] = graph[iCellPartition] + nProcs * iRegion
            else:
                combinedGraph[cellid[iCellPartition]] = graph[iCellPartition]

    # output the cell partition file
    cellPartitionFile = open("%s.part.%i" %(decompName,nBlocks), "w")
    for iCell in range(0,nCells):
        cellPartitionFile.write("%i\n" %(combinedGraph[iCell]))
    cellPartitionFile.close()

    # output block partition file
    if (cullEquatorialRegion):
        blockPartitionFile = open("%s.part.%i" %(decompName,nProcs), "w")
        for iRegion in range(0,nRegions):
            for iProc in range(0,nProcs):
                blockPartitionFile.write("%i\n" %(iProc))
        blockPartitionFile.close()


    # diagnostics
    if (plotting):
        plottingFile = Dataset(plotFilename,"a")
        partitionVariable = plottingFile.createVariable("partition_%i" %(nProcs),"i4",("nCells"))
        partitionVariable[:] = combinedGraph
        plottingFile.close()

    os.system("rm *tmp*")
