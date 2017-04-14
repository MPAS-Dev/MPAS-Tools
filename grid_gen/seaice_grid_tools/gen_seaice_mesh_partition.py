#!/usr/bin/env python

# tool for creating better than naive grid partitions for MPAS-Seaice

# requires python/anaconda-2.7
# The MPAS-Tools mesh_conversion_tools installed

from netCDF4 import Dataset
import os, math, string
import numpy as np
import argparse

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

    #meshToolsDir = "/Users/akt/Work/MPAS-CICE/MPAS-Tools/grid_gen/mesh_conversion_tools_v2/"
    
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

parser.add_argument('-m', dest="meshFilename",        required=True,  help='MPAS mesh file')
parser.add_argument('-c', dest="mpasCullerLocation",  required=True,  help='location of cell culler')
parser.add_argument('-n', dest="nProcs",              required=True,  help='number of processors', type=int)

args = parser.parse_args()

meshFilename = args.meshFilename
meshToolsDir = args.mpasCullerLocation
nProcsArray = [args.nProcs]

#./grid_partition.py -m seaice_QU_120km.nc -c /Users/akt/Work/MPAS-CICE/MPAS-Tools/grid_gen/mesh_conversion_tools_v2/ -n 32
#./grid_partition.py -m seaice_QU_120km.nc -c /Users/akt/Work/MPAS-CICE/MPAS-Tools_dev/feature_branches/mapping_output/MPAS-Tools/grid_gen/mesh_conversion_tools/ -n 32

# input name
#meshFilename = "/Users/akt/Work/MPAS-CICE/Standalone_sim_configs/domains/domain_RRS.18-6km/seaice.RRS18to6v3.170111.nc"
#meshFilename = "/Users/akt/Work/MPAS-CICE/Standalone_sim_configs/domains/domain_RRS.30-10km/seaice.RSS.30-10km.161221.nc"
#meshFilename = "seaice_QU_120km.nc"

#nProcsArray = [32]
#nProcsArray = [36000]

plotting  = True

cullEquatorialRegion = False
decompName = "lblat"

if (cullEquatorialRegion):
    minLatitudeLimits = [-100.0, -35.0,  35.0]
    maxLatitudeLimits = [ -35.0,  35.0, 100.0]
else:
    minLatitudeLimits = [-100.0, -60.0,  60.0]
    maxLatitudeLimits = [ -60.0,  60.0, 100.0]

nRegions = len(minLatitudeLimits)

# diagnostics
if (plotting):
    os.system("cp %s plotting.nc" %(meshFilename))

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
            if (latCell[iCell] >= degree_to_radian(minLatitudeLimits[iRegion]) and \
                latCell[iCell] <  degree_to_radian(maxLatitudeLimits[iRegion])):
                cullCell[iCell] = 0

        # cull the mesh
        tmpFilenamesPostcull = tmp+"_postcull.nc"
        cull_mesh(meshToolsDir, tmpFilenamesPrecull, tmpFilenamesPostcull, cullCell)

        # preserve the initial graph file
        os.system("mv culled_graph.info culled_graph_%i_tmp.info" %(iRegion))
    
        # partition the culled grid
        os.system("gpmetis culled_graph_%i_tmp.info %i" %(iRegion,nProcs))

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
    cellPartitionFile = open("graph.info.%s_%i.part.%i" %(decompName,nRegions,nBlocks), "w")
    for iCell in range(0,nCells):
        cellPartitionFile.write("%i\n" %(combinedGraph[iCell]))
    cellPartitionFile.close()

    # output block partition file
    if (cullEquatorialRegion):
        blockPartitionFile = open("block.info.%s_%i.part.%i" %(decompName,nRegions,nProcs), "w")
        for iRegion in range(0,nRegions):
            for iProc in range(0,nProcs):
                blockPartitionFile.write("%i\n" %(iProc))
        blockPartitionFile.close()

    
    # diagnostics
    if (plotting):
        plottingFile = Dataset("plotting.nc","a")
        partitionVariable = plottingFile.createVariable("partition_%i" %(nProcs),"i4",("nCells"))
        partitionVariable[:] = combinedGraph
        plottingFile.close()

    os.system("rm *tmp*")
