from netCDF4 import Dataset
import os
import math
import errno
import numpy as np
import subprocess
import argparse
import shutil

from .regrid import regrid_to_other_mesh
from .mask import extend_seaice_mask
from .regions import make_regions_file


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
    executable = "MpasCellCuller.x"
    if meshToolsDir is not None:
        executable = os.path.join(meshToolsDir, executable)

    subprocess.run([executable, filenameIn, filenameOut, "-c"], check=True)

#-------------------------------------------------------------------

def load_partition(graphFilename):

    graphFile = open(graphFilename,"r")
    lines = graphFile.readlines()
    graphFile.close()

    partition = []
    for line in lines:
        partition.append(int(line))

    return partition

#-------------------------------------------------------------------

def get_cell_ids(culledFilename, originalFilename):

    culledFile = Dataset(culledFilename,"r")
    nCellsCulled = len(culledFile.dimensions["nCells"])

    originalFile = Dataset(originalFilename,"r")
    nCellsOriginal = len(originalFile.dimensions["nCells"])

    cellid = np.zeros(nCellsCulled, dtype=int)

    cellMapFile = open("cellMapForward.txt","r")
    cellMapLines = cellMapFile.readlines()
    cellMapFile.close()

    iCellOriginal = 0
    for cellMapLine in cellMapLines:

        if (iCellOriginal % 1000 == 0):
            print(iCellOriginal, " of ", nCellsOriginal)

        cellMap = int(cellMapLine)

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

    cellid = np.zeros(nCellsCulled, dtype=int)

    for iCellCulled in range(0,nCellsCulled):

        if (iCellCulled % 1000 == 0):
            print("iCellCulled: ", iCellCulled, "of ", nCellsCulled)

        for iCellOriginal in range(0,nCellsOriginal):

            if (latCellCulled[iCellCulled] == latCellOriginal[iCellOriginal] and
                lonCellCulled[iCellCulled] == lonCellOriginal[iCellOriginal]):

                cellid[iCellCulled] = iCellOriginal
                break

    return cellid

#-------------------------------------------------------------------

def gen_seaice_mesh_partition(meshFilename, regionFilename, nProcs, mpasCullerLocation, outputPrefix, plotting, metis, cullEquatorialRegion):

    # arguments
    meshToolsDir = mpasCullerLocation
    if meshToolsDir is None:
        culler = shutil.which("MpasCellCuller.x")
        if culler is not None:
            meshToolsDir = os.path.dirname(culler)
        else:
            # no directory was provided and none
            this_dir = os.path.dirname(os.path.realpath(__file__))
            meshToolsDir = os.path.abspath(os.path.join(
                this_dir, "..", "..", "..", "mesh_tools",
                "mesh_conversion_tools"))
            culler = os.path.join(meshToolsDir, "MpasCellCuller.x")
        if not os.path.exists(culler):
            raise FileNotFoundError(
                "MpasCellCuller.x does not exist at the requested location.")

    plotFilename = "partition_diag.nc"

    # get regions
    regionFile = Dataset(regionFilename,"r")
    nRegions = regionFile.nRegions
    region = regionFile.variables["region"][:]
    regionFile.close()

    # diagnostics
    if plotting:
        shutil.copyfile(meshFilename, plotFilename)

    # load mesh file
    mesh = Dataset(meshFilename,"r")
    nCells = len(mesh.dimensions["nCells"])
    latCell = mesh.variables["latCell"][:]
    mesh.close()

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
        shutil.copyfile(meshFilename, tmpFilenamesPrecull)

        # make cullCell variable
        cullCell = np.ones(nCells)

        for iCell in range(0,nCells):
            if (region[iCell] == iRegion):
                cullCell[iCell] = 0

        # cull the mesh
        tmpFilenamesPostcull = tmp+"_postcull.nc"
        cull_mesh(meshToolsDir, tmpFilenamesPrecull, tmpFilenamesPostcull, cullCell)

        # preserve the initial graph file
        os.rename("culled_graph.info", f"culled_graph_{iRegion}_tmp.info")

        # partition the culled grid
        try:
            graphFilename = "culled_graph_%i_tmp.info" %(iRegion)
            subprocess.call([metis, graphFilename, str(nProcs)])
        except OSError as e:
            if e.errno == errno.ENOENT:
                raise FileNotFoundError("metis program %s not found" %(metis))
            else:
                print("metis error")
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
    cellPartitionFile = open("%s.part.%i" %(outputPrefix,nBlocks), "w")
    for iCell in range(0,nCells):
        cellPartitionFile.write("%i\n" %(combinedGraph[iCell]))
    cellPartitionFile.close()

    # output block partition file
    if (cullEquatorialRegion):
        blockPartitionFile = open("%s.part.%i" %(outputPrefix,nProcs), "w")
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

    subprocess.run("rm *tmp*", shell=True)


def prepare_partitions():
    # parsing input
    parser = argparse.ArgumentParser(description='Perform preparatory work for making seaice partitions.')

    parser.add_argument('-i', '--inputmesh', dest="meshFilenameSrc", required=True,
                        help='MPAS mesh file for source regridding mesh.')
    parser.add_argument('-p', '--presence', dest="filenameData", required=True,
                        help='Input ice presence file for source mesh.')
    parser.add_argument('-m', '--outputmesh', dest="meshFilenameDst", required=True,
                        help='MPAS mesh file for destination regridding mesh.')
    parser.add_argument('-o', '--outputDir', dest="outputDir", required=True,
                        help='Output directory for temporary files and partition files.')

    args = parser.parse_args()

    # Check if output directory exists
    if (not os.path.isdir(args.outputDir)):
        raise FileNotFoundError("ERROR: Output directory does not exist.")

    # 1) Regrid the ice presence from the input data mesh to the grid of choice
    print("Regrid to desired mesh...")
    filenameOut = args.outputDir + "/icePresent_regrid.nc"
    regrid_to_other_mesh(args.meshFilenameSrc, args.filenameData, args.meshFilenameDst, filenameOut)

    # 2) create icePresence variable
    print("fix_regrid_output...")

    # check executable exists
    if shutil.which("fix_regrid_output.exe") is not None:
        # it's in the system path
        executable = "fix_regrid_output.exe"
    elif os.path.exists("./fix_regrid_output.exe"):
        # found in local path
        executable = "./fix_regrid_output.exe"
    else:
        raise FileNotFoundError("fix_regrid_output.exe could not be found.")

    inputFile = args.outputDir + "/icePresent_regrid.nc"
    outputFile = args.outputDir + "/icePresent_regrid_modify.nc"
    subprocess.call([executable, inputFile, args.meshFilenameDst, outputFile])

    # 3) create variable icePresenceExtended
    print("extend_seaice_mask...")
    filenamePresence = args.outputDir + "/icePresent_regrid_modify.nc"
    extend_seaice_mask(args.meshFilenameDst, filenamePresence, 0.0, False)

    # 4) Make the regions file from the icePresenceExtended variable
    print("make_regions_file...")
    filenameIcePresent = args.outputDir + "/icePresent_regrid_modify.nc"
    filenameOut = args.outputDir + "/regions.nc"
    make_regions_file(filenameIcePresent, args.meshFilenameDst, "two_region_eq", "icePresenceExtended", 0.5,
                      filenameOut)


def create_partitions():

    # parsing input
    parser = argparse.ArgumentParser(description='Create sea-ice partitions.')

    parser.add_argument('-m', '--outputmesh', dest="meshFilename", required=True,
                        help='MPAS mesh file for destination regridding mesh.')
    parser.add_argument('-o', '--outputDir', dest="outputDir", required=True,
                        help='Output directory for temporary files and partition files.')
    parser.add_argument('-c', '--cullerDir', dest="mpasCullerLocation", required=False,
                        help='Location of MPAS MpasCellCuller.x executable.')
    parser.add_argument('-p', '--prefix', dest="outputPrefix", required=False,
                        help='prefix for output partition filenames.', default="graph.info")
    parser.add_argument('-x', '--plotting', dest="plotting", required=False,
                        help='create diagnostic plotting file of partitions', action='store_false')
    parser.add_argument('-g', '--metis', dest="metis", required=False, help='name of metis utility', default="gpmetis")
    parser.add_argument('-n', '--nProcs', dest="nProcs", required=False,
                        help='number of processors to create partition for.', type=int)
    parser.add_argument('-f', '--nProcsFile', dest="nProcsFile", required=False,
                        help='number of processors to create partition for.')

    args = parser.parse_args()

    # number of processors
    nProcsArray = []
    if (args.nProcs is not None and args.nProcsFile is None):
        nProcsArray.append(args.nProcs)
    elif (args.nProcs is None and args.nProcsFile is not None):
        fileNProcs = open(args.nProcsFile, "r")
        nProcsLines = fileNProcs.readlines()
        fileNProcs.close()
        for line in nProcsLines:
            nProcsArray.append(int(line))
    elif (args.nProcs is None and args.nProcsFile is None):
        raise ValueError("Must specify nProcs or nProcsFile")
    elif (args.nProcs is not None and args.nProcsFile is not None):
        raise ValueError("Can't specify both nProcs or nProcsFile")

    # create partitions
    regionFilename = args.outputDir + "/regions.nc"
    outputPrefix = args.outputDir + "/" + args.outputPrefix

    for nProcs in nProcsArray:
        gen_seaice_mesh_partition(args.meshFilename, regionFilename, nProcs, args.mpasCullerLocation, outputPrefix,
                                  args.plotting, args.metis, False)
