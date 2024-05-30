from netCDF4 import Dataset
import os
import math
import errno
import numpy as np
import subprocess
import argparse
import shutil
import glob

from .regrid import regrid_to_other_mesh
from .mask import extend_seaice_mask
from .regions import make_regions_file


def gen_seaice_mesh_partition(meshFilename, regionFilename, nProcsArray,
                              mpasCullerLocation, outputPrefix, plotting,
                              metis, cullEquatorialRegion):
    """
    Generate graph partition(s) for the given MPAS-Seaice mesh and the given
    number(s) of processors and a file defining regions that each processor
    should own part of (typically a polar region and an equatorial region)

    Parameters
    ----------
    meshFilename : str
        The name of a file containing the MPAS-Seaice mesh

    regionFilename : str
        The name of a file containing a ``region`` field defining different
        regions that each processor should own part of

    nProcsArray : list or int
        The number(s) of processors to create graph partitions for

    mpasCullerLocation : str or None
        The directory for the ``MpasCellCuller.x`` tool or ``None`` to look in
        the user's path

    outputPrefix : str
        The prefix to prepend to each graph partition file

    plotting : bool
        Whether to create a NetCDF file ``partition_diag.nc`` to use for
        plotting the partitions

    metis : str
        The exectable to use for partitioning in each region

    cullEquatorialRegion : bool
        Whether to remove the equatorial region from the paritions
    """

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
    regionFile = Dataset(regionFilename, "r")
    nRegions = regionFile.nRegions
    region = regionFile.variables["region"][:]
    regionFile.close()

    # diagnostics
    if plotting:
        shutil.copyfile(meshFilename, plotFilename)

    # load mesh file
    mesh = Dataset(meshFilename, "r")
    nCells = len(mesh.dimensions["nCells"])
    mesh.close()

    cellidsInRegion = []

    for iRegion in range(0, nRegions):

        # tmp file basename
        tmp = "%s_%2.2i_tmp" % (meshFilename, iRegion)

        # create precull file
        tmpFilenamesPrecull = tmp + "_precull.nc"
        shutil.copyfile(meshFilename, tmpFilenamesPrecull)

        # make cullCell variable
        cullCell = np.ones(nCells)

        for iCell in range(0, nCells):
            if region[iCell] == iRegion:
                cullCell[iCell] = 0

        # cull the mesh
        tmpFilenamesPostcull = tmp + "_postcull.nc"
        _cull_mesh(meshToolsDir, tmpFilenamesPrecull, tmpFilenamesPostcull,
                   cullCell)

        # get the cell IDs for this partition
        cellid = _get_cell_ids(tmpFilenamesPostcull, meshFilename)
        cellidsInRegion.append(cellid)

        # preserve the initial graph file
        os.rename("culled_graph.info", f"culled_graph_{iRegion}_tmp.info")

    if not isinstance(nProcsArray, (list, tuple, set)):
        # presumably, it's a single integer
        nProcsArray = [nProcsArray]

    for nProcs in nProcsArray:
        if cullEquatorialRegion:
            nBlocks = nRegions * nProcs
        else:
            nBlocks = nProcs

        combinedGraph = np.zeros(nCells)

        for iRegion in range(0, nRegions):

            # partition the culled grid
            try:
                graphFilename = "culled_graph_%i_tmp.info" % iRegion
                subprocess.call([metis, graphFilename, str(nProcs)])
            except OSError as e:
                if e.errno == errno.ENOENT:
                    raise FileNotFoundError(
                        "metis program %s not found" % metis)
                else:
                    print("metis error")
                    raise

            cellid = cellidsInRegion[iRegion]

            # load this partition
            graph = _load_partition(
                "culled_graph_%i_tmp.info.part.%i" %
                (iRegion, nProcs))

            # add this partition to the combined partition
            for iCellPartition in range(0, len(graph)):
                if cullEquatorialRegion:
                    combinedGraph[cellid[iCellPartition]] = \
                        graph[iCellPartition] + nProcs * iRegion
                else:
                    combinedGraph[cellid[iCellPartition]] = \
                        graph[iCellPartition]

        # output the cell partition file
        cellPartitionFile = open("%s.part.%i" % (outputPrefix, nBlocks), "w")
        for iCell in range(0, nCells):
            cellPartitionFile.write("%i\n" % (combinedGraph[iCell]))
        cellPartitionFile.close()

        # output block partition file
        if cullEquatorialRegion:
            blockPartitionFile = open(
                "%s.part.%i" %
                (outputPrefix, nProcs), "w")
            for iRegion in range(0, nRegions):
                for iProc in range(0, nProcs):
                    blockPartitionFile.write("%i\n" % iProc)
            blockPartitionFile.close()

        # diagnostics
        if plotting:
            plottingFile = Dataset(plotFilename, "a")
            partitionVariable = plottingFile.createVariable(
                "partition_%i" % nProcs, "i4", ("nCells",))
            partitionVariable[:] = combinedGraph
            plottingFile.close()

    subprocess.run("rm *tmp*", shell=True)


def prepare_partitions():
    """
    An entry point for performing preparatory work for making seaice partitions
    """
    # parsing input
    parser = argparse.ArgumentParser(
        description="Perform preparatory work for making seaice partitions.")

    parser.add_argument("-i", "--inputmesh", dest="meshFilenameSrc",
                        required=True,
                        help="MPAS mesh file for source regridding mesh.")
    parser.add_argument("-p", "--presence", dest="filenameData",
                        required=True,
                        help="Input ice presence file for source mesh.")
    parser.add_argument("-m", "--outputmesh", dest="meshFilenameDst",
                        required=True,
                        help="MPAS mesh file for destination regridding mesh.")
    parser.add_argument("-o", "--outputDir", dest="outputDir",
                        required=True,
                        help="Output directory for temporary files and "
                             "partition files.")
    parser.add_argument("-w", "--weightsFilename", dest="weightsFilename",
                        required=False,
                        help="A mapping file between the input and output "
                             "MPAS meshes.  One will be generated if it is "
                             "not supplied.")

    args = parser.parse_args()

    # Check if output directory exists
    if not os.path.isdir(args.outputDir):
        raise FileNotFoundError("ERROR: Output directory does not exist.")

    # 1) Regrid the ice presence from the input data mesh to the grid of choice
    print("Regrid to desired mesh...")
    filenameOut = args.outputDir + "/icePresent_regrid.nc"

    regrid_to_other_mesh(
        meshFilenameSrc=args.meshFilenameSrc,
        filenameData=args.filenameData,
        meshFilenameDst=args.meshFilenameDst,
        filenameOut=filenameOut,
        generateWeights=(args.weightsFilename is None),
        weightsFilename=args.weightsFilename)

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
    make_regions_file(filenameIcePresent=filenameIcePresent,
                      filenameMesh=args.meshFilenameDst,
                      regionType="two_region_eq",
                      varname="icePresenceExtended",
                      limit=0.5, filenameOut=filenameOut)


def create_partitions():
    """
    An entry point for creating sea-ice partitions
    """

    # parsing input
    parser = argparse.ArgumentParser(description='Create sea-ice partitions.')

    parser.add_argument(
        '-m', '--outputmesh', dest="meshFilename", required=True,
        help='MPAS mesh file for destination regridding mesh.')
    parser.add_argument(
        '-o', '--outputDir', dest="outputDir", required=True,
        help='Output directory for temporary files and partition files.')
    parser.add_argument(
        '-c', '--cullerDir', dest="mpasCullerLocation", required=False,
        help='Location of MPAS MpasCellCuller.x executable.')
    parser.add_argument(
        '-p', '--prefix', dest="outputPrefix", required=False,
        help='prefix for output partition filenames.', default="graph.info")
    parser.add_argument(
        '-x', '--plotting', dest="plotting", required=False,
        help='create diagnostic plotting file of partitions',
        action='store_true')
    parser.add_argument(
        '-g', '--metis', dest="metis", required=False,
        help='name of metis utility', default="gpmetis")
    parser.add_argument(
        '-n', '--nProcs', dest="nProcsArray", nargs='*', required=False,
        help='list of the number of processors to create partition for.',
        type=int)
    parser.add_argument(
        '-f', '--nProcsFile', dest="nProcsFile", required=False,
        help='number of processors to create partition for.')

    args = parser.parse_args()

    # number of processors
    if args.nProcsArray is None and args.nProcsFile is None:
        raise ValueError("Must specify nProcs or nProcsFile")
    if args.nProcsArray is not None and args.nProcsFile is not None:
        raise ValueError("Can't specify both nProcs or nProcsFile")

    if args.nProcsFile is not None:
        with open(args.nProcsFile, "r") as fileNProcs:
            nProcsLines = fileNProcs.readlines()
            nProcsArray = [int(line) for line in nProcsLines
                           if line.split() != '']
    else:
        nProcsArray = args.nProcsArray

    # create partitions
    regionFilename = args.outputDir + "/regions.nc"
    outputPrefix = args.outputDir + "/" + args.outputPrefix

    gen_seaice_mesh_partition(args.meshFilename, regionFilename, nProcsArray,
                              args.mpasCullerLocation, outputPrefix,
                              args.plotting, args.metis,
                              cullEquatorialRegion=False)


def simple_partitions():
    """
    An entry point for creating sea-ice partitions on LCRC (Anvil and
    Chrysalis)
    """

    data_dir = '/lcrc/group/e3sm/public_html/mpas_standalonedata/' \
               'mpas-seaice/partition'

    # parsing input
    parser = argparse.ArgumentParser(
        description='Create sea-ice partitions on LCRC.')

    parser.add_argument(
        '-m', '--mesh', dest="meshFilename", required=True,
        help='MPAS-Seaice mesh file.')
    parser.add_argument(
        '-p', '--prefix', dest="outputPrefix", required=True,
        help='prefix for output partition filenames.')
    parser.add_argument(
        '-n', '--nprocs', dest="nProcsArray", nargs='*', required=True,
        help='list of the number of processors to create partition for.',
        type=int)
    parser.add_argument(
        '-d', '--datadir', dest="dataDir", required=False,
        default=data_dir,
        help='Directory with seaice_QU60km_polar.nc and '
             'icePresent_QU60km_polar.nc.')

    args = parser.parse_args()

    meshFilenameDst = os.path.abspath(args.meshFilename)

    tmpdir = 'tmp_seaice_part_dir'
    try:
        shutil.rmtree(tmpdir)
    except FileNotFoundError:
        pass

    os.makedirs(tmpdir)

    cwd = os.getcwd()

    os.chdir(tmpdir)

    # make a local link to the mesh file
    basename = os.path.basename(meshFilenameDst)
    command = ['ln', '-s', meshFilenameDst, basename]
    subprocess.run(command, check=True)
    meshFilenameDst = basename

    # 1) Regrid the ice presence from the input data mesh to the grid of choice
    print("Regrid to desired mesh...")
    filenameOut = "icePresent_regrid.nc"

    meshFilenameSrc = os.path.join(args.dataDir, 'seaice_QU60km_polar.nc')
    filenameData = os.path.join(args.dataDir, 'icePresent_QU60km_polar.nc')

    regrid_to_other_mesh(
        meshFilenameSrc=meshFilenameSrc,
        filenameData=filenameData,
        meshFilenameDst=meshFilenameDst,
        filenameOut=filenameOut,
        generateWeights=True,
        weightsFilename=None)

    # 2) create icePresence variable
    print("fix_regrid_output...")

    inputFile = "icePresent_regrid.nc"
    outputFile = "icePresent_regrid_modify.nc"
    subprocess.call(["fix_regrid_output.exe", inputFile, meshFilenameDst,
                     outputFile])

    # 3) create variable icePresenceExtended
    print("extend_seaice_mask...")
    filenamePresence = "icePresent_regrid_modify.nc"
    extend_seaice_mask(meshFilenameDst, filenamePresence, 0.0, False)

    # 4) Make the regions file from the icePresenceExtended variable
    print("make_regions_file...")
    filenameIcePresent = "icePresent_regrid_modify.nc"
    filenameOut = "regions.nc"
    make_regions_file(filenameIcePresent=filenameIcePresent,
                      filenameMesh=meshFilenameDst,
                      regionType="two_region_eq",
                      varname="icePresenceExtended",
                      limit=0.5,
                      filenameOut=filenameOut)

    nProcsArray = args.nProcsArray

    # create partitions
    regionFilename = "regions.nc"
    outputPrefix = os.path.join(cwd, args.outputPrefix)

    gen_seaice_mesh_partition(meshFilename=meshFilenameDst,
                              regionFilename=regionFilename,
                              nProcsArray=nProcsArray,
                              mpasCullerLocation=None,
                              outputPrefix=outputPrefix,
                              plotting=False,
                              metis="gpmetis",
                              cullEquatorialRegion=False)

    for file in glob.glob(f'{outputPrefix}*'):
        command = ['chmod', 'ug+rw', file]
        subprocess.run(command, check=True)
        command = ['chmod', 'o+r', file]
        subprocess.run(command, check=True)

    os.chdir(cwd)
    shutil.rmtree(tmpdir)


# ---------------------------------------------------------------------
# Private functions
# ---------------------------------------------------------------------

def _degree_to_radian(degree):

    return (degree * math.pi) / 180.0


def _add_cell_cull_array(filename, cullCell):

    mesh = Dataset(filename, "a")

    cullCellVariable = mesh.createVariable("cullCell", "i4", ("nCells",))

    cullCellVariable[:] = cullCell

    mesh.close()


def _cull_mesh(meshToolsDir, filenameIn, filenameOut, cullCell):

    _add_cell_cull_array(filenameIn, cullCell)
    executable = "MpasCellCuller.x"
    if meshToolsDir is not None:
        executable = os.path.join(meshToolsDir, executable)

    subprocess.run([executable, filenameIn, filenameOut, "-c"], check=True)


def _load_partition(graphFilename):

    graphFile = open(graphFilename, "r")
    lines = graphFile.readlines()
    graphFile.close()

    partition = []
    for line in lines:
        partition.append(int(line))

    return partition


def _get_cell_ids(culledFilename, originalFilename):

    culledFile = Dataset(culledFilename, "r")
    nCellsCulled = len(culledFile.dimensions["nCells"])

    originalFile = Dataset(originalFilename, "r")
    nCellsOriginal = len(originalFile.dimensions["nCells"])

    cellid = np.zeros(nCellsCulled, dtype=int)

    cellMapFile = open("cellMapForward.txt", "r")
    cellMapLines = cellMapFile.readlines()
    cellMapFile.close()

    iCellOriginal = 0
    for cellMapLine in cellMapLines:

        if iCellOriginal % 1000 == 0:
            print(iCellOriginal, " of ", nCellsOriginal)

        try:
            cellMap = int(cellMapLine)
        except ValueError:
            # There are blank lines to skip
            continue

        if cellMap != -1:

            cellid[cellMap] = iCellOriginal

        iCellOriginal = iCellOriginal + 1

    return cellid


def _get_cell_ids_orig(culledFilename, originalFilename):

    culledFile = Dataset(culledFilename, "r")
    latCellCulled = culledFile.variables["latCell"][:]
    lonCellCulled = culledFile.variables["lonCell"][:]
    nCellsCulled = len(culledFile.dimensions["nCells"])

    originalFile = Dataset(originalFilename, "r")
    latCellOriginal = originalFile.variables["latCell"][:]
    lonCellOriginal = originalFile.variables["lonCell"][:]
    nCellsOriginal = len(originalFile.dimensions["nCells"])

    cellid = np.zeros(nCellsCulled, dtype=int)

    for iCellCulled in range(0, nCellsCulled):

        if iCellCulled % 1000 == 0:
            print("iCellCulled: ", iCellCulled, "of ", nCellsCulled)

        for iCellOriginal in range(0, nCellsOriginal):

            if (latCellCulled[iCellCulled] == latCellOriginal[iCellOriginal] and
                    lonCellCulled[iCellCulled] == lonCellOriginal[iCellOriginal]):

                cellid[iCellCulled] = iCellOriginal
                break

    return cellid
