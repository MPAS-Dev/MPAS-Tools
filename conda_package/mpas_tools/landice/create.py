import sys
from datetime import datetime
from optparse import OptionParser

import numpy
from netCDF4 import Dataset

from mpas_tools.cime.constants import constants


def create_from_generic_mpas_grid():  # noqa C901
    """
    Script to create a grid with land ice variables from an MPAS grid.
    Currently variable attributes are not copied.
    This script could be modified to copy them (looping over dir(var), skipping
    over variable function names "assignValue", "getValue", "typecode").
    """

    # earth radius, if needed
    sphere_radius = constants['SHR_CONST_REARTH']

    print(
        '** Gathering information.  (Invoke with --help for more details. '
        'All arguments are optional)'
    )
    parser = OptionParser()
    parser.add_option(
        '-i',
        '--in',
        dest='fileinName',
        help="input filename.  Defaults to 'grid.nc'",
        metavar='FILENAME',
    )
    parser.add_option(
        '-o',
        '--out',
        dest='fileoutName',
        help="output filename.  Defaults to 'landice_grid.nc'",
        metavar='FILENAME',
    )
    parser.add_option(
        '-l',
        '--level',
        dest='levels',
        help='Number of vertical levels to use in the output file.  '
        'Defaults to the number in the input file',
        metavar='FILENAME',
    )
    parser.add_option(
        '-v',
        '--vert',
        dest='vertMethod',
        help='Method of vertical layer spacing: uniform, glimmer.  '
        'Glimmer spacing follows Eq. 35 of Rutt, I. C., M. Hagdorn, '
        'N. R. J. Hulton, and A. J. Payne (2009), The Glimmer community '
        'ice sheet model, J. Geophys. Res., 114, F02004, '
        'doi:10.1029/2008JF001015',
        default='glimmer',
        metavar='FILENAME',
    )
    parser.add_option(
        '--beta', dest='beta', action='store_true', help='DEPRECATED'
    )
    parser.add_option(
        '--effecpress',
        dest='effecpress',
        action='store_true',
        help='DEPRECATED',
    )
    parser.add_option(
        '--diri',
        dest='dirichlet',
        action='store_true',
        help="Use this flag to include the fields 'dirichletVelocityMask', "
        "'uReconstructX', 'uReconstructY' needed for specifying Dirichlet "
        'velocity boundary conditions in the resulting file.',
    )
    parser.add_option(
        '--thermal',
        dest='thermal',
        action='store_true',
        help="Use this flag to include the fields 'temperature', "
        "'surfaceAirTemperature', 'basalHeatFlux' needed for specifying "
        'thermal initial conditions in the resulting file.',
    )
    parser.add_option(
        '--hydro',
        dest='hydro',
        action='store_true',
        help="Use this flag to include the fields 'waterThickness', "
        "'tillWaterThickness', 'basalMeltInput', 'externalWaterInput', "
        "'frictionAngle', 'waterPressure', 'waterFluxMask' needed for "
        'specifying hydro initial conditions in the resulting file.',
    )
    parser.add_option(
        '--obs',
        dest='obs',
        action='store_true',
        help='Use this flag to include the observational fields '
        'observedSurfaceVelocityX, observedSurfaceVelocityY, '
        'observedSurfaceVelocityUncertainty, observedThicknessTendency, '
        'observedThicknessTendencyUncertainty, thicknessUncertainty '
        'needed doing optimizations constrained by obs velocities.',
    )
    options, args = parser.parse_args()

    if not options.fileinName:
        print("No input filename specified, so using 'grid.nc'.")
        options.fileinName = 'grid.nc'
    else:
        print(f'Input file is: {options.fileinName}')
    if not options.fileoutName:
        print("No output filename specified, so using 'landice_grid.nc'.")
        options.fileoutName = 'landice_grid.nc'
    # make a space in stdout before further output
    print('')

    # Get the input file
    filein = Dataset(options.fileinName, 'r')

    # Define the new file to be output
    fileout = Dataset(options.fileoutName, 'w', format=filein.file_format)

    # ============================================
    # Copy over all of the netcdf global attributes
    # ============================================
    # Do this first as doing it last is slow for big files since adding
    # attributes forces the contents to get reorganized.
    print('---- Copying global attributes from input file to output file ----')
    for name in filein.ncattrs():
        # sphere radius needs to be set to that of the earth if on a sphere
        if (
            name == 'sphere_radius'
            and filein.on_a_sphere == 'YES             '
        ):
            fileout.sphere_radius = sphere_radius
            print(f'Set global attribute   sphere_radius = {sphere_radius}')
        elif name == 'history':
            # Update history attribute of netCDF file
            newhist = '\n'.join([filein.history, ' '.join(sys.argv[:])])
            fileout.history = newhist
        else:
            # Otherwise simply copy the attr
            setattr(fileout, name, getattr(filein, name))
            print(f'Copied global attribute  {name} = {getattr(filein, name)}')

    # Update history attribute of netCDF file if we didn't above
    if not hasattr(fileout, 'history'):
        fileout.history = sys.argv[:]

    fileout.sync()
    # make a space in stdout before further output
    print('')

    # ============================================
    # Copy over all the dimensions to the new file
    # ============================================
    # Note: looping over dimensions seems to result in them being written in
    #       seemingly random order.  don't think this matters but it is not
    #       aesthetically pleasing. It may be better to list them explicitly
    #       as I do for the grid variables, but this way ensures they all get
    #       included and is easier.
    # Note: The UNLIMITED time dimension will return a dimension value of None
    #       with Scientific.IO.  This is what is supposed to happen.  See below
    #       for how to deal with assigning values to a variable with a
    #       unlimited dimension.  Special handling is needed with the netCDF
    #       module.
    print('---- Copying dimensions from input file to output file ----')

    create_dims = list(filein.dimensions.keys())

    # We will add variables with Time as a dimension so we need it
    if 'Time' not in create_dims:
        create_dims.append('Time')

    for dim in create_dims:
        if dim == 'nTracers':
            pass  # Do nothing - we don't want this dimension
        elif dim == 'nVertInterfaces':
            pass  # Do nothing - this dimension will be handled below
        else:  # Copy over all other dimensions
            if dim == 'Time':
                # netCDF4 won't properly get this with the command below (you
                # need to use the isunlimited method)
                dimvalue = None
            elif dim == 'nVertLevels':
                if options.levels is None:
                    # If nVertLevels is in the input file, and a value for it
                    # was not specified on the command line, then use the value
                    # from the file (do nothing here)
                    print(
                        f'Using nVertLevels from the intput file: '
                        f'{len(filein.dimensions[dim])}'
                    )
                    dimvalue = len(filein.dimensions[dim])
                else:
                    # if nVertLevels is in the input file, but a value WAS
                    # specified on the command line, then use the command line
                    # value
                    print(
                        f'Using nVertLevels specified on the command line: '
                        f'{int(options.levels)}'
                    )
                    dimvalue = int(options.levels)
            else:
                dimvalue = len(filein.dimensions[dim])
            fileout.createDimension(dim, dimvalue)
    # There may be input files that do not have nVertLevels specified, in
    # which case it has not been added to the output file yet.  Treat those
    # here.
    if 'nVertLevels' not in fileout.dimensions:
        if options.levels is None:
            print(
                'nVertLevels not in input file and not specified.  Using '
                'default value of 10.'
            )
            fileout.createDimension('nVertLevels', 10)
        else:
            print(
                f'Using nVertLevels specified on the command line: '
                f'{int(options.levels)}'
            )
            fileout.createDimension('nVertLevels', int(options.levels))
    # Also create the nVertInterfaces dimension, even if none of the variables
    # require it.
    # nVertInterfaces = nVertLevels + 1
    fileout.createDimension(
        'nVertInterfaces', len(fileout.dimensions['nVertLevels']) + 1
    )
    print(
        f'Added new dimension nVertInterfaces to output file with value of '
        f'{len(fileout.dimensions["nVertInterfaces"])}'
    )

    fileout.sync()
    # include an extra blank line here
    print('Finished creating dimensions in output file.\n')

    # ============================================
    # Copy over all of the required grid variables to the new file
    # ============================================
    print('Beginning to copy mesh variables to output file.')
    vars2copy = [
        'latCell',
        'lonCell',
        'xCell',
        'yCell',
        'zCell',
        'indexToCellID',
        'latEdge',
        'lonEdge',
        'xEdge',
        'yEdge',
        'zEdge',
        'indexToEdgeID',
        'latVertex',
        'lonVertex',
        'xVertex',
        'yVertex',
        'zVertex',
        'indexToVertexID',
        'cellsOnEdge',
        'nEdgesOnCell',
        'nEdgesOnEdge',
        'edgesOnCell',
        'edgesOnEdge',
        'weightsOnEdge',
        'dvEdge',
        'dcEdge',
        'angleEdge',
        'areaCell',
        'areaTriangle',
        'cellsOnCell',
        'verticesOnCell',
        'verticesOnEdge',
        'edgesOnVertex',
        'cellsOnVertex',
        'kiteAreasOnVertex',
    ]
    # Add these optional fields if they exist in the input file
    for optionalVar in [
        'meshDensity',
        'gridSpacing',
        'cellQuality',
        'triangleQuality',
        'triangleAngleQuality',
        'obtuseTriangle',
    ]:
        if optionalVar in filein.variables:
            vars2copy.append(optionalVar)

    for _ in vars2copy:
        print('- ', end='')
    print('|')
    for varname in vars2copy:
        thevar = filein.variables[varname]
        datatype = thevar.dtype
        newVar = fileout.createVariable(varname, datatype, thevar.dimensions)
        if filein.on_a_sphere == 'YES             ':
            if varname in (
                'xCell',
                'yCell',
                'zCell',
                'xEdge',
                'yEdge',
                'zEdge',
                'xVertex',
                'yVertex',
                'zVertex',
                'dvEdge',
                'dcEdge',
            ):
                newVar[:] = thevar[:] * sphere_radius / filein.sphere_radius
            elif varname in ('areaCell', 'areaTriangle', 'kiteAreasOnVertex'):
                newVar[:] = (
                    thevar[:] * (sphere_radius / filein.sphere_radius) ** 2
                )
            else:
                newVar[:] = thevar[:]
        else:
            # not on a sphere
            newVar[:] = thevar[:]
        del newVar, thevar
        sys.stdout.write('* ')
        sys.stdout.flush()
    fileout.sync()
    print('|')
    print('Finished copying mesh variables to output file.\n')

    # ============================================
    # Create the land ice variables (all the shallow water vars in the input
    # file can be ignored)
    # ============================================
    nVertLevels = len(fileout.dimensions['nVertLevels'])
    # Get the datatype for double precision float
    datatype = filein.variables['xCell'].dtype
    # Get the datatype for integers
    datatypeInt = filein.variables['indexToCellID'].dtype
    #  Note: it may be necessary to make sure the Time dimension has size 1,
    #        rather than the 0 it defaults to.  For now, letting it be 0 which
    #        seems to be fine.

    # layerThicknessFractions
    layerThicknessFractions = fileout.createVariable(
        'layerThicknessFractions', datatype, ('nVertLevels',)
    )
    layerThicknessFractionsData = numpy.zeros(layerThicknessFractions.shape)
    # Assign default values to layerThicknessFractions.  By default they will
    # be uniform fractions.  Users can modify them in a subsequent step, but
    # doing this here ensures the most likely values are already assigned.
    # (Useful for e.g. setting up Greenland where the state variables are
    # copied over but the grid variables are not modified.)
    if options.vertMethod == 'uniform':
        layerThicknessFractionsData[:] = 1.0 / nVertLevels
    elif options.vertMethod == 'glimmer':
        nInterfaces = nVertLevels + 1
        layerInterfaces = numpy.zeros((nInterfaces,))
        for k in range(nInterfaces):
            layerInterfaces[k] = (
                4.0
                / 3.0
                * (1.0 - ((k + 1.0 - 1.0) / (nInterfaces - 1.0) + 1.0) ** -2)
            )
        for k in range(nVertLevels):
            layerThicknessFractionsData[k] = (
                layerInterfaces[k + 1] - layerInterfaces[k]
            )
        print(
            f'Setting layerThicknessFractions to: '
            f'{layerThicknessFractionsData}'
        )
    else:
        raise ValueError(
            f'Unknown method for vertical spacing method '
            f'(--vert): {options.vertMethod}'
        )

    # explictly specify layer fractions
    # layerThicknessFractionsData[:] = [
    #     0.1663,
    #     0.1516,
    #     0.1368,
    #     0.1221,
    #     0.1074,
    #     0.0926,
    #     0.0779,
    #     0.0632,
    #     0.0484,
    #     0.0337,
    # ]

    layerThicknessFractions[:] = layerThicknessFractionsData[:]

    # With Scientific.IO.netCDF, entries are appended along the unlimited
    # dimension one at a time by assigning to a slice. Therefore we need to
    # assign to time level 0, and what we need to assign is a zeros array that
    # is the shape of the new variable, exluding the time dimension!
    newvar = fileout.createVariable('thickness', datatype, ('Time', 'nCells'))
    newvar[0, :] = numpy.zeros(newvar.shape[1:])
    # These landice variables are stored in the mesh currently, and therefore
    # do not have a time dimension. It may make sense to eventually move them
    # to state.
    newvar = fileout.createVariable(
        'bedTopography', datatype, ('Time', 'nCells')
    )
    newvar[:] = numpy.zeros(newvar.shape)
    newvar = fileout.createVariable('sfcMassBal', datatype, ('Time', 'nCells'))
    newvar[:] = numpy.zeros(newvar.shape)
    newvar = fileout.createVariable(
        'floatingBasalMassBal', datatype, ('Time', 'nCells')
    )
    newvar[:] = numpy.zeros(newvar.shape)
    print(
        'Added default variables: thickness, temperature, bedTopography, '
        'sfcMassBal, floatingBasalMassBal'
    )

    newvar = fileout.createVariable('beta', datatype, ('Time', 'nCells'))
    # Give a default beta that won't have much sliding.
    newvar[:] = 1.0e8
    print('Added variable: beta')

    newvar = fileout.createVariable('muFriction', datatype, ('Time', 'nCells'))
    # Give a default mu that won't have much sliding.
    newvar[:] = 1.0e8
    print('Added variable: muFriction')

    newvar = fileout.createVariable(
        'effectivePressure', datatype, ('Time', 'nCells')
    )
    # Give a default effective pressure of 1.0 so that, for the linear sliding
    # law, beta = mu*effecpress = mu.
    newvar[:] = 1.0
    print('Added variable: effectivePressure')

    newvar = fileout.createVariable(
        'stiffnessFactor', datatype, ('Time', 'nCells')
    )
    # Give default value
    newvar[:] = 1.0
    print('Added variable: stiffnessFactor')

    newvar = fileout.createVariable(
        'eigencalvingParameter', datatype, ('Time', 'nCells')
    )
    # Give default value for eigencalvingParameter
    newvar[:] = 3.14e16
    print('Added variable: eigencalvingParameter')

    newvar = fileout.createVariable(
        'groundedVonMisesThresholdStress', datatype, ('Time', 'nCells')
    )
    # Give default value
    newvar[:] = 1.0e6
    print('Added variable: groundedVonMisesThresholdStress')

    newvar = fileout.createVariable(
        'floatingVonMisesThresholdStress', datatype, ('Time', 'nCells')
    )
    # Give default value
    newvar[:] = 1.0e6
    print('Added variable: floatingVonMisesThresholdStress')

    newvar = fileout.createVariable('iceMask', datatype, ('Time', 'nCells'))
    newvar[:] = 0
    print('Added variable: iceMask')

    if options.dirichlet:
        newvar = fileout.createVariable(
            'dirichletVelocityMask',
            datatypeInt,
            ('Time', 'nCells', 'nVertInterfaces'),
        )
        # default: no Dirichlet b.c.
        newvar[:] = 0
        newvar = fileout.createVariable(
            'uReconstructX',
            datatype,
            (
                'Time',
                'nCells',
                'nVertInterfaces',
            ),
        )
        newvar[:] = 0.0
        newvar = fileout.createVariable(
            'uReconstructY',
            datatype,
            (
                'Time',
                'nCells',
                'nVertInterfaces',
            ),
        )
        newvar[:] = 0.0
        print(
            'Added optional dirichlet variables: dirichletVelocityMask, '
            'uReconstructX, uReconstructY'
        )

    if options.thermal:
        newvar = fileout.createVariable(
            'temperature', datatype, ('Time', 'nCells', 'nVertLevels')
        )
        # Give default value for temperate ice
        newvar[:] = 273.15
        newvar = fileout.createVariable(
            'surfaceAirTemperature', datatype, ('Time', 'nCells')
        )
        # Give default value for temperate ice
        newvar[:] = 273.15
        newvar = fileout.createVariable(
            'basalHeatFlux', datatype, ('Time', 'nCells')
        )
        # Default to none (W/m2)
        newvar[:] = 0.0
        print(
            'Added optional thermal variables: temperature, '
            'surfaceAirTemperature, basalHeatFlux'
        )

    if options.hydro:
        newvar = fileout.createVariable(
            'waterThickness', datatype, ('Time', 'nCells')
        )
        newvar[:] = 0.0
        newvar = fileout.createVariable(
            'tillWaterThickness', datatype, ('Time', 'nCells')
        )
        newvar[:] = 0.0
        newvar = fileout.createVariable(
            'basalMeltInput', datatype, ('Time', 'nCells')
        )
        newvar[:] = 0.0
        newvar = fileout.createVariable(
            'externalWaterInput', datatype, ('Time', 'nCells')
        )
        newvar[:] = 0.0
        newvar = fileout.createVariable(
            'frictionAngle', datatype, ('Time', 'nCells')
        )
        newvar[:] = 0.0
        newvar = fileout.createVariable(
            'waterPressure', datatype, ('Time', 'nCells')
        )
        newvar[:] = 0.0
        newvar = fileout.createVariable(
            'waterFluxMask', 'i', ('Time', 'nEdges')
        )
        newvar[:] = 0.0
        print(
            'Added optional hydro variables: waterThickness, '
            'tillWaterThickness, meltInput, frictionAngle, waterPressure, '
            'waterFluxMask'
        )

    if options.obs:
        newvar = fileout.createVariable(
            'observedSurfaceVelocityX', datatype, ('Time', 'nCells')
        )
        newvar[:] = 0.0
        newvar = fileout.createVariable(
            'observedSurfaceVelocityY', datatype, ('Time', 'nCells')
        )
        newvar[:] = 0.0
        newvar = fileout.createVariable(
            'observedSurfaceVelocityUncertainty', datatype, ('Time', 'nCells')
        )
        newvar[:] = 0.0
        newvar = fileout.createVariable(
            'observedThicknessTendency', datatype, ('Time', 'nCells')
        )
        newvar[:] = 0.0
        newvar = fileout.createVariable(
            'observedThicknessTendencyUncertainty',
            datatype,
            ('Time', 'nCells'),
        )
        newvar[:] = 0.0
        newvar = fileout.createVariable(
            'thicknessUncertainty', datatype, ('Time', 'nCells')
        )
        newvar[:] = 0.0
        print(
            'Added optional velocity optimization variables: '
            'observedSurfaceVelocityX, observedSurfaceVelocityY, '
            'observedSurfaceVelocityUncertainty, observedThicknessTendency, '
            'observedThicknessTendencyUncertainty, thicknessUncertainty'
        )

    # Update history attribute of netCDF file
    thiscommand = (
        datetime.now().strftime('%a %b %d %H:%M:%S %Y')
        + ': '
        + ' '.join(sys.argv[:])
    )
    if hasattr(fileout, 'history'):
        newhist = '\n'.join([thiscommand, fileout.history])
    else:
        newhist = thiscommand
    fileout.history = newhist

    print(
        'Completed creating land ice variables in new file. Now syncing to '
        'file.'
    )
    fileout.sync()

    filein.close()
    fileout.close()

    print(f'\n** Successfully created {options.fileoutName}.**')
