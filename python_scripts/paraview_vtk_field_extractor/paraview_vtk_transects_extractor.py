#!/usr/bin/env python
"""
Name: paraview_vtk_transects_extractor.py
Authors: Xylar Asay-Davis
Date: 04/17/2016

This script is used to extract a field from a time series of NetCDF files as
VTK files for plotting in paraview.

It can extract a field across multiple files by passing in a regular expression
for the filename patter. As an example, one can run the script using:

`./paraview_vtk_field_extractor.py -v areaCell,latVertex -f "hist.comp.*.nc"`

To extract a time series of areaCell,latVertex that spans multiple files.
By default, time-independent fields on cells are written to a file
vtk_files/staticFieldsOnCells.vtp
and time-dependent fields on cells are written to
vtk_files/timeDependentFieldsOnCells.pvd
vtk_files/time_series/timeDependentFieldsOnCells.0.vtp
vtk_files/time_series/timeDependentFieldsOnCells.1.vtp
...
and similarly for edges and vertices.  Time-independent fields can be
included in each time step of the time-dependent fields for with
the --combine flag.  This allows time-dependent and -independent fields
to be combined in filters within Paraview at the expense of considerable
additional storage space.


Requirements:
This script requires access to the following non standard modules:
pyevtk
netCDF4
numpy

Optional modules:
progressbar
"""
from pyevtk.vtk import VtkFile, VtkPolyData

import sys, os, glob
import numpy as np

from netCDF4 import Dataset as NetCDFFile
import argparse

try:
    from progressbar import ProgressBar, Percentage, Bar, ETA
    use_progress_bar = True
except:
    use_progress_bar = False
    
class Transect(object):
    def __init__(self, cellIndices, mesh_file, nLevels, minLevel, maxLevel, dtype):

        self.cellIndices = cellIndices
        self.nLevels = nLevels

        nCells = len(cellIndices)
        self.nCells = nCells
        
        self.minLevel = minLevel
        self.maxLevel = maxLevel

        cellMask = np.ones((nCells,nLevels), bool)
        for iLevel in range(nLevels):
            if minLevel is not None:
                cellMask[:,iLevel] = np.logical_and(cellMask[:,iLevel], iLevel >= minLevel)
            if maxLevel is not None:
                cellMask[:,iLevel] = np.logical_and(cellMask[:,iLevel], iLevel <= maxLevel)

        edgeMask = np.logical_or(cellMask[0:-1,:], cellMask[1:,:])
        
        
        self.mask = np.zeros((2*self.nCells-1,self.nLevels), dtype=bool)
        self.mask[0::2,:] = cellMask
        self.mask[1::2,:] = edgeMask


        cellsOnCell = mesh_file.variables['cellsOnCell'][cellIndices,:]-1
        edgesOnCell = mesh_file.variables['edgesOnCell'][cellIndices,:]-1

        edgeIndices = np.zeros(nCells-1,int)
        for indexInTransect in range(nCells-1):
            iCellNext = cellIndices[indexInTransect+1]
            iNeighbor = np.nonzero(cellsOnCell[indexInTransect,:] == iCellNext)[0][0]
            edgeIndices[indexInTransect] = edgesOnCell[indexInTransect,iNeighbor]


        x = np.zeros(2*nCells-1, dtype=dtype)
        y = np.zeros(2*nCells-1, dtype=dtype)
        z = np.zeros(2*nCells-1, dtype=dtype)
        transectCellIndices = np.zeros(2*nCells-2, dtype=int)

        x[0::2] = mesh_file.variables['xCell'][cellIndices]
        x[1::2] = mesh_file.variables['xEdge'][edgeIndices]
        y[0::2] = mesh_file.variables['yCell'][cellIndices]
        y[1::2] = mesh_file.variables['yEdge'][edgeIndices]
        z[0::2] = mesh_file.variables['zCell'][cellIndices]
        z[1::2] = mesh_file.variables['zEdge'][edgeIndices]
        transectCellIndices[0::2] = np.arange(0,nCells-1)
        transectCellIndices[1::2] = np.arange(1,nCells)

        pointCellIndices = []
        pointLevelIndices = []
        pointHorizIndices = []
        pointVertIndices = []
        X = []
        Y = []
        Z = []
        for index in range(2*nCells-2):
            iCell = transectCellIndices[index]
            for iLevel in range(nLevels):
                if not (self.mask[index,iLevel] and self.mask[index+1,iLevel]):
                    continue
                horizIndices = [index, index+1, index+1, index]
                vertIndices = [iLevel, iLevel, iLevel+1, iLevel+1]
                X.extend([x[i] for i in horizIndices])
                Y.extend([y[i] for i in horizIndices])
                Z.extend([z[i] for i in horizIndices])
                pointCellIndices.extend([iCell,iCell,iCell,iCell])
                pointLevelIndices.extend([iLevel,iLevel,iLevel,iLevel])
                pointHorizIndices.extend(horizIndices)
                pointVertIndices.extend(vertIndices)

        self.nPoints = len(X)
        self.nPolygons = self.nPoints/4
        
        self.points = (np.array(X, dtype=dtype), np.array(Y, dtype=dtype), np.array(Z, dtype=dtype))
        self.pointCellIndices = np.array(pointCellIndices, dtype=int)
        self.pointLevelIndices = np.array(pointLevelIndices, dtype=int)
        self.pointHorizIndices = np.array(pointHorizIndices, dtype=int)
        self.pointVertIndices = np.array(pointVertIndices, dtype=int)
        
        self.connectivity = np.arange(self.nPoints)
        self.offsets = 4 + 4*np.arange(self.nPolygons)
    
    def computeZInterface(self, layerThicknessCell, zMin, zMax, dtype):
        
        cellMask = self.mask[0::2,:]
        edgeMask = self.mask[1::2,:]
        
        if zMin is not None:
            zInterfaceCell = np.zeros((self.nCells,self.nLevels+1))
            zInterfaceCell[:,0] = zMin
            for iLevel in range(self.nLevels):
                zInterfaceCell[:,iLevel+1] = (zInterfaceCell[:,iLevel] 
                    + cellMask[:,iLevel]*layerThicknessCell[:,iLevel])
        else:
            zInterfaceCell = np.zeros((self.nCells,self.nLevels+1))
            zInterfaceCell[:,-1] = zMax
            for iLevel in range(self.nLevels-1,-1,-1):
                zInterfaceCell[:,iLevel] = (zInterfaceCell[:,iLevel+1]
                    - cellMask[:,iLevel]*layerThicknessCell[:,iLevel])
                    
        validCell = np.zeros(zInterfaceCell.shape,bool)
        validCell[:,0:-1] = cellMask
        validCell[:,1:] = np.logical_or(validCell[:,1:],cellMask)
        
                   
        denom = (1.0*validCell[0:-1,:] + 1.0*validCell[1:,:])
        
        zInterfaceEdge = (validCell[0:-1,:]*zInterfaceCell[0:-1,:]
                              + validCell[1:,:]*zInterfaceCell[1:,:])
        
        validEdge = np.zeros(zInterfaceEdge.shape,bool)
        validEdge[:,0:-1] = edgeMask
        validEdge[:,1:] = np.logical_or(validEdge[:,1:],edgeMask)

        zInterfaceEdge[validEdge] /= denom[validEdge]

        zInterface = np.zeros((2*self.nCells-1,self.nLevels+1), dtype=dtype)
        zInterface[0::2,:] = zInterfaceCell
        zInterface[1::2,:] = zInterfaceEdge
        
        return zInterface[self.pointHorizIndices, self.pointVertIndices]

def setup_time_indices(fn_pattern):#{{{
    # Build file list and time indices
    file_list = sorted(glob.glob(fn_pattern))

    local_indices = []
    file_names = []

    if len(file_list) == 0:
        print "No files to process."
        print "Exiting..."
        sys.exit(0)

    if use_progress_bar:
        widgets = ['Build time indices: ', Percentage(), ' ', Bar(), ' ', ETA()]
        time_bar = ProgressBar(widgets=widgets, maxval=len(file_list)).start()
    else:
        print "Build time indices..."

    i_file = 0
    for file in file_list:
        nc_file = NetCDFFile(file, 'r')

        try:
            times = len(nc_file.dimensions['Time'])
        except:
            times = 1

        times = max(times, 1)

        for time_idx in np.arange(0, times):
            local_indices.append(time_idx)
            file_names.append(file)

        i_file = i_file + 1
        nc_file.close()
        if use_progress_bar:
            time_bar.update(i_file)

    if use_progress_bar:
        time_bar.finish()

    return (local_indices, file_names)
#}}}

def parse_extra_dim(dimName, indexString, time_series_file, mesh_file):#{{{
    if indexString == '':
        return np.zeros(0,int)
        
    if (mesh_file is not None) and (dimName in mesh_file.dimensions):
        nc_file = mesh_file
    else:
        nc_file = time_series_file
    dimSize = len(nc_file.dimensions[dimName])
    if ',' in indexString:
        indices = []
        for index in indexString.split(','):
            indices.append(int(index))
        indices = np.array(indices,int)
    elif ':' in indexString:
        index_list = indexString.split(':')
        if len(index_list) in [2,3]:
            if index_list[0] == '':
                first = 0
            else:
                first = int(index_list[0])
            if index_list[1] == '':
                last = dimSize
            else:
                last = int(index_list[1])
            if (len(index_list) == 2) or (index_list[2] == ''):
                step = 1
            else:
                step = int(index_list[2])
            indices = np.arange(first,last,step)
        else:
            print "Improperly formatted extra dimension:", dimName, indexString
            return None
    else:
        indices = np.array([int(indexString)])

    valid = np.all(np.logical_and(indices >= 0,indices <= dimSize-1))
    if not valid:
        print "Index (or indices) out of bounds for extra dimension:", dimName, indexString
        return None

    return indices

#}}}

def parse_extra_dims(dimension_list, time_series_file, mesh_file):#{{{

    extra_dims = {}
    
    if dimension_list:
        for dim_item in dimension_list:
            (dimName,indexString) = dim_item.split('=')
            indices = parse_extra_dim(dimName, indexString, time_series_file, mesh_file)
            if indices is not None and len(indices) <= 1:
                extra_dims[dimName] = indices

    return extra_dims

#}}}

def setup_dimension_values_and_sort_vars(time_series_file, mesh_file, variable_list, 
                                         transect_dim_name, extra_dims):#{{{

    if variable_list == 'all':
        variables = []
        files = [time_series_file]
        if mesh_file is not None:
            files.append(mesh_file)
        for nc_file in files:
            for var in nc_file.variables:
                dims = nc_file.variables[var].dimensions
                if 'nCells' in dims:
                    variables.append(str(var))
    else:
        variables = variable_list.split(',')

    # make sure the variables are unique
    variables = list(set(variables))
    variables.sort()

    for var in variables:
        if (mesh_file is not None) and (var in mesh_file.variables):
            nc_file = mesh_file
        else:
            nc_file = time_series_file
        field_dims = nc_file.variables[var].dimensions
        localExtraDims = []
        for dim in field_dims:
            if dim not in ['Time', 'nCells', 'nEdges', 'nVertices', transect_dim_name]:
                localExtraDims.append(str(dim))
                        
        for dim in localExtraDims:
            if dim not in extra_dims:
                # ignored variables with this dimension
                extra_dims[dim] = []

    cellVars = []

    for var in variables:
        # first try the mesh file, then the time series
        if (mesh_file is not None) and (var in mesh_file.variables):
            nc_file = mesh_file
        else:
            nc_file = time_series_file
        field_var = nc_file.variables[var]

        field_dims = field_var.dimensions
        
        skip = False
        for dim in field_dims:
            if dim in extra_dims and len(extra_dims[dim]) == 0:
                # this variable has an empty dimension, so skip it
                skip = True
                break

        if skip:
            continue

        if "nCells" in field_dims:
            cellVars.append(var)

    return (extra_dims, cellVars)
#}}}

def summarize_extraction(mesh_file, transects_file, time_indices, cellVars):#{{{
    print ""
    print "Extracting a total of %d time levels."%(len(time_indices))
    print "Using file '%s' as the mesh file and '%s' "%(mesh_file, transects_file)
    print "as the transects file for this extraction."
    print ""
    print "Transects of the following variables will be extracted from the input file(s)."
    print ""

    if len(cellVars) > 0:
        print "   Variables with 'nCells' as a dimension:"
        for var in cellVars:
            print "      name: %s"%(var)

    print ""
#}}}

def build_transects( mesh_file, transects_file, time_series_file, 
                     transect_dim_name, min_level_name, max_level_name, 
                     output_32bit ):#{{{

    if output_32bit:
        outType = 'f4'
    else:
        outType = 'f8'

    nTransects = len(transects_file.dimensions['nTransects'])
    nTransectLevels = len(time_series_file.dimensions[transect_dim_name])

    if (mesh_file is not None):
        nc_file = mesh_file
    else:
        nc_file = time_series_file        

    # list of points starts out empty; will add points from each transect
    transects = []

    for iTransect in range(nTransects):
        # get the global indices of this transect
        cellIndices = transects_file.variables['transectCellGlobalIDs'][iTransect,:]-1
        cellIndices = cellIndices[cellIndices >= 0]
        if min_level_name is not None:
            minLevelCell = time_series_file.variables[min_level_name][cellIndices]-1
        else:
            minLevelCell = None
        if max_level_name is not None:
            maxLevelCell = time_series_file.variables[max_level_name][cellIndices]-1
        else:
            maxLevelCell = None
        
        transect = Transect(cellIndices, nc_file, nTransectLevels, minLevelCell, 
                            maxLevelCell, outType)
                            
        transects.append(transect)    
        
    return transects

#}}}

def build_transects_time_series( local_time_indices, file_names, mesh_file, extra_dims,
                                 transect_dim_name, transects, variable_names,
                                 layer_thickness_name, z_min_name, z_max_name,
                                 output_32bit, combine_output, 
                                 append ):#{{{
    def write_pvd_header(prefix):
        pvd_file = open('vtk_files/%s.pvd'%(prefix), 'w')
        pvd_file.write('<?xml version="1.0"?>\n')
        pvd_file.write('<VTKFile type="Collection" version="0.1"\n')
        pvd_file.write('\tbyte_order="LittleEndian"\n')
        pvd_file.write('\tcompressor="vtkZLibDataCompressor">\n')
        return pvd_file

    def write_vtp_header(prefix, activeVarIndex, varIndices, transect):
        vtkFile = VtkFile("vtk_files/%s"%prefix, VtkPolyData)
        vtkFile.openElement(vtkFile.ftype.name)
        vtkFile.openPiece(npoints=transect.nPoints,npolys=transect.nPolygons)

        vtkFile.openElement("Points")
        vtkFile.addData("points", transect.points)
        vtkFile.closeElement("Points")

        vtkFile.openElement("Polys")
        vtkFile.addData("connectivity", transect.connectivity)
        vtkFile.addData("offsets", transect.offsets)
        vtkFile.closeElement("Polys")

        vtkFile.openData("Point", scalars = variable_names[activeVarIndex])
        for iVar in varIndices:
            vtkFile.addHeader(variable_names[iVar], outType, transect.nPoints, 1)
        vtkFile.addHeader('zInterface', outType, transect.nPoints, 1)
        vtkFile.closeData("Point")

        vtkFile.closePiece()
        vtkFile.closeElement(vtkFile.ftype.name)

        vtkFile.appendData(transect.points)
        vtkFile.appendData(transect.connectivity)
        vtkFile.appendData(transect.offsets)

        return vtkFile
        
    def read_transect_field(transect, field_name):
        if field_name[0] == '-':
            field_name = field_name[1:]
            sign = -1
        else:
            sign = 1
            
        field_var = variables[field_name]
        field_ndims = len(field_var.dimensions)
        if transect_dim_name in field_var.dimensions:
            field = np.zeros((transect.nCells,transect.nLevels), dtype=outType)
            
            for iCellTransect in range(transect.nCells):
                iCell = transect.cellIndices[iCellTransect]

                dim_vals = []
                for iDim in range(field_ndims):
                    dim = field_var.dimensions[iDim]
                    if dim == 'Time':
                        dim_vals.append(local_time_indices[time_index])
                    elif dim == 'nCells':
                        dim_vals.append(iCell)
                    elif dim == transect_dim_name:
                        dim_vals.append(np.arange(transect.nLevels))
                    else:
                        dim_vals.append(extra_dims[dim][0])


                if field_ndims == 1:
                    field[iCellTransect,:] = field_var[dim_vals[0]]
                elif field_ndims == 2:
                    field[iCellTransect,:]= field_var[dim_vals[0], dim_vals[1]]
                elif field_ndims == 3:
                    field[iCellTransect,:]= field_var[dim_vals[0], dim_vals[1], dim_vals[2]]
                elif field_ndims == 4:
                    field[iCellTransect,:]= field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3]]
                elif field_ndims == 5:
                    field[iCellTransect,:]= field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3], dim_vals[4]]

        else:
            field = np.zeros((transect.nCells), dtype=outType)

            dim_vals = []
            for iDim in range(field_ndims):
                dim = field_var.dimensions[iDim]
                if dim == 'Time':
                    dim_vals.append(local_time_indices[time_index])
                elif dim == 'nCells':
                    dim_vals.append(transect.cellIndices)
                else:
                    dim_vals.append(extra_dims[dim][0])

            if field_ndims == 1:
                field[:] = field_var[dim_vals[0]]
            elif field_ndims == 2:
                field[:]= field_var[dim_vals[0], dim_vals[1]]
            elif field_ndims == 3:
                field[:]= field_var[dim_vals[0], dim_vals[1], dim_vals[2]]
            elif field_ndims == 4:
                field[:]= field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3]]
            elif field_ndims == 5:
                field[:]= field_var[dim_vals[0], dim_vals[1], dim_vals[2], dim_vals[3], dim_vals[4]]


        return sign*field

    if len(variable_names) == 0:
        return

    if output_32bit:
        outType = 'float32'
    else:
        outType = 'float64'


    # Get dimension info to allocate the size of Colors
    time_series_file = NetCDFFile(file_names[0], 'r')

    nTransects = len(transects)
    nTimes = len(local_time_indices)

        
    nVars = len(variable_names)
    var_has_time_dim = np.zeros(nVars,bool)
    variables = {}
    variableExtraDims = {}
    for iVar in range(nVars):
        var_name = variable_names[iVar]
        if (mesh_file is not None) and (var_name in mesh_file.variables):
            var = mesh_file.variables[var_name]
            if 'Time' in var.dimensions:
                # we can't support time dependence in the mesh file
                var = time_series_file.variables[var_name]

        else:
            var = time_series_file.variables[var_name]

        variables[var_name] = var

        dims = var.dimensions

        var_has_time_dim[iVar] = 'Time' in dims
 
        variableExtraDims[var_name] = []
        for dim in dims:
            if dim not in ['nCells','nVertices','nEdges','Time']:
                variableExtraDims[var_name].append(dim)
        
    time_series_file.close()
    
    any_var_has_time_dim = np.any(var_has_time_dim)

    try:
        os.makedirs('vtk_files')
    except OSError:
        pass

    if any_var_has_time_dim:
        try:
            os.makedirs('vtk_files/time_series')
        except OSError:
            pass
    else:
        # there is no point in combining output if no fields have Time dim
        combine_output = False
        nTimes = 1

    # Output time series
    if use_progress_bar:
        widgets = ['Writing time series: ', Percentage(), ' ', Bar(), ' ', ETA()]
        field_bar = ProgressBar(widgets=widgets, maxval=nTimes*nVars*nTransects).start()
    else:
        print "Writing time series...."
        
    pad = np.array(np.floor(np.log10(nTransects)),int)+1
    template = '%%0%dd'%pad


    for iTransect in range(nTransects):
        transect = transects[iTransect]
        
        if any_var_has_time_dim:
            if combine_output or np.all(var_has_time_dim):
                out_prefix = "transect_%s"%(template%iTransect)
            else:
                out_prefix = "timeDependentTransect_%s"%(template%iTransect)
            # start the pvd file
            pvd_file = write_pvd_header(out_prefix)
            pvd_file.write('<Collection>\n')

        if not combine_output and not np.all(var_has_time_dim):
            out_prefix = "staticTransect_%s"%(template%iTransect)
            varIndices = np.arange(nVars)[var_has_time_dim == False]
            timeIndependentFile = write_vtp_header(out_prefix, varIndices[0], varIndices, transect)

        prev_file = ""
        for time_index in range(nTimes):
    
            if prev_file != file_names[time_index]:
                if prev_file != "":
                    time_series_file.close()
                time_series_file = NetCDFFile(file_names[time_index], 'r')
                prev_file = file_names[time_index]

            if any_var_has_time_dim:
                # write the header for the vtp file
                vtp_file_prefix = "time_series/%s.%d"%(out_prefix, time_index)
                file_name = 'vtk_files/%s.vtp'%(vtp_file_prefix)
                if append and os.path.exists(file_name):
                    continue
    
                if combine_output:
                    varIndices = np.arange(nVars)
                else:
                    varIndices = np.arange(nVars)[var_has_time_dim]
                timeDependentFile = write_vtp_header(vtp_file_prefix, varIndices[0], varIndices, transect)
    
                # add time step to pdv file
                pvd_file.write('<DataSet timestep="%d" group="" part="0"\n'%(time_index))
                pvd_file.write('\tfile="%s.vtp"/>\n'%(vtp_file_prefix))
            
            for iVar in range(nVars):
                var_name = variable_names[iVar]
                field_var = variables[var_name]

                has_time = 'Time' in field_var.dimensions
                if not combine_output and not has_time and time_index > 0:
                    continue
                
                if has_time or combine_output:
                    vtkFile = timeDependentFile
                else:
                    vtkFile = timeIndependentFile
                    
                field = read_transect_field(transect, var_name)
                
                if(len(field.shape) == 2):
                    field = field[transect.pointCellIndices,transect.pointLevelIndices]
                else:
                    field = field[transect.pointCellIndices]

                vtkFile.appendData(field)
                
                if use_progress_bar:
                    field_bar.update((iTransect*nTimes + time_index)*nVars + iVar)

                del field

            # build a cum sum of layer thickness for transect vertical coord
                
            vtkFile = timeDependentFile

            layerThickness = read_transect_field(transect, layer_thickness_name)
            zMin = None
            zMax = None
            if z_min_name is not None:
                zMin = read_transect_field(transect, z_min_name)
            elif z_max_name is not None:
                zMax = read_transect_field(transect, z_max_name)
            else:
                zMax = np.zeros(transect.nCells)

            zInterface = transect.computeZInterface(layerThickness, zMin, zMax, dtype=outType)
                            
            vtkFile.appendData(zInterface)
                
            if any_var_has_time_dim:
                timeDependentFile.save()
                del timeDependentFile
    
            if time_index == 0 and not combine_output and not np.all(var_has_time_dim):
                timeIndependentFile.save()
                del timeIndependentFile
    
        time_series_file.close()

    
        if any_var_has_time_dim:
            # finish the pdv file
            pvd_file.write('</Collection>\n')
            pvd_file.write('</VTKFile>\n')

    if use_progress_bar:
        field_bar.finish()
#}}}


if __name__ == "__main__":
    if use_progress_bar:
        print " -- Using progress bars --"
    else:
        print " -- Progress bars are not available--"
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--file_pattern", dest="filename_pattern", help="MPAS Filename pattern.", metavar="FILE", required=True)
    parser.add_argument("--transects_file", dest="transects_filename", help="MPAS transects filename.", required=True)
    parser.add_argument("--transect_dim", dest="transect_dim_name", help="A dimension for transect layers.", required=True)
    parser.add_argument("--layer_thickness", dest="layer_thickness_name", help="Variable for layer thickness, used to compute zInterface", required=True)
    parser.add_argument("--min_level_cell", dest="min_level_cell_name", help="Index array indicating the minimum valid layer in a cell (default is 0 for all cells)")
    parser.add_argument("--max_level_cell", dest="max_level_cell_name", help="Index array indicating the maximum valid layer in a cell (default is the transect_dim-1 for all cells)")
    parser.add_argument("--z_min", dest="z_min_name", help="Variable specifying the depth of the lower interface (in index space) of the first valid layer on cells, used to compute zInterface")
    parser.add_argument("--z_max", dest="z_max_name", help="Variable specifying the depth of the upper interface (in index space) of the last valid layer on cells, used to compute zInterface")
    parser.add_argument("-d", "--dim_list", dest="dimension_list", nargs="+", help="A list of dimensions and associated indices.")
    parser.add_argument("-m", "--mesh_file", dest="mesh_filename", help="MPAS Mesh filename. If not set, it will use the first file in the -f flag as the mesh file.")
    parser.add_argument("-v", "--variable_list", dest="variable_list", help="List of variables to extract ('all' for all variables, 'allOnCells' for all variables on cells, etc.)", metavar="VAR", required=True)
    parser.add_argument("-3", "--32bit", dest="output_32bit", help="If set, the vtk files will be written using 32bit floats.", action="store_true")
    parser.add_argument("-c", "--combine", dest="combine_output", help="If set, time-independent fields are written to each file along with time-dependent fields.", action="store_true")
    parser.add_argument("-a", "--append", dest="append", help="If set, only vtp files that do not already exist are written out.", action="store_true")
    args = parser.parse_args()

    if not args.output_32bit:
        use_32bit = False
    else:
        use_32bit = True

    (time_indices, time_file_names) = setup_time_indices(args.filename_pattern)

    separate_mesh_file = True
    if not args.mesh_filename:
        args.mesh_filename = time_file_names[0]
        separate_mesh_file = False
        
    transects_file = NetCDFFile(args.transects_filename, 'r')

    # Setting dimension values:
    time_series_file = NetCDFFile(time_file_names[0], 'r')
    if separate_mesh_file:
        mesh_file = NetCDFFile(args.mesh_filename, 'r')
    else:
        mesh_file = None
    extra_dims = parse_extra_dims(args.dimension_list, time_series_file, mesh_file)
                                  
    (extra_dims, cellVars) = setup_dimension_values_and_sort_vars(
            time_series_file, mesh_file, args.variable_list, 
            args.transect_dim_name, extra_dims)
    
    if len(cellVars) == 0:
        print "No variables to extract."
        exit(0)
            
    print " -- Building transects --"
    
    transects = build_transects( mesh_file, transects_file, time_series_file,
                                 args.transect_dim_name, args.min_level_cell_name,
                                 args.max_level_cell_name, use_32bit )


    time_series_file.close()

    summarize_extraction(args.mesh_filename, args.transects_filename, time_indices, cellVars)

    print " -- Extracting cell fields on transects --"

    build_transects_time_series( time_indices, time_file_names, mesh_file,
                                 extra_dims, args.transect_dim_name, transects,
                                 cellVars, args.layer_thickness_name, 
                                 args.z_min_name, args.z_max_name, use_32bit,
                                 args.combine_output, args.append )
    if separate_mesh_file:
        mesh_file.close()

    transects_file.close()

# vim: set expandtab:
