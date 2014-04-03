function [avgVertVelocityTop, refBottomDepth, latCell,lonCell, areaCell] ...
    = load_vertical_velocity ...
   (wd,dir,netcdf_file)

% load vertical velcoity
%
% Mark Petersen, MPAS-Ocean Team, LANL, May 2012
%
%%%%%%%%%% input arguments %%%%%%%%%
% The text string [wd '/' dir '/' netcdf_file ] is the file path,
% where wd is the working directory and dir is the run directory.
%
%%%%%%%%%% output arguments %%%%%%%%%
% refBottomDepth(nVertLevels)                         depth of center of each layer, for plotting
% vertVelocityTop(nVertLevels,nCells)       vertical velocity at cell center, top

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Read data from file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['** compute u_acc, simulation: ' dir '\n'])

filename = [wd '/' dir '/' netcdf_file ]
ncid = netcdf.open(filename,'nc_nowrite');

avgVertVelocityTop = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'avgVertVelocityTop'));
refBottomDepth = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'refBottomDepth'));
latCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'latCell'))*180/pi;
lonCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lonCell'))*180/pi;
areaCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'areaCell'));

netcdf.close(ncid)

fprintf(['\n'])
