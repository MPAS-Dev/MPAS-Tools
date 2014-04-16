function [avgVertVelocityTop, refBottomDepth, latCell,lonCell, areaCell, nVertLevels] ...
    = load_vertical_velocity ...
   (wd,dir,netcdf_file,vert_var_name)

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

temp=char(vert_var_name(1));
if temp(1:7)=='avgEddy'
  fprintf('Computing eddy-induced vertical velocity')
  avgVertVelocityTopEulerian = netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(vert_var_name(2))));
  avgVertVelocityTopTransport = netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(vert_var_name(3))));
  avgVertVelocityTop = avgVertVelocityTopEulerian - avgVertVelocityTopTransport;
else
  avgVertVelocityTop = netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(vert_var_name(1))));
end  

refBottomDepth = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'refBottomDepth'));
latCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'latCell'))*180/pi;
lonCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lonCell'))*180/pi;
areaCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'areaCell'));
[dimname,nVertLevels]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nVertLevels'));

netcdf.close(ncid)

fprintf(['\n'])
