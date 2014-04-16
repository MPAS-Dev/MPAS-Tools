function [cellsOnVertexSection, cellWeightsSection, latSection,lonSection, ...
	  refLayerThickness, refMidDepth, refBottomDepth, maxLevelCellSection, sphere_radius] = find_cell_weights ...
   (wd,dir,netcdf_file,sectionText,coord,nLat,nLon)

% This function reads grid data from an MPAS-Ocean grid or restart
% netCDF file, and finds a path of cells that connect the endpoints
% specified in coord.  The path is forced to travel through cells
% that are closest to the line connecting the beginning and end
% cells. 
%
% Mark Petersen, MPAS-Ocean Team, LANL, May 2012
%
%%%%%%%%%% input arguments %%%%%%%%%
% The text string [wd '/' dir '/' netcdf_file ] is the file path,
% where wd is the working directory and dir is the run directory.
% sectionText         a cell array with text describing each section
% coord(nSections,4)  endpoints of sections, with one section per row as
%                     [startlat startlon endlat endlon]
%
%%%%%%%%%% output arguments %%%%%%%%%
% cellsOnVertexSection(vertexDegree,nLat,nLon,nSections)  cells neighboring nearest vertex
% cellWeightsSection(vertexDegree,nLat,nLon,nSections)    weights for each cell
% latSection(nLat,nSections) lat coordinates of each section
% lonSection(nLon,nSections) lon coordinates of each section
% maxLevelCellSection(nLat,nLon,nSections)   min of maxLevelCell of cells surrounding vertex
% refMidDepth(nVertLevels)                         depth of center of each layer, for plotting
% latCellDeg(nCells)                         lat arrays for all cells
% lonCellDeg(nCells)                         lon arrays for all cells  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Read data from file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['** find_cell_sections, simulation: ' dir '\n'])
 
filename = [wd '/' dir '/' netcdf_file ]
ncid = netcdf.open(filename,'nc_nowrite');

xCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'xCell'));
yCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'yCell'));
zCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zCell'));
latVertex = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'latVertex'));
lonVertex = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lonVertex'));
xVertex = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'xVertex'));
yVertex = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'yVertex'));
zVertex = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'zVertex'));
cellsOnVertex = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'cellsOnVertex'));
refLayerThickness = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'refLayerThickness'));
refBottomDepth = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'refBottomDepth'));
maxLevelCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'maxLevelCell'));
sphere_radius = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'sphere_radius');
[dimname,nVertices]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nVertices'));
[dimname,vertexDegree]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'vertexDegree'));
[dimname,nVertLevels]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nVertLevels'));
netcdf.close(ncid)

nSections = size(coord,1);

% Compute depth of center of each layer, for plotting
refMidDepth(1) = refLayerThickness(1)/2;
for i=2:nVertLevels
  refMidDepth(i) = refMidDepth(i-1) + 0.5*(refLayerThickness(i) + refLayerThickness(i-1));
end
% depth(1)=0;  % make top layer plot to surface

latSection = zeros(nLat,nSections);
lonSection = zeros(nLon,nSections);
maxLevelCellSection = zeros(nLat,nLon,nSections);
nearestVertexSection = zeros(nLat,nLon,nSections);
cellsOnVertexSection = zeros(vertexDegree,nLat,nLon,nSections);
cellWeightsSection   = zeros(vertexDegree,nLat,nLon,nSections);
margin=.5;

for iSection=1:nSections
   fprintf('Finding nearest vertices for Section %g \n',iSection)
   latSection(:,iSection) = linspace(coord(iSection,1),coord(iSection,3),nLat);
   lonSection(:,iSection) = linspace(coord(iSection,2),coord(iSection,4),nLon);

   maxLon = (max(lonSection(:,iSection))+margin)*pi/180;
   minLon = (min(lonSection(:,iSection))-margin)*pi/180;

%   maxLat = (max(latSection(:,iSection))+margin)*pi/180;
%   minLat = (min(latSection(:,iSection))-margin)*pi/180;
   
%   vInd = find(latVertex>minLat&latVertex<maxLat ...
%	     & lonVertex>minLon&lonVertex<maxLon);

%     distToVertex = zeros(length(vInd),1);
   
   for iLat=1:nLat
   
   maxLat = (latSection(iLat,iSection)+margin)*pi/180;
   minLat = (latSection(iLat,iSection)-margin)*pi/180;
   
   vInd = find(latVertex>minLat&latVertex<maxLat ...
	     & lonVertex>minLon&lonVertex<maxLon);

     distToVertex = zeros(length(vInd),1);

     for iLon=1:nLon

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %
     %  Find x,y,z coordinates of each section point
     %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     latPt = latSection(iLat,iSection)*pi/180; % lat of point, in radians
     lonPt = lonSection(iLon,iSection)*pi/180; % lon of point, in radians
     zPt = sphere_radius*sin(latPt);
     r   = sphere_radius*cos(latPt);
     xPt = r*cos(lonPt);
     yPt = r*sin(lonPt);

     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %
     %  Find nearest vertex to each section point
     %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     for i=1:length(vInd)
         distToVertex(i) = sqrt(...
             (xPt - xVertex(vInd(i)))^2 ...
           + (yPt - yVertex(vInd(i)))^2 ...
           + (zPt - zVertex(vInd(i)))^2 );
     end

     [mindist,i] = min(distToVertex);
     iVertex = vInd(i);
     
     % Nonexistant cells ("land cells") have a cell index of zero.
     % In this case, throw out this sampling point.
     if min(cellsOnVertex(:,iVertex))==0

       cellsOnVertexSection(:,iLat,iLon,iSection) = 0;
       nearestVertexSection(iLat,iLon,iSection) = 0;
       maxLevelCellSection(iLat,iLon,iSection) = 0;
     
     else % yes, vertex is surrounded by three ocean cells
       
       cellsOnVertexSection(:,iLat,iLon,iSection) = cellsOnVertex(:,iVertex);
       nearestVertexSection(iLat,iLon,iSection) = iVertex;
       maxLevelCellSection(iLat,iLon,iSection) = min(... 
         maxLevelCell(cellsOnVertexSection(:,iLat,iLon,iSection)));
     
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %
       %  Compute area weights for each cell
       %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       vPt = [xPt,yPt,zPt];
       for iCell=1:3
	 vCell(:,iCell) = [...
	     xCell(cellsOnVertexSection(iCell,iLat,iLon,iSection)) ...
	     yCell(cellsOnVertexSection(iCell,iLat,iLon,iSection)) ...
	     zCell(cellsOnVertexSection(iCell,iLat,iLon,iSection))];
       end
       
       for iCell=1:3
	 cellWeightsSection(iCell,iLat,iLon,iSection) = triArea(vPt,...
						  vCell(:,mod(iCell,3)+1),vCell(:,mod(iCell+1,3)+1),sphere_radius);
       end
       cellWeightsSection(:,iLat,iLon,iSection) = ...
	   cellWeightsSection(:,iLat,iLon,iSection) ...
	   /sum(cellWeightsSection(:,iLat,iLon,iSection));

     end % vertex is surrounded by three ocean cells

    end % iLon
   end % iLat

end % iSections


fprintf(['\n'])
