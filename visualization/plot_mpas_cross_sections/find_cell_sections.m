function [sectionCellIndex, nCellsInSection, latSection,lonSection, ...
	  depth, latCellDeg,lonCellDeg] = find_cell_sections ...
   (wd,dir,netcdf_file,sectionText,coord)

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
% sectionCellIndex(maxCells,nSections)       cell index of each section
% nCellsInSection(nSections)                 number of cells in each section
% latSection(max(nCellsInSection),nSections) lat coordinates of each section
% lonSection(max(nCellsInSection),nSections) lon coordinates of each section
% depth(nVertLevels)                         depth of center of each layer, for plotting
% latCellDeg(nCells)                         lat arrays for all cells
% lonCellDeg(nCells)                         lon arrays for all cells  

%%%%%%%%%% parameters internal to this function %%%%%%%%%%

% maxCells specifies the maximum number of Cells attempted along
% the path to the end-cell before stopping with a warning.
maxCells = 1500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Read cell and edge data from grid file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['** find_cell_sections, simulation: ' dir '\n'])
 
filename = [wd '/' dir '/' netcdf_file ];
ncid = netcdf.open(filename,'nc_nowrite');

latCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'latCell'));
lonCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lonCell'));
cellsOnCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'cellsOnCell'));
nEdgesOnCell = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'nEdgesOnCell'));
hZLevel = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'hZLevel'));
[dimname,nCells]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nCells'));
[dimname,nVertLevels]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nVertLevels'));
netcdf.close(ncid)

nSections = size(coord,1);

% Compute depth of center of each layer, for plotting
depth(1) = hZLevel(1)/2;
for i=2:nVertLevels
  depth(i) = depth(i-1) + 0.5*(hZLevel(i) + hZLevel(i-1));
end

% The native grid variables use:
% lat varies from -pi/2:pi/2
% lon varies from -pi:pi
% convert to degrees for plotting:
latCellDeg = latCell*180/pi;
lonCellDeg = lonCell*180/pi;

sectionCellIndex = zeros(maxCells,nSections);
nCellsInSection = zeros(1,nSections);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Find cells that connect beginning and ending points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSection=1:nSections
   latCoord = [coord(iSection,1) coord(iSection,3)]/180*pi;
   lonCoord = [coord(iSection,2) coord(iSection,4)]/180*pi;

   % Find cell closest to start and end coordinates.
   % The seed cell array simply stores start and end index.
   minDist = 1e10*ones(1,2);
   seedCellIndex = zeros(2,nSections);
   for iCell = 1:nCells
     for i=1:2
       dist = sqrt( ...
	  (lonCoord(i) - lonCell(iCell))^2 ...
	+ (latCoord(i) - latCell(iCell))^2);
       if (dist<minDist(i))
	 minDist(i) = dist;
	 seedCellIndex(i) = iCell;
       end
     end
   end

   % Rename first cell on route:
   lonBeginCell = lonCell(seedCellIndex(1));
   latBeginCell = latCell(seedCellIndex(1));
   beginCellIndex = seedCellIndex(1);

   % Assign first cell on route to cell index array:
   sectionCellIndex(1,iSection) = beginCellIndex;

   % Rename last cell on route:
   lonEndCell = lonCell(seedCellIndex(2));
   latEndCell = latCell(seedCellIndex(2));
   endCellIndex = seedCellIndex(2);
   
   % Assign distance from beginCell to endCell.  This is the
   % distance to 'beat' by candidate neighbor cells.
   distLastCell = sqrt( ...
       (lonBeginCell - lonEndCell)^2 ...
     + (latBeginCell - latEndCell)^2 );

   % loop through cells in search for section
   for i=1:maxCells
     additionalCellFound = 0;
     minDistToLine = 1e10;

     for j = 1:nEdgesOnCell(sectionCellIndex(i,iSection))
       iCell = cellsOnCell(j,sectionCellIndex(i,iSection));

       if (iCell<nCells+1 && iCell>0)
	  % I am using lat/lon Cartesian distance.
	  % This is distance to the final cell location.
	  dist = sqrt( ...
	     (lonCell(iCell) - lonEndCell)^2 ...
	   + (latCell(iCell) - latEndCell)^2 );

	  % check if this cell is closer to the end cell than the
          % last cell.  If so, it is a candidate, and we can continue.
	  if (dist<distLastCell)
	    
	    % distance formula from http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
	    % distToLine = abs((x2-x1)*(y1-y0) - (x1-x0)*(y2-y1)) ...
	    %   / sqrt((x2-x1)^2 + (y2-y1)^2);
    
	    distToLine = abs(...
                 (lonEndCell  -lonBeginCell  )*(latBeginCell-latCell(iCell)) ...
               - (lonBeginCell-lonCell(iCell))*(latEndCell  -latBeginCell) ) ...
	       / sqrt((lonEndCell-lonBeginCell)^2 + (latEndCell-latBeginCell)^2);

            % Find the cell that is closest to the line connecting
            % beginning and ending cells.
            if (distToLine<minDistToLine)
              additionalCellFound = 1;
	      distChosenCell = dist;
	      minDistToLine = distToLine;
	      minDistCellIndex = iCell;
	    end
	    
	  end
       end
       
     end  % nEdgesOnCell

     distLastCell = distChosenCell;

     sectionCellIndex(i+1,iSection) = minDistCellIndex;
     nCellsInSection(iSection) = i+1;

     if (minDistCellIndex==endCellIndex)
        break
     end

   end % maxCells

   temptext = char(sectionText(iSection));
   fprintf(['Section %3i, ' temptext ' '],iSection)
   if (minDistCellIndex==endCellIndex)
     fprintf(' complete with %g cells.\n',nCellsInSection(iSection))
   else
     fprintf(['WARNING: Did not complete path to end' ...
	      ' cell with %g maxCells. \n'],maxCells)
   end

end % iSections

maxNumberCells = max(nCellsInSection);

% reduce size of sectionCellIndex array
sectionCellIndex = sectionCellIndex(1:maxNumberCells,:);

% assign lat/lon on section
latSection = zeros(maxNumberCells,nSections);
lonSection = zeros(maxNumberCells,nSections);
for iSection = 1:nSections
  for i=1:nCellsInSection(iSection)
    iCell = sectionCellIndex(i,iSection);
    latSection(i,iSection) = latCellDeg(iCell);
    lonSection(i,iSection) = lonCellDeg(iCell);
  end
end

fprintf(['\n'])
