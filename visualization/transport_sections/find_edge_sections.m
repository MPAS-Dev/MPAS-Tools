function [sectionEdgeIndex, sectionEdgeSign, nEdgesInSection, ...
	  latSectionVertex,lonSectionVertex, ...
	  latVertexDeg,lonVertexDeg] = find_edge_sections ...
   (wd,dir,netcdf_file,sectionText,sectionCoord)

% This function reads grid data from an MPAS-Ocean grid or restart
% netCDF file, and finds a path of edges that connect the endpoints
% specified in sectionCoord.  The path is forced to travel through edges
% that are closest to the line connecting the beginning and end
% edges. 
%
% Mark Petersen, MPAS-Ocean Team, LANL, May 2012
%
%%%%%%%%%% input arguments %%%%%%%%%
% The text string [wd '/' dir '/' netcdf_file ] is the file path,
% where wd is the working directory and dir is the run directory.
% sectionText         a cell array with text describing each section
% sectionCoord(nSections,4)  endpoints of sections, with one section per row as
%                     [startlat startlon endlat endlon]
%
%%%%%%%%%% output arguments %%%%%%%%%
% sectionEdgeIndex(maxNEdgesInSection,nSections) edge index of each section
% sectionEdgeSign(maxNEdgesInSection,nSections)  sign of each
%    section, positive is to right of path direction.
% nEdgesInSection(nSections)                     number of edges in each section
% latSectionVertex(maxNEdgesInSection,nSections) lat coordinates of each section
% lonSectionVertex(maxNEdgesInSection,nSections) lon coordinates of each section
% latVertexDeg(nEdges)                           lat arrays for all edges
% lonVertexDeg(nEdges)                           lon arrays for all edges  

%%%%%%%%%% parameters internal to this function %%%%%%%%%%

% maxEdges specifies the maximum number of Edges attempted along
% the path to the end-edge before stopping with a warning.
maxEdges = 1500;

% Make sure sectionCoord traverse from south to north, and from east to west.
%   [startlat startlon endlat endlon]
nSections = size(sectionCoord,1);
for j=1:nSections
  latChange = sectionCoord(j,3) - sectionCoord(j,1);
  lonChange = sectionCoord(j,4) - sectionCoord(j,2);
  if abs(lonChange)>abs(latChange) % zonal section
    if lonChange>0
      fprintf(['Warning: Zonal sections should go from east to west.  ' ...
         'For section %g start and end longitudes are %g, %g \n'], ...
	  j,sectionCoord(j,2),sectionCoord(j,4))
    end
  else
    if latChange<0
      fprintf(['Warning: Meridional sections should go from south to north.  ' ...
         'For section %g start and end latitudes are %g, %g \n'], ...
	  j,sectionCoord(j,1),sectionCoord(j,3))
    end
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Read edge and edge data from grid file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['** find_edge_sections, simulation: ' dir '\n'])
 
filename = [wd '/' dir '/' netcdf_file ];
ncid = netcdf.open(filename,'nc_nowrite');

latVertex = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'latVertex'));
lonVertex = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lonVertex'));
verticesOnEdge = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'verticesOnEdge'));
edgesOnVertex = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'edgesOnVertex'));
[dimname,nEdges]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nEdges'));
[dimname,nVertices]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nVertices'));
[dimname,vertexDegree]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'vertexDegree'));
netcdf.close(ncid)

% Grid variables should be:
% lat varies from -pi/2:pi/2
% lon varies from 0:2*pi
if (min(lonVertex)<-1e-8)
  lonVertex = mod(lonVertex,2*pi);
end
% convert to degrees for plotting:
latVertexDeg = latVertex*180/pi;
lonVertexDeg = lonVertex*180/pi;

sectionVertexIndex = zeros(maxEdges,nSections);
sectionEdgeIndex = zeros(maxEdges,nSections);
sectionEdgeSign = zeros(maxEdges,nSections);
nEdgesInSection = zeros(1,nSections);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Find edges that connect beginning and ending points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSection=1:nSections
   latCoord = [sectionCoord(iSection,1) sectionCoord(iSection,3)]/180*pi;
   lonCoord = [sectionCoord(iSection,2) sectionCoord(iSection,4)]/180*pi;

   % Find vertex closest to start and end coordinates.
   % The seed vertex array simply stores start and end index.
   minDist = 1e10*ones(1,2);
   seedVertexIndex = zeros(1,2);
   for iVertex = 1:nVertices
     for i=1:2
       dist = sqrt( ...
	  (lonCoord(i) - lonVertex(iVertex))^2 ...
	+ (latCoord(i) - latVertex(iVertex))^2);
       if (dist<minDist(i))
	 minDist(i) = dist;
	 seedVertexIndex(i) = iVertex;
%	 fprintf([' iVertex %i, i %i, dist %g lonVertex %g latVertex' ...
%		  ' %g \n'],iVertex,i,dist,latVertex(iVertex),lonVertex(iVertex))
       end
     end
   end

   % Rename first vertex on route:
   lonBeginVertex = lonVertex(seedVertexIndex(1));
   latBeginVertex = latVertex(seedVertexIndex(1));
   beginVertexIndex = seedVertexIndex(1);

   % Assign first vertex on route to vertex index array:
   sectionVertexIndex(1,iSection) = beginVertexIndex;

   % Rename last vertex on route:
   lonEndVertex = lonVertex(seedVertexIndex(2));
   latEndVertex = latVertex(seedVertexIndex(2));
   endVertexIndex = seedVertexIndex(2);
   
   % Assign distance from beginVertex to endVertex.  This is the
   % distance to 'beat' by candidate neighbor vertexs.
   distLastVertex = sqrt( ...
       (lonBeginVertex - lonEndVertex)^2 ...
     + (latBeginVertex - latEndVertex)^2 );

   % loop through edges in search for section
   for i=1:maxEdges
     additionalEdgeFound = 0;
     minDistToLine = 1e10;

     for j = 1:vertexDegree
       iEdge = edgesOnVertex(j,sectionVertexIndex(i,iSection));
       if (iEdge>0)
	  % Find the vertex on the other side of iEdge
	  if (verticesOnEdge(1,iEdge)==sectionVertexIndex(i,iSection))
	     iVertex = verticesOnEdge(2,iEdge);
	     % Going from vertex 1 to vertex 2.  Leave positive.
	     edgeSign = 1;
	  else
	     iVertex = verticesOnEdge(1,iEdge);
	     % Going from vertex 2 to vertex 1.  Make negative.
	     edgeSign = -1;
	  end	

          % I am using lat/lon Cartesian distance.
	  % This is distance to the final vertex location.
	  dist = sqrt( ...
	     (lonVertex(iVertex) - lonVertex(endVertexIndex))^2 ...
	   + (latVertex(iVertex) - latVertex(endVertexIndex))^2 );

%fprintf('%6i %6i %8.4f %8.4f h1=plot(%g,%g);  h2=plot(%g,%g); \n',...
%i,j,dist,distLastVertex,...
%	lonVertex(iVertex)*180/pi,latVertex(iVertex)*180/pi,...
%	lonVertex(endVertexIndex)*180/pi,latVertex(endVertexIndex)*180/pi)
          % check if this vertex is closer to the end vertex than the
          % last vertex.  If so, it is a candidate, and we can continue.
	  if (dist<distLastVertex)
	    
	    % distance formula from http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
	    % distToLine = abs((x2-x1)*(y1-y0) - (x1-x0)*(y2-y1)) ...
	    %   / sqrt((x2-x1)^2 + (y2-y1)^2);
    
	    distToLine = abs(...
                 (lonEndVertex  -lonBeginVertex  )*(latBeginVertex-latVertex(iVertex)) ...
               - (lonBeginVertex-lonVertex(iVertex))*(latEndVertex  -latBeginVertex) ) ...
	       / sqrt((lonEndVertex-lonBeginVertex)^2 + (latEndVertex-latBeginVertex)^2);

            % Find the vertex that is closest to the line connecting
            % beginning and ending vertexs.
            if (distToLine<minDistToLine)
               additionalEdgeFound = 1;
               distChosenVertex = dist;
               minDistToLine = distToLine;
               minDistVertexIndex = iVertex;
               minDistEdgeIndex = iEdge;
               minDistEdgeSign = edgeSign;
	    end
	    
	  end
       end
       
     end  % nEdgesOnEdge

     distLastVertex = distChosenVertex;

     sectionVertexIndex(i+1,iSection) = minDistVertexIndex;
     sectionEdgeIndex  (i,iSection)   = minDistEdgeIndex;
     sectionEdgeSign   (i,iSection)   = minDistEdgeSign ;
     nEdgesInSection(iSection) = i;

     if (minDistVertexIndex==endVertexIndex)
        break
     end

   end % maxEdges

   temptext = char(sectionText(iSection));
%   fprintf(['Section %3i, ' temptext(1:49) ' '],iSection)
   if (minDistVertexIndex==endVertexIndex)
%     fprintf(' complete with %g edges.\n',nEdgesInSection(iSection))
   else
     fprintf(['WARNING: Did not complete path to end' ...
	      ' vertex with %g maxEdges. \n'],maxEdges)
   end

end % iSections

maxNEdgesInSection = max(nEdgesInSection);

% reduce size of sectionEdgeIndex array
sectionEdgeIndex = sectionEdgeIndex(1:maxNEdgesInSection,:);
sectionEdgeSign = sectionEdgeSign(1:maxNEdgesInSection,:);

% assign lat/lon on section
latSection = zeros(maxNEdgesInSection+1,nSections);
lonSection = zeros(maxNEdgesInSection+1,nSections);
for iSection = 1:nSections
  for i=1:nEdgesInSection(iSection)+1
    iVertex = sectionVertexIndex(i,iSection);
    latSectionVertex(i,iSection) = latVertexDeg(iVertex);
    lonSectionVertex(i,iSection) = lonVertexDeg(iVertex);
  end
end

fprintf(['\n'])
