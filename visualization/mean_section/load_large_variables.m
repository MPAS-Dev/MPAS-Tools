function [sectionData] = load_large_variables ...
   (wd,dir,netcdf_file, var_name, var_conv_factor, ...
    cellsOnVertexSection, cellWeightsSection, maxLevelCellSection,latSection,meanDirection)

% Load large variables from netcdf file
%
% Mark Petersen, MPAS-Ocean Team, LANL, May 2012
%
%%%%%%%%%% input arguments %%%%%%%%%
% The text string [wd '/' dir '/' netcdf_file ] is the file path,
% where wd is the working directory and dir is the run directory.
% var_name(nVars)    a cell array with text for each variable to
%                    load or compute.
% var_conv_factor    multiply each variable by this unit conversion.
% cellsOnVertexSection(vertexDegree,nLat,nLon,nSections)  cells neighboring nearest vertex
% cellWeightsSection(vertexDegree,nLat,nLon,nSections)    weights for each cell
% maxLevelCellSection(nLat,nLon,nSections)   min of maxLevelCell of cells surrounding vertex
% latSection(nLat,nSections) lat coordinates of each section
% meanDirection    Direction to take mean: zonal (z) or meridional (m)
%
%%%%%%%%%% output arguments %%%%%%%%%
% sectionData(nVertLevels,nLat,nSections,nVars)
%   data in each cross-section for each variable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Load large variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
fprintf(['** load_large_variables simulation: ' dir '\n'])

filename = [wd '/' dir '/' netcdf_file ];
ncid = netcdf.open(filename,'nc_nowrite');

nAverage = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'nAverage')); 

[dimname,nVertLevels]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nVertLevels'));
[dimname,nCells]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nCells'));

% use only when averaging in matlab:
%[dimname,nTimeSlices]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'Time'));

nVars = length(var_name);
vertexDegree = size(cellsOnVertexSection,1);
nLat         = size(cellsOnVertexSection,2);
nLon         = size(cellsOnVertexSection,3);
nSections    = size(cellsOnVertexSection,4);

if meanDirection == 'z' % zonal mean
  sectionData = zeros(nVertLevels,nLat,nSections,nVars);
elseif meanDirection == 'm' % meridional mean
  sectionData = zeros(nVertLevels,nLon,nSections,nVars);
end

for iVar=1:nVars
  temptext = char(var_name(iVar));
  fprintf(['loading: ' temptext '\n'])
  if (temptext(1:3)=='keA')
    % assume zonal and meridional velocity are in index 1 and 2.
    sectionData(:,:,:,iVar) = sqrt(sectionData(:,:,:,1).^2 + sectionData(:,:,:,2).^2)/2;
  else

    mean_var = netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(var_name(iVar))))*var_conv_factor(iVar); 

    % use only when averaging in matlab:
    %acc_var = netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(var_name(iVar)))); 
    %mean_var = zeros(nVertLevels, nCells);
    %for i=1:nTimeSlices
    %  mean_var = mean_var + nAverage(i)*squeeze(acc_var(:,:,i));
    %end
    %mean_var = mean_var/sum(nAverage)*var_conv_factor(iVar);
			       
    for iSection = 1:nSections

      if meanDirection == 'z' % zonal mean
	for iLat=1:nLat
	  counter = zeros(1,nVertLevels);

	  for iLon=1:nLon
	    for k=1:maxLevelCellSection(iLat,iLon,iSection)
	      sectionData(k,iLat,iSection,iVar) = sectionData(k,iLat,iSection,iVar) + ...
		  mean_var(k,cellsOnVertexSection(1,iLat,iLon,iSection))*cellWeightsSection(1,iLat,iLon,iSection)...
		  + mean_var(k,cellsOnVertexSection(2,iLat,iLon,iSection))*cellWeightsSection(2,iLat,iLon,iSection)...
		  + mean_var(k,cellsOnVertexSection(3,iLat,iLon,iSection))*cellWeightsSection(3,iLat,iLon,iSection);
	      counter(k) = counter(k) + 1;
	    end
	  end % iLon

	  for k=1:nVertLevels
	    % If no data was recorded at this level, sectionData
	    % is assigned a nan.  Otherwise, take the mean.  
	    if counter(k)~=0
	      sectionData(k,iLat,iSection,iVar) = sectionData(k,iLat,iSection,iVar) /counter(k);
	    else
	      sectionData(k,iLat,iSection,iVar) = nan;
	    end
	  end
	  
	end % iLat
      
      elseif meanDirection == 'm' % meridional mean
        % This average is weighted by area of the lat/lon
	% rectangular cell.  For Meridional means,
	% we weight by cos(lat*pi/180)
	latWeight = cos(latSection(:,iSection)*pi/180);
      
	for iLon=1:nLon
	  latWeightSum = zeros(1,nVertLevels);

	  for iLat=1:nLat
	    for k=1:maxLevelCellSection(iLat,iLon,iSection)
	      sectionData(k,iLon,iSection,iVar) = sectionData(k,iLon,iSection,iVar) + latWeight(iLat)*...
		  mean_var(k,cellsOnVertexSection(1,iLat,iLon,iSection))*cellWeightsSection(1,iLat,iLon,iSection)...
		  + mean_var(k,cellsOnVertexSection(2,iLat,iLon,iSection))*cellWeightsSection(2,iLat,iLon,iSection)...
		  + mean_var(k,cellsOnVertexSection(3,iLat,iLon,iSection))*cellWeightsSection(3,iLat,iLon,iSection);
	      latWeightSum(k) = latWeightSum(k) + latWeight(iLat);
	    end
	  end % iLat

	  for k=1:nVertLevels
	    % If no data was recorded at this level, sectionData
	    % is assigned a nan.  Otherwise, take the mean.  
	    if latWeightSum(k)~=0.0
	      sectionData(k,iLon,iSection,iVar) = sectionData(k,iLon,iSection,iVar) /latWeightSum(k);
	    else
	      sectionData(k,iLon,iSection,iVar) = nan;
	    end
	  end
	  
	end % iLon
      
      end % meanDirection
      
    end % iSection
  end % iVar

end
netcdf.close(ncid)

fprintf('\n')

