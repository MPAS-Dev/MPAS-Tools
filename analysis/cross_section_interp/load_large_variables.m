function [sectionData] = load_large_variables ...
   (wd,dir,netcdf_file, var_name, var_conv_factor, ...
    cellsOnVertexSection, cellWeightsSection)

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
% cellsOnVertexSection(vertexDegree,nPoints,nSections)  cells neighboring nearest vertex
% cellWeightsSection(vertexDegree,nPoints,nSections)    weights for each cell
%
%%%%%%%%%% output arguments %%%%%%%%%
% sectionData(nVertLevels,nPoints,nSections,nVars)
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

nAccumulate = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'nAccumulate')); 

[dimname,nVertLevels]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nVertLevels'));
[dimname,nCells]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nCells'));

% use only when averaging in matlab:
%[dimname,nTimeSlices]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'Time'));

nVars = length(var_name);
vertexDegree = size(cellsOnVertexSection,1);
nPoints      = size(cellsOnVertexSection,2);
nSections    = size(cellsOnVertexSection,3);

sectionData = zeros(nVertLevels,nPoints,nSections,nVars);

for iVar=1:nVars
  temptext = char(var_name(iVar));
  fprintf(['loading: ' temptext '\n'])
  if (temptext(1:6)=='ke_acc')
    % assume zonal and meridional velocity are in index 1 and 2.
    sectionData(:,:,:,iVar) = sqrt(sectionData(:,:,:,1).^2 + sectionData(:,:,:,2).^2)/2;
  else

    mean_var = netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(var_name(iVar))))*var_conv_factor(iVar); 

    % use only when averaging in matlab:
    %acc_var = netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(var_name(iVar)))); 
    %mean_var = zeros(nVertLevels, nCells);
    %for i=1:nTimeSlices
    %  mean_var = mean_var + nAccumulate(i)*squeeze(acc_var(:,:,i));
    %end
    %mean_var = mean_var/sum(nAccumulate)*var_conv_factor(iVar);
			       
    for iSection = 1:nSections
      for iPt = 1:nPoints
        for k=1:nVertLevels
          sectionData(k,iPt,iSection,iVar) = ...
	      mean_var(k,cellsOnVertexSection(1,iPt,iSection))*cellWeightsSection(1,iPt,iSection)...
	    + mean_var(k,cellsOnVertexSection(2,iPt,iSection))*cellWeightsSection(2,iPt,iSection)...
	    + mean_var(k,cellsOnVertexSection(3,iPt,iSection))*cellWeightsSection(3,iPt,iSection);
        end
      end
    end
  end

end
netcdf.close(ncid)

fprintf('\n')

