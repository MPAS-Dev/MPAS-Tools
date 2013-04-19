function [sectionData] = load_large_variables_avg ...
   (wd,dir,netcdf_file, var_name, var_conv_factor, ...
    sectionCellIndex, nCellsInSection)

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
% sectionCellIndex(maxCells,nSections)       cell index of each section
% nCellsInSection(nSections)                 number of cells in each section
%
%%%%%%%%%% output arguments %%%%%%%%%
% sectionData(nVertLevels,max(nCellsInSection),nSections,nVars)
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
[dimname,nTimeSlices]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'Time'));
nVars = length(var_name);
nSections = length(nCellsInSection);

maxNumberCells = max(nCellsInSection);
sectionData = zeros(nVertLevels,maxNumberCells,nSections,nVars);

for iVar=1:nVars
  temptext = char(var_name(iVar));
  fprintf(['loading: ' temptext '\n'])
  if (temptext(1:6)=='ke_acc')
    % assume zonal and meridional velocity are in index 1 and 2.
    sectionData(:,:,:,iVar) = sqrt(sectionData(:,:,:,1).^2 + sectionData(:,:,:,2).^2)/2;
  else
  
    acc_var = netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(var_name(iVar)))); 
    mean_var = zeros(nVertLevels, nCells);
    for i=1:nTimeSlices
      mean_var = mean_var + nAccumulate(i)*squeeze(acc_var(:,:,i));
    end
    mean_var = mean_var/sum(nAccumulate)*var_conv_factor(iVar);
			       
    for iSection = 1:nSections
      for i=1:nCellsInSection(iSection)
        iCell = sectionCellIndex(i,iSection);
        for k=1:nVertLevels
          sectionData(k,i,iSection,iVar) = mean_var(k,iCell);
        end
      end
    end
  end

end
netcdf.close(ncid)

fprintf('\n')

