function [sectionData] = load_large_variables_edge ...
   (wd,dir,netcdf_file, var_name, var_conv_factor, ...
    sectionEdgeIndex, nEdgesInSection)

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
% sectionEdgeIndex(maxEdges,nSections)       cell index of each section
% nEdgesInSection(nSections)                 number of cells in each section
%
%%%%%%%%%% output arguments %%%%%%%%%
% sectionData(nVertLevels,max(nEdgesInSection),nSections,nVars)
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
[dimname,nEdges]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nEdges'));
[dimname,nTimeSlices]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'Time'));
nVars = length(var_name);
nSections = length(nEdgesInSection);

maxNumberEdges = max(nEdgesInSection);
sectionData = zeros(nVertLevels,maxNumberEdges,nSections,nVars,nTimeSlices);

for iVar=1:nVars
  temptext = char(var_name(iVar));
  fprintf(['loading: ' temptext '\n'])
  
    avg_var = netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(var_name(iVar)))); 
    mean_var = zeros(nVertLevels, nEdges);
    for iTime=1:nTimeSlices
       for iSection = 1:nSections
         for i=1:nEdgesInSection(iSection)
           iEdge = sectionEdgeIndex(i,iSection);
           for k=1:nVertLevels
             sectionData(k,i,iSection,iVar,iTime) = avg_var(k,iEdge,iTime);
           end
         end
       end
    
    end
end

netcdf.close(ncid)

fprintf('\n')

