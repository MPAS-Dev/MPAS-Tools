function write_edge_sections_netcdf ...
     (dir, ...
      sectionEdgeIndex, sectionEdgeSign, nEdgesInSection,...
      sectionText,sectionAbbreviation,sectionCoord)

% Write section edge information to the netcdf file
% netcdf_files/your_dir_transport_section_edges.nc
%
% Mark Petersen, MPAS-Ocean Team, LANL, May 2012
%
%%%%%%%%%% input arguments %%%%%%%%%
% dir                string with run directory name
% sectionEdgeIndex(maxNEdgesInSection,nSections) edge index of each section
% sectionEdgeSign(maxNEdgesInSection,nSections)  sign of each
%    section, positive is to right of path direction.
% nEdgesInSection(nSections)                 number of cells in each section
% sectionText        a cell array with text describing each section
% sectionAbbreviation an 8-character title for each section
% sectionCoord(nSections,4)  endpoints of sections, with one section per row as
%                     [startlat startlon endlat endlon]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Write section edge information to a netcdf file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
fprintf(['** Write edge information to file: ' dir '\n'])

nSections = length(nEdgesInSection);
maxNEdgesInSection = max(nEdgesInSection);

dir_name1 =  regexprep(dir,'/','_');
filename = ['netcdf_files/' dir_name1 '_transport_section_edges.nc'];
ncid = netcdf.create(filename,'nc_clobber');

% Define the dimensions of the variable.
dimid_nSections = netcdf.defDim(ncid,'nSections',nSections);
dimid_maxNEdgesInSection = netcdf.defDim(ncid,'maxNEdgesInSection',maxNEdgesInSection);
dimid_latLonPairs = netcdf.defDim(ncid,'latLonPairs',4);
dimid_CharLength8 = netcdf.defDim(ncid,'CharLength8',8);
dimid_CharLength120 = netcdf.defDim(ncid,'CharLength120',120);

% Define a new variable in the file.
sectionEdgeIndex_varID = netcdf.defVar(ncid,'sectionEdgeIndex',...
   'int',[dimid_maxNEdgesInSection dimid_nSections]);

sectionEdgeSign_varID = netcdf.defVar(ncid,'sectionEdgeSign',...
   'int',[dimid_maxNEdgesInSection dimid_nSections]);

nEdgesInSection_varID = netcdf.defVar(ncid,'nEdgesInSection',...
  'int', [dimid_nSections]);

sectionText_varID = netcdf.defVar(ncid,'sectionText',...
  'char',[dimid_CharLength120 dimid_nSections]);
sectionAbbreviation_varID = netcdf.defVar(ncid,'sectionAbbreviation',...
  'char',[dimid_CharLength8 dimid_nSections]);
sectionCoord_varID = netcdf.defVar(ncid,'sectionCoord',...
   'double',[dimid_latLonPairs dimid_nSections]);


% Leave define mode and enter data mode to write data.
netcdf.endDef(ncid)

% Write data to variable.
netcdf.putVar(ncid,sectionEdgeIndex_varID,sectionEdgeIndex);
netcdf.putVar(ncid,sectionEdgeSign_varID,sectionEdgeSign);
netcdf.putVar(ncid,nEdgesInSection_varID,nEdgesInSection);
netcdf.putVar(ncid,sectionText_varID,char(sectionText)');
netcdf.putVar(ncid,sectionAbbreviation_varID,sectionAbbreviation');
netcdf.putVar(ncid,sectionCoord_varID,sectionCoord);

netcdf.close(ncid)

fprintf('\n')

