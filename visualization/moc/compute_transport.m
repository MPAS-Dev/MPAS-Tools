function transport = compute_transport ...
   (wd,dir,netcdf_file, var_name, ...
    sectionEdgeIndex, sectionEdgeSign, ...
    nEdgesInSection, sectionData,sectionText,sectionAbbreviation)

% Load large variables from netcdf file

% Mark Petersen, MPAS-Ocean Team, LANL, May 2012

%%%%%%%%%% input arguments %%%%%%%%%
% The text string [wd '/' dir '/' netcdf_file ] is the file path,
% where wd is the working directory and dir is the run directory.
% var_name(nVars)    a cell array with text for each variable to
%                    load or compute.
% var_conv_factor    multiply each variable by this unit conversion.
% sectionEdgeIndex(maxEdges,nSections)       cell index of each section
% nEdgesInSection(nSections)                 number of cells in each section
% sectionData(nVertLevels,max(nEdgesInSection),nSections,nVars)
%   data in each cross-section for each variable
% sectionText        a cell array with text describing each section
% sectionAbbreviation an 8-character title for each section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Compute transport through each section
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
fprintf(['** Compute transport: ' dir '\n'])

filename = [wd '/' dir '/' netcdf_file ];
ncid = netcdf.open(filename,'nc_nowrite');

refLayerThickness = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'refLayerThickness'));
dvEdge = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'dvEdge'));
[dimname,nVertLevels]= netcdf.inqDim(ncid,netcdf.inqDimID(ncid,'nVertLevels'));
netcdf.close(ncid)

nSections = length(nEdgesInSection);
maxNEdgesInSection = max(nEdgesInSection);

m3ps_to_Sv = 1e-6; % m^3/sec flux to Sverdrups

% the volume transport
tr = zeros(nVertLevels,maxNEdgesInSection,nSections);
tr_total = zeros(nSections,1);
transport = zeros(nVertLevels,nSections);
header='  ';
data_str='  ';

for iSection = 1:nSections
   for i=1:nEdgesInSection(iSection)
      iEdge = sectionEdgeIndex(i,iSection);
      for k=1:nVertLevels
	 % Compute transport.
	 % I am assuming here that sectionData(:,:,:,1) contains acc_u
	 tr(k,i,iSection) = sectionEdgeSign(i,iSection)...
	     *sectionData(k,i,iSection,1)*dvEdge(iEdge)* ...
	     refLayerThickness(k)*m3ps_to_Sv;
	 tr_total(iSection) = tr_total(iSection) + tr(k,i,iSection);
	 transport(k,iSection) = transport(k,iSection) + tr(k,i,iSection);

      end
   end
   
   % Optional, for plotting the flow across a cross-section.
   % This plots u on edges, so columns oscillate as edges change
   % direction.  The best way to view a cross-section is to use the
   % uMeridional and uZonal at the cell center.
   %figure(iSection+1)
   %imagesc(log(abs(tr(:,1:nEdgesInSection(iSection),iSection))))
   %imagesc(tr(:,1:nEdgesInSection(iSection),iSection))
   %colorbar  

   % note: flow computed in matlab only matches that computed in
   % MPAS-O if they both use refLayerThickness.  To do a verification check,
   % replace the line 
   %     * h_edge(k,iEdge)*m3ps_to_Sv;
   % in mpas_ocn_time_average.F with the line
   %     * refLayerThickness(k,iEdge)*m3ps_to_Sv;

   temptext = char(sectionText(iSection));
   fprintf(['Section %3i, ' temptext(1:22) ' observed flow:' ...
	    temptext(63:75) ' mpas flow: %20.15f Sv\n'],iSection,tr_total(iSection))

%   fprintf(['Section %3i, ' temptext(1:20) 'observed flow:' ...
%	    temptext(63:75) ' mpas flow: %6.1f Sv\n'],iSection,tr_total(iSection))

   header = [header sectionAbbreviation(iSection,:) '   '];
   data_str = [data_str num2str_fixed(tr_total(iSection),'%4.1f',8)...
	      '   '];
end

fprintf(['\n Summary, in Sv: \n' header ' Simulation: \n' data_str '  ' dir '\n'])


fprintf('\n')

