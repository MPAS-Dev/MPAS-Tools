% function moc

% Plot cross-sections of means of MPAS fields.
% 
% This is the main function, where the user can specify data files,
% coordinates and text, then call functions to find sections, load
% data, and plot cross-sections. 
%
% The final product is a set of plots as jpg files, a latex file,
% and a compiled pdf file of the plots, if desired.
%
% Mark Petersen, MPAS-Ocean Team, LANL, January 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify data files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all plots are placed in the f directory.  Comment out if not needed.
unix('mkdir -p f docs');

% The text string [wd '/' sim(i).dir '/' sim(i).netcdf_file ] is the file path,
% where wd is the working directory and dir is the run directory.

wd = '/var/tmp/mpeterse/runs/';

% These files only need to contain a small number of variables.
% You may need to reduce the file size before copying to a local
% machine using:
% ncks -v avgNormalVelocityReconstructMeridional,avgNormalVelocityReconstructZonal, \
% nAverage,latVertex,lonVertex,verticesOnEdge,edgesOnVertex,refLayerThickness,\
% dvEdge,latCell,lonCell,cellsOnCell,nEdgesOnCell \
% file_in.nc file_out.nc

dir='m91';
abc='klmnop';

for j=1:length(abc)
  sim(j).dir = [dir abc(j)];
  sim(j).netcdf_file = ['output_total_avg.nc'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify section coordinates and text
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose Atlantic or global MOC
%region='Atlant' 
 region='global'

% Compute MOC using section
%   [startlat startlon endlat endlon]
if region=='Atlant' 
  sectionText = {'Atlantic MOC'};
  sectionCoord = [-34.5 19.9 -34.5  -55] % '34.5S, South America to Africa -63 to 22
  mocLat = [-34.5:.5:70];
elseif region=='global'
  sectionText = {'Global MOC'};
  mocLat = [-80:.5:85];
else
  fprintf('Incorrect region name')
  return
end
  
% For plotting, only four plots are allowed per row.
% Choose sections above for each page.
% page.name          name of this group of sections 
% sectionID          section numbers for each row of this page
page(1).name = 'NA';
page(1).sectionID = [1:1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify variables to view
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% var_name(nVars)    a cell array with text for each variable to
%                    load or compute.
% var_conv_factor    multiply each variable by this unit conversion.
% var_lims(nVars,3)  contour line definition: min, max, interval 

% Typical variables used for plotting:
% Eulerian velocity from prognostic momentum equation
hor_var_name = {'avgNormalVelocity'};vert_var_name = {'avgVertVelocityTop'};fign=1;
% total transport velocity
%hor_var_name = {'avgNormalTransportVelocity'}; vert_var_name = {'avgVertTransportVelocityTop'};fign=2
% remaining: eddy-induced transport
%hor_var_name = {'avgNormalGMBolusVelocity'}; vert_var_name = {'avgVertGMBolusVelocityTop'};fign=3

var_conv_factor = [1 1 1]; % No conversion here.
if region=='Atlant' 
  contour_lims = [-10:2:10];
elseif region=='global'
  contour_lims = [-40:4:40];
  %contour_lims = [-20:4:-4 -2 2 4:4:20];  % for Bolus velocity MOS
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify latex command
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This matlab script will invoke a latex compiler in order to
% produce a pdf file.  Specify the unix command-line latex
% executable, or 'none' to not compile the latex document.

latex_command = 'latex';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify actions to be taken
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if region=='Atlant' 
  find_edge_sections_flag         = true ;
  plot_edge_sections_flag         = false ;
  compute_transport_flag          = true ;
elseif region=='global'
  find_edge_sections_flag         = true ;
  plot_edge_sections_flag         = true ;
  compute_transport_flag          = false ;
end
load_vertical_velocity_flag     = true ;
compute_moc_flag                = true ;
plot_moc_flag                   = true ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Begin main code.  Normally this does not need to change.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSections = 1;

for iSim = 1:length(sim)

  fprintf(['**** simulation: ' sim(iSim).dir '\n'])
  fid_latex = fopen('temp.tex','w');
  fprintf(fid_latex,['%% file created by plot_mpas_cross_sections, ' date '\n\n']);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  load vertical velocity
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if load_vertical_velocity_flag
    [sim(iSim).avgVertVelocityTop, sim(iSim).botDepth, ...
     sim(iSim).latCell,sim(iSim).lonCell, sim(iSim).areaCell,nVertLevels] = ...
        load_vertical_velocity(wd,sim(iSim).dir,sim(iSim).netcdf_file,vert_var_name);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Find edges that connect beginning and end points of section
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if find_edge_sections_flag 
    [sim(iSim).sectionEdgeIndex, sim(iSim).sectionEdgeSign, sim(iSim).nEdgesInSection, ...
     sim(iSim).latSectionVertex,sim(iSim).lonSectionVertex, ...
     sim(iSim).latVertexDeg,sim(iSim).lonVertexDeg] ...
       = find_edge_sections(wd,sim(iSim).dir,sim(iSim).netcdf_file, ...
       sectionText,sectionCoord);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Plot edge section locations on world map
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if plot_edge_sections_flag
    sub_plot_edge_sections(sim(iSim).dir,sectionCoord, ...
       sim(iSim).latSectionVertex,sim(iSim).lonSectionVertex, ...
       sim(iSim).latVertexDeg,sim(iSim).lonVertexDeg, ...
       sim(iSim).sectionEdgeIndex, sim(iSim).nEdgesInSection,...
       fid_latex);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Load large variables from netcdf file for section
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if compute_transport_flag
    [sim(iSim).sectionData] = load_large_variables_edge ...
       (wd,sim(iSim).dir,sim(iSim).netcdf_file, hor_var_name, var_conv_factor, ...
        sim(iSim).sectionEdgeIndex, sim(iSim).nEdgesInSection);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Compute transport through each section
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if compute_transport_flag
    transport = compute_transport ...
     (wd,sim(iSim).dir,sim(iSim).netcdf_file, hor_var_name,  ...
      sim(iSim).sectionEdgeIndex, sim(iSim).sectionEdgeSign, sim(iSim).nEdgesInSection, ...
      sim(iSim).sectionData,sectionText,sectionAbbreviation);
  else
    transport = zeros(nVertLevels,nSections);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Compute MOC
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if compute_moc_flag
    if region=='Atlant' 
      [sim(iSim).landMask] = land_mask_Atlantic(sim(iSim).latCell,sim(iSim).lonCell);
    elseif region=='global'
      [sim(iSim).landMask] = land_mask_global(sim(iSim).latCell,sim(iSim).lonCell);
    end
    
    [sim(iSim).mocTop] = compute_moc_from_w ...
       (sim(iSim).avgVertVelocityTop, sim(iSim).botDepth, ...
     sim(iSim).latCell,sim(iSim).lonCell, sim(iSim).areaCell,transport,mocLat,sim(iSim).landMask);
  
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Plot MOC
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if plot_moc_flag
    plot_moc(sim(iSim).dir,sectionText,sim(iSim).mocTop,mocLat, ...
	     sim(iSim).botDepth,contour_lims,vert_var_name(1),fign)
  end
  
%  doc_dir = ['docs/' regexprep(sim(iSim).dir,'/','_') '_' ...
%	     sim(iSim).netcdf_file '_dir' ];
%  unix(['mkdir -p ' doc_dir '/f']);
%  unix(['mv f/*jpg ' doc_dir '/f']);

%  filename = [ regexprep(sim(iSim).dir,'/','_') '_' sim(iSim).netcdf_file '.tex'];
%  unix(['cat mpas_sections.head.tex temp.tex > ' doc_dir '/' filename ]);

%  if not(strcmp(latex_command,'none'))
%    fprintf('*** Compiling latex document \n')
%    cd(doc_dir);
%    unix([latex_command ' ' filename]);
%    cd('../..');
%  end
  
end % iSim


