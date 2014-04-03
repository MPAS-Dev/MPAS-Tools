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

abc='abcdef'

for j=1:length(abc)
  sim(j).dir = ['m91' abc(j)];
  sim(j).netcdf_file = ['output_total_avg.nc'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify section coordinates and text
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sectionText        a cell array with text describing each section
sectionText = {
'N Atlantic zonal mean',...
'N Atlantic EUC zonal mean',...
'Eq Pacific 140W lon',...
	      };

% coord(nSections,4)  endpoints of sections, with one section per row as
%                     [startLat startLon endLat endLon]
% Traverse from south to north, and from east to west.
% Then positive velocities are eastward and northward.
coord = [...
  -35 -97 70 -1;...   % N Atlantic zonal mean
  ];

% number of points to plot for each figure
nLat = 100;
nLon = 100;

% Direction to take mean: zonal (z) or meridional (m)
meanDirection = 'z';

% plotDepth(nSections) depth to which to plot each section
plotDepth = 5000*ones(1,size(coord,1));

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

var_conv_factor = [100 100 1]; % convert m/s to cm/s for velocities

%var_lims = [-20 20 2.0; -10 10 1.0; 0 20 2.5];
var_lims = [-1 1 .1; -10 10 1.0; 0 20 2.5];
%var_lims = [-5 5 .5; -10 10 1.0; 0 20 2.5];

%var_lims = [-110 110 10.0; -10 10 1.0; 0 20 2.5];

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

load_vertical_velocity_flag     = true ;
find_edge_sections_flag         = true ;
plot_edge_sections_flag         = false ;
compute_transport_flag          = true ;
compute_moc_flag                = true ;
plot_moc_flag                   = true ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Begin main code.  Normally this does not need to change.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change the coordinate range to be 0 to 360.
coord(:,2) = mod(coord(:,2),360);
coord(:,4) = mod(coord(:,4),360);

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
     sim(iSim).latCell,sim(iSim).lonCell, sim(iSim).areaCell] = ...
        load_vertical_velocity(wd,sim(iSim).dir,sim(iSim).netcdf_file);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Find edges that connect beginning and end points of section
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if find_edge_sections_flag 
    sectionText = {'35S, Atlantic, South America to Africa		                                                 '};
    sectionAbbreviation = [...
    '35S Atlc'];

    %   [startlat startlon endlat endlon]
    sectionCoord = [...
     -34.5 19.9 -34.5  -55] % '35S, South America to Africa -63 to 22
    var_name = {'avgNormalVelocity'};
    var_conv_factor = [1 1 1]; % No conversion here.

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
       (wd,sim(iSim).dir,sim(iSim).netcdf_file, var_name, var_conv_factor, ...
        sim(iSim).sectionEdgeIndex, sim(iSim).nEdgesInSection);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Compute transport through each section
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if compute_transport_flag
    transport = compute_transport ...
     (wd,sim(iSim).dir,sim(iSim).netcdf_file, var_name,  ...
      sim(iSim).sectionEdgeIndex, sim(iSim).sectionEdgeSign, sim(iSim).nEdgesInSection, ...
      sim(iSim).sectionData,sectionText,sectionAbbreviation);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Compute MOC
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if compute_moc_flag
    mocLat = [-34.5 -34:.5:70];
    [sim(iSim).landMask] = land_mask_Atlantic(sim(iSim).latCell,sim(iSim).lonCell);
    [sim(iSim).mocTop] = compute_moc_from_w ...
       (sim(iSim).avgVertVelocityTop, sim(iSim).botDepth, ...
     sim(iSim).latCell,sim(iSim).lonCell, sim(iSim).areaCell,transport,mocLat,sim(iSim).landMask);
  
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Plot MOC
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  var_lims_moc = [-50 50 10];
  sectionText = {
  'moc streamfunction',...
     };
  if plot_moc_flag
    plot_moc(sim(iSim).dir,sectionText,sim(iSim).mocTop,mocLat, ...
	     sim(iSim).botDepth)
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


