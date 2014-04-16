% function mean_section

% Plot cross-sections of means of MPAS fields.
% 
% This is the main function, where the user can specify data files,
% coordinates and text, then call functions to find sections, load
% data, and plot cross-sections. 
%
% The final product is a set of plots as jpg files, a latex file,
% and a compiled pdf file of the plots, if desired.
%
% Mark Petersen, MPAS-Ocean Team, LANL, March 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify data files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all plots are placed in the f directory.  Comment out if not needed.
unix('mkdir -p f docs');

% The text string [wd '/' sim(i).dir '/' sim(i).netcdf_file ] is the file path,
% where wd is the working directory and dir is the run directory.

wd = '/var/tmp/mpeterse/runs';

% These files only need to contain a small number of variables.
% You may need to reduce the file size before copying to a local
% machine using, and also add the variables you want to plot.
%    ncks -v nAverage,latVertex,lonVertex,verticesOnEdge,edgesOnVertex,refLayerThickness,dvEdge,latCell,lonCell,refBottomDepth,areaCell,xCell,yCell,zCell,xVertex,yVertex,zVertex,cellsOnVertex,maxLevelCell \
% file_in.nc file_out.nc

abc='o'

for j=1:length(abc)
  sim(j).dir = ['m91' abc(j)];
  sim(j).netcdf_file = ['output_total_avg.nc'];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify section coordinates and text
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see exampleSections.m for more example sections.

% sectionText        a cell array with text describing each section
sectionText = {
'S Ocean Atlantic zonal & 10 yr mean',...
'S Ocean Indian zonal & 10 yr mean',...
'S Ocean Pacific zonal & 10 yr mean',...
	      };

% coord(nSections,4)  endpoints of sections, with one section per row as
%                     [startLat startLon endLat endLon]
% Traverse from south to north, and from east to west.
% Then positive velocities are eastward and northward.

coord = [...
 -65   -45   -35    0;...    % S Ocean Atlantic
 -65    30   -35  120;...    % S Ocean Indian  
 -65  -180   -35  -90;...    % S Ocean Pacific 
  ];

nSections = size(coord,1);

% number of points to plot for each figure
nLat =  100;
nLon =  100;

% Direction to take mean: zonal (z) or meridional (m)
meanDirection = 'z';

% plotDepth(nSections) depth to which to plot each section, in m
plotDepth = 5000*ones(1,size(coord,1));

% For plotting, only four plots are allowed per row.
% Choose sections above for each page.
% page.name          name of this group of sections 
% sectionID          section numbers for each row of this page
page(1).name = 'SO';
page(1).sectionID = [1:nSections];

% coord range may need alteration to match lonVertex:
% If lonVertex is between 0 and 2*pi, ensure the coordinate range is 0 to 360.
%coord(:,2) = mod(coord(:,2),360);
%coord(:,4) = mod(coord(:,4),360);
% If lonVertex is between -pi and pi, ensure the coordinate range is -180 to 180.
coord(:,2) = mod(coord(:,2)+180,360)-180;
coord(:,4) = mod(coord(:,4)+180,360)-180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify variables to view
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see exampleSections.m for more example variables

% var_name(nVars)    a cell array with text for each variable to
%                    load or compute.
% var_conv_factor    multiply each variable by this unit conversion.
% var_lims(nVars,3)  contour line definition: min, max, interval 

var_name = {...
'avgMeridionalVelocity',...
'avgZonalVelocity',...
'keAvgVelocity',...
'temperature',...
'salinity',...
'potentialDensity',...
    }

var_name = {...
'avgZonalGMBolusVelocity',...
'avgMeridionalGMBolusVelocity',...
'zonalRelativeSlope',...
'meridionalRelativeSlope',...
   }


var_conv_factor = [100 100 1 1 1 1]; % convert m/s to cm/s for velocities

var_lims = [-1.2 1.2 .1; -12 12 1.0; 0 4 0.5; -1 14 1; 34 35.2 .1; 1025.8 1028 .1];

var_lims = [-0.5 0.5 .05; -0.5 0.5 .05; 1e-6*[-20 20 2; -140 140 10]];
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify actions to be taken
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

find_cell_weights_flag          = false ;
plot_section_locations_flag     = false ;
load_large_variables_flag       = true ;
plot_cross_sections_flag        = true ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Begin main code.  Normally this does not need to change.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSim = 1:length(sim)

  fprintf(['**** simulation: ' sim(iSim).dir '\n'])
  fid_latex = fopen('temp.tex','w');
  fprintf(fid_latex,['%% file created by plot_mpas_cross_sections, ' date '\n\n']);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Find cells that connect beginning and end points of section
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if find_cell_weights_flag
    [sim(iSim).cellsOnVertexSection, sim(iSim).cellWeightsSection, ...
     sim(iSim).latSection,sim(iSim).lonSection, ...
     refLayerThickness, refMidDepth, refBottomDepth, sim(iSim).maxLevelCellSection, sphere_radius] ...
       = find_cell_weights(wd,sim(iSim).dir,sim(iSim).netcdf_file, ...
       sectionText,coord,nLat,nLon);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Plot cell section locations on world map
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if plot_section_locations_flag
    sub_plot_section_locations(sim(iSim).dir,coord, ...
       sim(iSim).latSection,sim(iSim).lonSection,fid_latex)
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Load large variables from netcdf file
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if load_large_variables_flag
    [sim(iSim).sectionData] = load_large_variables ...
       (wd,sim(iSim).dir,sim(iSim).netcdf_file, var_name, var_conv_factor, ...
        sim(iSim).cellsOnVertexSection, sim(iSim).cellWeightsSection, ...
	sim(iSim).maxLevelCellSection,sim(iSim).latSection,meanDirection);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Plot data on cross-sections
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if plot_cross_sections_flag

    for iPage = 1:length(page)

    sub_plot_cross_sections(sim(iSim).dir,sim(iSim).netcdf_file,sectionText, ...
       page(iPage).name, page(iPage).sectionID,sim(iSim).sectionData,refMidDepth,refBottomDepth,...
       sim(iSim).latSection,sim(iSim).lonSection,sim(iSim).maxLevelCellSection,coord,plotDepth,...
       var_name,var_lims,meanDirection,fid_latex)

    end % iPage
    fprintf(fid_latex,['\n\\end{document}\n\n']);
    fclose(fid_latex);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Latex Compilation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This matlab script will invoke a latex compiler in order to
% produce a pdf file.  Specify the unix command-line latex
% executable, or 'none' to not compile the latex document.

% latex_command = 'latex';

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


