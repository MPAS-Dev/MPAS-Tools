%function cross_section

% Plot cross-sections of MPAS fields.
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

sim(1).dir = 'm91a';
sim(1).netcdf_file = ['output.0010-02-01_00.00.00.nc_transport_vars.nc'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify section coordinates and text
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see exampleSections.m for more example sections.

% sectionText        a cell array with text describing each section
sectionText = {
'ACC  0E lon',...
'ACC 30E lon',...
'ACC 60E lon',...
'ACC 90E lon',...
'Drake Pass 65W lon',...
'ACC Tasman 147E lon',...
	      };

% coord(nSections,4)  endpoints of sections, with one section per row as
%                     [startlat startlon endlat endlon]
% Traverse from south to north, and from east to west.
% Then positive velocities are eastward and northward.
coord = [...
 -67     0   -35    0;...   % S Oc Tas
 -67    30   -35   30;...   % S Oc Tas
 -67    60   -35   60;...   % S Oc Tas
 -67    90   -35   90;...   % S Oc Tas
 -65   -65   -55    -65;...   % Drake
 -67   147   -43.5  147;...   % S Oc Tas
  ];
nSections = size(coord,1);

% number of points to plot for each figure
nPoints = 300;

% plotDepth(nSections) depth to which to plot each section, in m
plotDepth = 5000*ones(1,size(coord,1));

% For plotting, only four plots are allowed per row.
% Choose sections above for each page.
% page.name          name of this group of sections 
% sectionID          section numbers for each row of this page
page(1).name = 'ACC';
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
'avgVelocityZonal',...
'avgVelocityMeridional',...
'keAvgVelocity'};
var_conv_factor = [100 100 1]; % convert m/s to cm/s for velocities
var_lims = [-20 20 2.0; -10 10 1.0; 0 20 2.5];

var_name = {...
'temperature',...
};
var_conv_factor = [1]; 
var_lims = [-1 14 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify actions to be taken
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

find_cell_weights_flag          = true ;
plot_section_locations_flag     = true ;
load_large_variables_flag       = true ;

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
     refMidDepth, refBottomDepth, sim(iSim).maxLevelCellSection] ...
       = find_cell_weights(wd,sim(iSim).dir,sim(iSim).netcdf_file, ...
       sectionText,coord,nPoints);
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
        sim(iSim).cellsOnVertexSection, sim(iSim).cellWeightsSection);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Plot data on cross-sections
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for iPage = 1:length(page)

    sub_plot_cross_sections(sim(iSim).dir,sim(iSim).netcdf_file,sectionText, ...
       page(iPage).name, page(iPage).sectionID,sim(iSim).sectionData,refMidDepth,refBottomDepth,...
       sim(iSim).latSection,sim(iSim).lonSection,sim(iSim).maxLevelCellSection,coord,plotDepth,...
       var_name,var_lims,fid_latex)

  end % iPage
  fprintf(fid_latex,['\n\\end{document}\n\n']);
  fclose(fid_latex);

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


