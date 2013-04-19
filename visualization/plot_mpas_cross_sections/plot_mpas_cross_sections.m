function plot_mpas_cross_sections

% Plot cross-sections of MPAS fields.
% 
% This is the main function, where the user can specify data files,
% coordinates and text, then call functions to find sections, load
% data, and plot cross-sections. 
%
% The final product is a set of plots as jpg files, a latex file,
% and a compiled pdf file of the plots, if desired.
%
% Mark Petersen, MPAS-Ocean Team, LANL, May 2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify data files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all plots are placed in the f directory.  Comment out if not needed.
unix('mkdir -p f docs');

% The text string [wd '/' sim(i).dir '/' sim(i).netcdf_file ] is the file path,
% where wd is the working directory and dir is the run directory.

wd = '/var/tmp/mpeterse/runs/todds_runs/P1';

% These files only need to contain a small number of variables.
% You may need to reduce the file size before copying to a local
% machine using:
% ncks -v acc_uReconstructMeridional,acc_uReconstructZonal, \
% nAccumulate,latVertex,lonVertex,verticesOnEdge,edgesOnVertex,hZLevel,\
% dvEdge,latCell,lonCell,cellsOnCell,nEdgesOnCell \
% file_in.nc file_out.nc

sim(1).dir = 'x5.NA.75km_15km';
sim(1).netcdf_file = ['total_avg_o.x5.NA.75km_15km.0029-02-01_00.00.00.nc'];

sim(2).dir = 'x5.NA.75km_15km/leith';
sim(2).netcdf_file = ['total_avg_o.x5.NA.75km_15km.0027-01-01_00.00.00.nc'];

sim(3).dir = 'x1.15km';
sim(3).netcdf_file = ['total_avg_o.x1.15km.0018-11-03_00.00.00.nc'];

sim(4).dir = 'x5.NA.37.5km_7.5km';
sim(4).netcdf_file = ['total_avg_o.x5.NA.37.5km_7.5km.0007-01-01_00.00.00.nc'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify section coordinates and text
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sectionText        a cell array with text describing each section
sectionText = {
'N Atlantic 26N lat',...
'N Atlantic 36N lat',...
'N Atlantic 41N lat',...
'N Atlantic 46N lat',...
'N Atlantic 56N lat',...
'N Atlantic 70W lon',...
'N Atlantic 65W lon',...
'N Atlantic 60W lon',...
'N Atlantic 50W lon',...
'N Atlantic 40W lon',...
'N Atlantic 30W lon',...
'Eq Pacific 140W lon',...
'Eq Pacific 0 lat   ',...
'Drake Pass 65W lon',...
'ACC Tasman 147E lon',...
	      };

% coord(nSections,4)  endpoints of sections, with one section per row as
%                     [startlat startlon endlat endlon]
% Traverse from south to north, and from east to west.
% Then positive velocities are eastward and northward.
coord = [...
  26   -80    26    -15;...   % N Atl Zonal
  36   -76    36    -10;...   % N Atl Zonal
  41   -72    41    -10;...   % N Atl Zonal
  46   -60    46    -10;...   % N Atl Zonal
  56   -60    56    -10;...   % N Atl Zonal
  20   -70    44    -70;...   % N Atl Meridional
  19   -65    44    -65;...   % N Atl Meridional
   8.5 -60    46    -60;...   % N Atl Meridional
   1.8 -50    62    -50;...   % N Atl Meridional
  -3   -40    65    -40;...   % N Atl Meridional
  -5   -30    68.2  -30;...   % N Atl Meridional
  -8  -140     8   -140;...   % Eq Pac Meridional
   0   140     0    -95;...   % Eq Pac Zonal
 -65   -65   -55    -65;...   % Drake
 -67   147   -43.5  147;...   % S Oc Tas
  ];

% plotDepth(nSections) depth to which to plot each section
plotDepth = 4000*ones(1,size(coord,1));
plotDepth(12:13) = 400;

% For plotting, only four plots are allowed per row.
% Choose sections above for each page.
% page.name          name of this group of sections 
% sectionID          section numbers for each row of this page
page(1).name = 'nAtlantic_Zon_p1';
page(1).sectionID = [1:4];
page(2).name = 'nAtlantic_Zon_p2';
page(2).sectionID = [5];

page(3).name = 'nAtlantic_Mer_p1';
page(3).sectionID = [6:9];
page(4).name = 'nAtlantic_Mer_p2';
page(4).sectionID = [10:11];
		    
page(5).name = 'eq_Pac_ACC';
page(5).sectionID = [12:15];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify variables to view
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% var_name(nVars)    a cell array with text for each variable to
%                    load or compute.
% var_conv_factor    multiply each variable by this unit conversion.
% var_lims(nVars,3)  contour line definition: min, max, interval 

var_name = {...
'acc_uReconstructZonal',...
'acc_uReconstructMeridional',...
'ke_acc_uReconstruct'};

var_conv_factor = [100 100 1]; % convert m/s to cm/s for velocities

var_lims = [-10 10 2.5; -10 10 2.5; 0 20 2.5];

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
  %  Find cells that connect beginning and end points of section
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [sim(iSim).sectionCellIndex, sim(iSim).nCellsInSection, sim(iSim).latSection,sim(iSim).lonSection, ...
   depth, sim(iSim).latCellDeg,sim(iSim).lonCellDeg] ...
     = find_cell_sections(wd,sim(iSim).dir,sim(iSim).netcdf_file, ...
     sectionText,coord);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Plot cell section locations on world map
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  sub_plot_cell_sections(sim(iSim).dir,coord, ...
     sim(iSim).latCellDeg,sim(iSim).lonCellDeg, ...
     sim(iSim).sectionCellIndex, sim(iSim).nCellsInSection,fid_latex)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Load large variables from netcdf file
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [sim(iSim).sectionData] = load_large_variables ...
     (wd,sim(iSim).dir,sim(iSim).netcdf_file, var_name, var_conv_factor, ...
      sim(iSim).sectionCellIndex, sim(iSim).nCellsInSection);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Plot data on cross-sections
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for iPage = 1:length(page)

    sub_plot_cross_sections(sim(iSim).dir,sectionText, ...
       page(iPage).name, page(iPage).sectionID,sim(iSim).sectionData,depth,...
       sim(iSim).sectionCellIndex, sim(iSim).nCellsInSection, ...
       sim(iSim).latSection,sim(iSim).lonSection,coord,plotDepth,...
       var_name,var_lims,fid_latex)

  end % iPage
  fprintf(fid_latex,['\n\\end{document}\n\n']);
  fclose(fid_latex);

  doc_dir = ['docs/' regexprep(sim(iSim).dir,'/','_') '_' ...
	     sim(iSim).netcdf_file '_dir' ];
  unix(['mkdir -p ' doc_dir '/f']);
  unix(['mv f/*jpg ' doc_dir '/f']);

  filename = [ regexprep(sim(iSim).dir,'/','_') '_' sim(iSim).netcdf_file '.tex'];
  unix(['cat mpas_sections.head.tex temp.tex > ' doc_dir '/' filename ]);

  if not(strcmp(latex_command,'none'))
    fprintf('*** Compiling latex document \n')
    cd(doc_dir);
    unix([latex_command ' ' filename]);
    cd('../..');
  end
  
end % iSim


