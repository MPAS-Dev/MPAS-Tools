%function transport_sections

% Specify data files, coordinates and text, then call functions
% to find edge sections, load data, and compute transport through
% each section.
%
% This script produces new netcdf files in the subdirectory
% netcdf_files which can then be merged with grid.nc or restart.nc
% files to collect transport data in MPAS-Ocean
%
% To merge the new *_section_edge_data.nc with an existing grid or
% restart file, use:
% ncks -A -v sectionEdgeIndex,sectionEdgeSign,nEdgesInSection,\
% sectionText,sectionAbbreviation,sectionCoord \
% your_file_section_edge_data.nc your_restart_file.nc
%
% Mark Petersen, MPAS-Ocean Team, LANL, March 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify data files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all plots are placed in the f directory.  Comment out if not needed.
unix('mkdir -p f netcdf_files docs text_files');

% The text string [wd '/' sim(i).dir '/' sim(i).netcdf_file ] is the file path,
% where wd is the working directory and dir is the run directory.

wd = '/var/tmp/mpeterse/runs/';
dir='m91';
abc = 'klmnop';

for letter=1:length(abc)

% These files only need to contain a small number of variables.
% You may need to reduce the file size before copying to a local
% machine using:
% ncks -v avgNormalVelocity,avgNormalTransportVelocity,nAverage,latVertex,lonVertex,verticesOnEdge,edgesOnVertex,refLayerThickness,dvEdge \
% file_in.nc file_out.nc

clear sim
for j=1:3
  sim(j).dir=[dir abc(letter)];
  sim(j).netcdf_file = ['output.00' num2str_fixed0(16+j,'%g',2) '-02-01_00.00.00.nc_transport_vars.nc'];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify section coordinates and text
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sectionText        a cell array with text describing each section
sectionText = {
'Drake Passage, S Ocean -56  to -63 lat,  68W lon, section A21,  140+/- 6 Sv in Ganachaud_Wunsch00n and Ganachaud99thesis',...
'Tasmania-Ant, S Ocean  -44  to -66 lat, 140E lon, section P12,  157+/-10 Sv in Ganachaud_Wunsch00n and Ganachaud99thesis',...
'Africa-Ant, S Ocean    -31.3to -70 lat,  30E lon, section  I6,           Sv in Ganachaud99thesis                        ',...
'Antilles Inflow, Carib.                                       -18.4+/-4.7Sv in Johns_ea02dsr                            '...
'Mona Passage, Caribbian                                        -2.6+/-1.2Sv in Johns_ea02dsr                            '...
'Windward Passage, Carib                                        -7.0      Sv in Nelepo_ea76sr, Roemmich81jgr             '...
'Florida-Cuba, Caribbian                                        31.5+/-1.5Sv in Johns_ea02dsr, 32.3+/-3.2Sv Larsen92rslpt'...
'Florida-Bahamas, Carib. 27 lat, -80  to -78.8lon,              31.5+/-1.5Sv in Johns_ea02dsr, 32.3+/-3.2Sv Larsen92rslpt'...
'Indonesian Throughflow, -9  to -18 lat, 116E lon, section J89,  -16+/- 5 Sv in Ganachaud_Wunsch00n and Ganachaud99thesis',...
'Agulhas                                                         -70+/-20 Sv in Bryden_Beal01dsr                         ',...
'Mozambique Channel,    -25 lat,  35  to  44E lon, section I4 ,  -14+/- 6 Sv in Ganachaud_Wunsch00n and Ganachaud99thesis',...
'Bering Strait, Arctic                                          0.83+/-0.66Sv in Roach_ea95jgr                           '...
    };
%'Lancaster Sound, Arctic                                        0.67+/-0.3Sv in Maltrud_McLean05om                       '...
%'Fram Strait, Arctic                                           -4.2+/-2.3Sv in Fahrbach_ea01pr                           '...
%'Robeson Channel, Arctic                                       -0.75+/-0.2Sv in Maltrud_McLean05om                       '...

% sectionAbbreviation an 8-character title for each section
sectionAbbreviation = [...
    'Drake Pa';...
    'Tasm-Ant';...
    'Afri-Ant';...
    'Antilles';...
    'Mona Pas';...
    'Wind Pas';...
    'FL-Cuba ';...
    'FL-Baham';...
    'Ind Thru';...
    'Agulhas ';...
    'Mozam Ch';...
    'Bering  ';...
];
%    'Lancastr';...
%   'Fram    ';...
%   'Robeson ';...

% sectionCoord(nSections,4)  endpoints of sections, with one section per row as
%   [startlat startlon endlat endlon]
% Traverse from south to north, and from east to west.
% Then positive velocities are eastward and northward.
sectionCoord = [...
 -64.5  -64    -55    -65.3;... % Drake
 -67    140    -43.5  147  ;... % Tasm-Ant
 -70.0   30    -31.3   30  ;... % Afri-Ant
  10.7  -63.2   18.0  -65.9;... % Antilles  
  18.4  -67.2   18.4  -68.5;... % Mona Passage
  19.8  -73.4   20.1  -74.3;... % Windward Passage
  23.1  -81.0   25.15 -81.0;... % Florida-Cuba
  26.52 -78.78  26.7  -80.1;... % Florida-Bahamas
 -21    116.0   -8.8  116  ;... % Ind Thru
 -32.4   32.0  -31.0   30.2;... % Agulhas
 -25     44.0  -25.0   34  ;... % Mozam Ch
  65.8 -167.7   66.1 -169.7;... % Bering St
  ];
%  73.7  -80.6   74.6  -81.0;... % Lancaster Sound- was not able to
%  get this to connect for all resolutions
% 79.7   10.7   79.7  -17.7;... % Fram St - crosses 0 lon.  This is not in code yet.
% 81.0  -63.5   82.0  -63.5;... % Robeson Ch - was not able to get this to connect 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% var_name(nVars)    a cell array with text for each variable to
%                    load or compute.
% var_conv_factor    multiply each variable by this unit conversion.
% var_lims(nVars,3)  contour line definition: min, max, interval 

% Eulerian velocity from prognostic momentum equation
var_name = {'avgNormalVelocity'};
% total transport velocity
%var_name = {'avgNormalTransportVelocity'}

var_conv_factor = [1 1 1]; % No conversion here.

var_lims = [-10 10 2.5; -10 10 2.5; 0 20 2.5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify actions to be taken
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

find_edge_sections_flag         = true ;
write_edge_sections_text_flag   = false ;
write_edge_sections_netcdf_flag = false ;
plot_edge_sections_flag         = true ;
compute_transport_flag          = true ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Begin main code.  Normally this does not need to change.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all

% change the coordinate range to be 0 to 360.
sectionCoord(:,2) = mod(sectionCoord(:,2),360);
sectionCoord(:,4) = mod(sectionCoord(:,4),360);

for iSim = 1:length(sim)

  fprintf(['**** simulation: ' sim(iSim).dir '  ' sim(iSim).netcdf_file '\n'])
  unix(['mkdir -p docs/' sim(iSim).netcdf_file '_dir/f']);
  fid_latex = fopen('temp.tex','w');
  fprintf(fid_latex,['%% file created by plot_mpas_cross_sections, ' date '\n\n']);

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
  %  Write section edge information to a netcdf file
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if write_edge_sections_text_flag 
    write_edge_sections_text...
       (sim(iSim).dir, sim(iSim).sectionEdgeIndex, ...
        sim(iSim).sectionEdgeSign, sim(iSim).nEdgesInSection, ...
        sectionText,sectionAbbreviation,sectionCoord)
  end
  
  if write_edge_sections_netcdf_flag 
    write_edge_sections_netcdf...
       (sim(iSim).dir, sim(iSim).sectionEdgeIndex, ...
        sim(iSim).sectionEdgeSign, sim(iSim).nEdgesInSection, ...
        sectionText,sectionAbbreviation,sectionCoord)
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
  %  Load large variables from netcdf file
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if compute_transport_flag
    [sim(iSim).sectionData] = load_large_variables_edge ...
       (wd,sim(iSim).dir,sim(iSim).netcdf_file, var_name,var_conv_factor, ...
        sim(iSim).sectionEdgeIndex, sim(iSim).nEdgesInSection);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Compute transport through each section
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if compute_transport_flag
    sim(iSim).tr_total = compute_transport ...
     (wd,sim(iSim).dir,sim(iSim).netcdf_file, ...
      sim(iSim).sectionEdgeIndex, sim(iSim).sectionEdgeSign, sim(iSim).nEdgesInSection, ...
      sim(iSim).sectionData,sectionText,sectionAbbreviation);
    
    if iSim==1
      tr_total = sim(iSim).tr_total';
    else
      tr_total = [tr_total; sim(iSim).tr_total'];
    end
    
  end

  fclose(fid_latex);
  
end % iSim

% tr_total
mean_transport = mean(tr_total,1);
%fprintf(['mean over time, ' sim(1).dir ' \n' ])
fprintf('%10.2f',mean_transport)
fprintf([' mean, ' sim(1).dir ' \n'])

var_transport = var(tr_total,1);
%fprintf(['variance over time, ' sim(1).dir ' \n' ])
fprintf('%10.2f',var_transport)
fprintf([' var,  ' sim(1).dir ' \n'])

std_transport = std(tr_total,1);
%fprintf(['stdev over time, ' sim(1).dir ' \n' ])
fprintf('%10.2f',std_transport)
fprintf([' std,  ' sim(1).dir ' \n'])

min_transport = min(tr_total,[],1);
%fprintf(['minimum over time, ' sim(1).dir ' \n' ])
fprintf('%10.2f',min_transport)
fprintf([' min,  ' sim(1).dir ' \n'])

max_transport = max(tr_total,[],1);
%fprintf(['maximum over time, ' sim(1).dir ' \n' ])
fprintf('%10.2f',max_transport)
fprintf([' max,  ' sim(1).dir ' \n'])

filename = ['data/' sim(1).dir '_' char(var_name) '_small_data_file.mat']
clear sim
save(filename)

end % letter
