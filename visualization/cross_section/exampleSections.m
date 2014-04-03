% example_sections.m

% This file simply contains example cross sections with text names.
%
% Mark Petersen, MPAS-Ocean Team, LANL, March 2014

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


% These are the cross-sections for the current plots:
% sectionText        a cell array with text describing each section
sectionText = {
'N Atlantic 70W lon',...
'N Atlantic 65W lon',...
'N Atlantic 60W lon',...
'N Atlantic 50W lon',...
'N Atlantic 40W lon',...
'N Atlantic 30W lon',...
	      };

% coord(nSections,4)  endpoints of sections, with one section per row as
%                     [startlat startlon endlat endlon]
% Traverse from south to north, and from east to west.
% Then positive velocities are eastward and northward.
coord = [...
  20   -70    44    -70;...   % N Atl Meridional
  19   -65    44    -65;...   % N Atl Meridional
   8.5 -60    46    -60;...   % N Atl Meridional
   1.8 -50    62    -50;...   % N Atl Meridional
  -3   -40    65    -40;...   % N Atl Meridional
  -5   -30    68.2  -30;...   % N Atl Meridional
  ];

sectionText = {
'N Atlantic 77W lon',...
'N Atlantic 76.4W lon',...
'N Atlantic 76W lon',...
'N Atlantic 75W lon',...
'N Atlantic 26N lat',...
'N Atlantic 26N lat',...
              };

coord = [...
  21   283.    32   283.;...   % DWBC N Atl meridional section
  21   283.6    32   283.6;...   % DWBC N Atl meridional section
  21   284    32   284;...   % DWBC N Atl meridional section
  21   285    32   285;...   % DWBC N Atl meridional section
  26.5   -77.1    26.5    -75;...   % DWBC N Atl zonal section
  26.5   -80    26.5    -14;...   % DWBC N Atl zonal section
  ];

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
 -67     0   -43.5    0;...   % S Oc Tas
 -67    30   -43.5   30;...   % S Oc Tas
 -67    60   -43.5   60;...   % S Oc Tas
 -67    90   -43.5   90;...   % S Oc Tas
 -65   -65   -55    -65;...   % Drake
 -67   147   -43.5  147;...   % S Oc Tas
  ];


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
'avgVelocityZonal',...
'avgVelocityMeridional',...
'ke_fromAvgVelocity'};
var_conv_factor = [100 100 1]; % convert m/s to cm/s for velocities
var_lims = [-20 20 2.0; -10 10 1.0; 0 20 2.5];

var_name = {...
'temperature',...
};
var_conv_factor = [1]; 
var_lims = [-2 12 .5];
