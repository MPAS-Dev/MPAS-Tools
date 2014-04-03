% example_sections.m

% This file simply contains example cross sections with text names.
%
% Mark Petersen, MPAS-Ocean Team, LANL, March 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify section coordinates and text
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see exampleSections.m for more example sections.

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
  0 -80 60 0;...   % N Atlantic zonal mean
  ];

coord = [...
  0 -90 60 -60;...   % N Atlantic zonal mean
  ];

coord = [...
  -8  -140     10   -140;...   % Eq Pac Meridional
  ];

coord = [...
  -8  -170     10   -95;...   % Eq Pac Meridional
  ];

coord = [...
  0 -80 60 0;...   % N Atlantic zonal mean
  21   -80    45   -60;...   % DWBC N Atl meridional section
  21   283    32   285;...   % DWBC N Atl meridional section
  ];

coord = [...
  -35 -80 70 -1;...   % N Atlantic zonal mean
  ];

coord = [...
  -35 -97 70 -1;...   % N Atlantic zonal mean
  ];
nSections = size(coord,1);

% number of points to plot for each figure
nLat = 100;
nLon = 100;

% Direction to take mean: zonal (z) or meridional (m)
meanDirection = 'z';

% plotDepth(nSections) depth to which to plot each section, in m
plotDepth = 5000*ones(1,size(coord,1));

% For plotting, only four plots are allowed per row.
% Choose sections above for each page.
% page.name          name of this group of sections 
% sectionID          section numbers for each row of this page
page(1).name = 'NA';
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
'ke_avgVelocity'};

var_name = {...
'avgVelocityMeridional',...
};

var_conv_factor = [100 100 1]; % convert m/s to cm/s for velocities

%var_lims = [-20 20 2.0; -10 10 1.0; 0 20 2.5];
var_lims = [-1 1 .1; -10 10 1.0; 0 20 2.5];
%var_lims = [-5 5 .5; -10 10 1.0; 0 20 2.5];

%var_lims = [-110 110 10.0; -10 10 1.0; 0 20 2.5];

