function [land_mask] = land_mask_Atlantic(lat,lon)

% Given latitude and longitude coordinates, produce a land mask
% land_mask = 1 in specified region
% land_mask = 0 elsewhere
%
% Mark Petersen, MPAS-Ocean Team, LANL, January 2013
%
%%%%%%%%%% input arguments %%%%%%%%%
% lat(nPoints)  latitude in degrees, ranging from 0 to 360 
%               or -180 to 180
% lon(nPoints)  longitude in degrees, ranging from -90 to 90
%
%%%%%%%%%% output arguments %%%%%%%%%
% land_mask(nPoints)

if size(lat) ~= size(lon)
  fprintf('Size of lat and lon must be the same.\n')
  return
end

land_mask = ones(size(lat));

% convert longitude to a range of -180 to 180
lon180 = mod(lon+180,360)-180;

% mask out south and north
land_mask(find( lat>=70)) = 0.0;
land_mask(find( lat<-35)) = 0.0;

% mask out eastern boundary
land_mask(find( lat>=-35 & lat<10 & lon180>22)) = 0.0;
land_mask(find( lat>= 10 & lat<49 & lon180> 0)) = 0.0;
land_mask(find( lat>= 49 & lat<66 & lon180>13)) = 0.0;
land_mask(find( lat>= 66 & lat<70 & lon180>30)) = 0.0;

% mask out western boundary
land_mask(find( lat>=-35 & lat< 9 & lon180<-63)) = 0.0;
land_mask(find( lat>=  9 & lat<14 & lon180<-84)) = 0.0;
land_mask(find( lat>= 14 & lat<17 & lon180<-89)) = 0.0;
land_mask(find( lat>= 17 & lat<50 & lon180<-98)) = 0.0;
land_mask(find( lat>= 50 & lat<70 & lon180<-70)) = 0.0;


