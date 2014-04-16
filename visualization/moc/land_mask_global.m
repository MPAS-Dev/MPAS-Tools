function [land_mask] = land_mask_global(lat,lon)

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

% for global, include all points
land_mask = ones(size(lat));
