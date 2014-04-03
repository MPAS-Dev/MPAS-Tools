function [mocTop] = compute_moc_from_w ...
       (vertVelocityTop, botDepth, ...
     latCell,lonCell, areaCell,transport,mocLat,landMask)

% Compute moc streamfunction from vertical velocity
%
% Mark Petersen, MPAS-Ocean Team, LANL, May 2012
%
%%%%%%%%%% input arguments %%%%%%%%%
%
%%%%%%%%%% output arguments %%%%%%%%%
% mocTop(nVertLevels,nLat)
%   data in each cross-section for each variable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Load large variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
fprintf(['** compute moc simulation: \n'])

nVertLevels = length(botDepth);
nCells = length(areaCell);
nLat = length(mocLat);

mocTop = zeros(nVertLevels+1,nLat);

for k=2:nVertLevels+1
  mocTop(k,1) = mocTop(k-1,1) + transport(k-1)*1e6;
end

for iLat = 2:nLat
  ind = find(landMask==1 & latCell>=mocLat(iLat-1) & latCell<mocLat(iLat));
  % if iLat<10;   figure(10+iLat); scatter(lonCell(ind),latCell(ind)); end
  for k=2:nVertLevels+1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %mrp temp warning: used 1 as time index on vertVelocityTop,
    %just to make sure.
    mocTop(k,iLat) = mocTop(k,iLat-1) + sum(vertVelocityTop(k,ind,1).*areaCell(ind)');
  end
end

% convert m^3/s to Sverdrup
mocTop = mocTop / 1e6;

fprintf('\n')

