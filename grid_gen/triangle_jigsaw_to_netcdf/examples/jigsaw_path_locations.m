function path = jigsaw_path_locations
% adds path locations to MATLAB path
% note need path to jigsaw-geo-matlab and jigsaw-geo-matlab/mesh-util in MATLAB path.
% need absolute paths here
base = '/Users/pwolfram/Documents/GridGen/MultiscaleMeshGen/';
addpath([base, 'dual-mesh/'])
addpath([base, 'jigsaw-geo-matlab/'])
addpath([base, 'jigsaw-geo-matlab/jigsaw/'])
addpath([base, 'jigsaw-geo-matlab/mesh-util/'])
path = [base, 'jigsaw-geo-matlab/'];
end
