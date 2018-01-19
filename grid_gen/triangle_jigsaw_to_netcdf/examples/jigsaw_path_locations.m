function path = jigsaw_path_locations
% adds path locations to MATLAB path
% note need path to jigsaw-geo-matlab and jigsaw-geo-matlab/mesh-util in MATLAB path.
addpath('../../../../dual-mesh/')
addpath('../../../../jigsaw-geo-matlab/')
addpath('../../../../jigsaw-geo-matlab/jigsaw/')
addpath('../../../../jigsaw-geo-matlab/mesh-util/')
path = '../../../../jigsaw-geo-matlab/';
end
