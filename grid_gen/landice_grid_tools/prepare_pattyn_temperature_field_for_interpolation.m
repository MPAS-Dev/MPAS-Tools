clear all
close all


%Tall = ncread('~/documents/antarctica_data/temperature/Antarctic_Temperature.nc','Temperature');
x=ncread('~/documents/antarctica_data/temperature/Antarctic_Temperature.nc','x');
y=ncread('~/documents/antarctica_data/temperature/Antarctic_Temperature.nc','y');
z=ncread('~/documents/antarctica_data/temperature/Antarctic_Temperature.nc','z');

load /Users/mhoffman/documents/antarctica_data/temperature/from_frank_2016_03_02/Temp_shelves.mat  % loads 'temp'.  It uses the same x/y as the above file

load /Users/mhoffman/documents/antarctica_data/temperature/zeta.mat  % loads 'zeta'.  These are the z-levels for the 21 levels, where 1 is the bottom and 0 the surface
z=(zeta');

nx=length(x);
ny=length(y);
nz = length(zeta);



zindex_to_zindex = [1:nz];


%z2index_to_zindex = [1:2:nz];
%z2=z(z2index_to_zindex); nz2=length(z2); % Use these to get 11 layers instead of 21.


%% set up the vertical levels specified by Frank - doesn't work right
% "The spacing is obtained from a polynomial of degree 2 through three points: 1, 1-0.015 and 0, where 1 is the bottom and 0 the surface."
% vals = [1, 1-0.015, 0];
% p = polyfit([0 1 2]/2.0, vals, 2);
% zs = [0:20]/20.0;
% f1 = polyval(p, zs);
% 
% figure(10); clf; hold all
% plot(zs, f1)
% plot([0 1 2]/2.0, vals, 'o')
% xlabel('normalized index')
% ylabel('normalized vertical height')

%%  setup output file

fname = '/Users/mhoffman/documents/antarctica_data/temperature/ais_temp_pattyn_cism_format.nc';

unix(['rm ' fname]);

% create dimension variables in file
nccreate(fname, 'x1', 'Dimensions', {'x1',nx}, 'Datatype', 'double')
nccreate(fname, 'y1', 'Dimensions', {'y1',ny}, 'Datatype', 'double')
nccreate(fname, 'time', 'Dimensions', {'time',1}, 'Datatype', 'double')  
nccreate(fname, 'stagwbndlevel', 'Dimensions', {'stagwbndlevel',nz}, 'Datatype', 'double')
%nccreate(fname, 'stagwbndlevel', 'Dimensions', {'stagwbndlevel',nz2}, 'Datatype', 'double')

% write dimension variables to file  
ncwrite(fname, 'x1', x)
ncwrite(fname, 'y1', y)
ncwrite(fname, 'stagwbndlevel', (z))
%ncwrite(fname, 'stagwbndlevel', flipud(1.0-z2))
ncwrite(fname, 'time', 0.0)

% create data variables in files
nccreate(fname, 'tempstag', 'Dimensions', {'x1',nx,'y1',ny,'stagwbndlevel',nz,'time',1}, 'Datatype', 'double')
%nccreate(fname, 'tempstag', 'Dimensions', {'x1',nx,'y1',ny,'stagwbndlevel',nz2,'time',1}, 'Datatype', 'double')



%% loop over levels
for kk=1:length(zindex_to_zindex)
    kk
    k=zindex_to_zindex(kk)
    z(k)

tic
% for kk=1:length(z2index_to_zindex)
%     kk
%     z2(kk)
%     k=z2index_to_zindex(kk)
%     z(k)
    
    Torig = temp(:,:,k)';
    
    Torig(Torig<(273.15-70.0))=NaN;  % there are some regions of bad data around Lambert Glacier.  Simply remove those spots.  The extrapolation mechanism will fill them back in with something reasonable.

%%  boxcar filter - create an extrapolated, smoothed version of this vertical level

w=8;  % window size  - number of grid cells to average in a boxcar sense

Tclean = zeros(size(Torig));
Torig_zeronan = Torig;
Torig_zeronan(isnan(Torig)) = 0.0;

Tnonan = ~isnan(Torig);
Tcntgood = zeros(size(Torig));
for i = -w:w;
   for j = -w:w;
      Tclean = Tclean + circshift(Torig_zeronan, [i,j]);
      Tcntgood = Tcntgood + circshift(Tnonan, [i,j]);
   end    
end
Tclean = Tclean ./ Tcntgood;


%% Create a combined version that uses the original data where we have and the extrapolated data elsewhere
Tcombo = Torig;
Tcombo(isnan(Torig)) = Tclean(isnan(Torig));

Tcombo(isnan(Tcombo)) = -1.0e36;  % put in an obvious bad value for nans so interpolation tool will pick it up 

%% now write it out
ncwrite(fname, 'tempstag', Tcombo-273.15, [1, 1, kk, 1] )  % 
%ncwrite(fname, 'tempstag', Tcombo-273.15, [1, 1, nz2-kk+1, 1] )  % need to reverse the vertical indexing, version for 11 levels


toc
end

%%


%%
c= [260 273.15];

figure(5); clf; 
ax1=subplot(1,3,1); hold all
imagesc(Torig); colorbar; axis equal
title('Torig')
caxis(c)

ax2=subplot(1,3,3); hold all
imagesc(Tclean); colorbar; axis equal
caxis(c)
title('Tclean')

ax3=subplot(1,3,2); hold all
imagesc(Tcombo); colorbar; axis equal
caxis(c)
title('Tcombo')


linkaxes([ax3,ax2,ax1],'xy');