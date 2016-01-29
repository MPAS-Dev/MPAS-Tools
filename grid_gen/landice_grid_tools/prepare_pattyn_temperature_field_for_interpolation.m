clear all
close all


Tall = ncread('~/documents/antarctica_data/temperature/Antarctic_Temperature.nc','Temperature');
x=ncread('~/documents/antarctica_data/temperature/Antarctic_Temperature.nc','x');
y=ncread('~/documents/antarctica_data/temperature/Antarctic_Temperature.nc','y');
z=ncread('~/documents/antarctica_data/temperature/Antarctic_Temperature.nc','z');

nx=length(x);
ny=length(y);
nz = size(Tall,3);
zindex_to_zindex = [1:nz];


z2index_to_zindex = [1:2:nz];
z2=z(z2index_to_zindex); nz2=length(z2); % Use these to get 11 layers instead of 21.




%%  setup output file

fname = 'ais_temp_pattyn_cism_format.nc';

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
ncwrite(fname, 'stagwbndlevel', flipud(1.0-z))
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


% for kk=1:length(z2index_to_zindex)
%     kk
%     z2(kk)
%     k=z2index_to_zindex(kk)
%     z(k)
    
    Torig = Tall(:,:,k);

%%

w=4;  % window size

Tclean = Torig;

for j=1+w:ny-w
    for i=1+w:nx-w
        Tclean(i,j) = nanmean(nanmean(Torig(i-w:i+w, j-w:j+w)));

    end
end


Tcombo = Torig;
Tcombo(isnan(Tcombo)) = Tclean(isnan(Tcombo));

Tcombo(isnan(Tcombo)) = -1.0e36;  % put in an obvious bad value for nans so interpolation tool will pick it up 

%% now write it out
ncwrite(fname, 'tempstag', Tcombo-273.15, [1, 1, nz-kk+1, 1] )  % need to reverse the vertical indexing
%ncwrite(fname, 'tempstag', Tcombo-273.15, [1, 1, nz2-kk+1, 1] )  % need to reverse the vertical indexing



end

%%


%%

figure(5); clf; 
ax1=subplot(1,3,1); hold all
imagesc(Torig); colorbar; axis equal

ax2=subplot(1,3,3); hold all
imagesc(Tclean); colorbar; axis equal

ax3=subplot(1,3,2); hold all
imagesc(Tcombo); colorbar; axis equal

linkaxes([ax3,ax2,ax1],'xy');