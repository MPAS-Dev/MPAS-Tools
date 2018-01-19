function refined_disk
%refined_disk setup provides a disk of high resolution

jigsaw_path_locations;

%-----------------------------------------------------------
%   Phillip Wolfram (pwolfram@lanl.gov)
%   01/19/2018
%-----------------------------------------------------------


%------------------------------------ setup files for JIGSAW
    name = mfilename('fullpath');

    opts.geom_file = ...                % GEOM file
        [name, '.msh'];

    opts.jcfg_file = ...                % JCFG file
        [name, '.jig'];

    opts.mesh_file = ...                % MESH file
        [name, '-MESH.msh'];

    opts.hfun_file = ...                % HFUN file
        [name, '-HFUN.msh'];

%------------------------------------ compute HFUN over GEOM

    dtor = pi/180. ;

    cellWidth = 240.0 ;
    cellWidthDisk = 60.0 ;
    disklon = -60.0 * dtor;
    disklat = 35.0 * dtor;
    radius = 10.0 * dtor;
    minCellWidth = ...
        min([cellWidth cellWidthDisk]);


    nlat = +180 ;
    nlon = +360 ;

    alat = linspace( -90., +90.,nlat);
    alon = linspace(-180.,+180.,nlon);

    alat = alat * dtor ;
    alon = alon * dtor ;

    density = (minCellWidth/cellWidth)^4;
    densityDisk = (minCellWidth/cellWidthDisk)^4;

    cell_dens = zeros(nlat,nlon) ;

    for i = +1 : nlon
      for j = +1 : nlat
        dist = sqrt((disklon - alon(i))^2.0 + (disklat - alat(j))^2.0) ;
        scal_temp = exp(-(dist/radius)^2.0) ;
        cell_dens(j,i) = (densityDisk-density) * scal_temp + density;
    end
  end

    cell_spac = minCellWidth./ cell_dens.^.25 ;


    figure('color','w');
    surf(alon*180/pi,alat*180/pi,cell_spac) ;
    view(2); axis image; hold on ;
    shading interp;
    title('JIGSAW HFUN data') ;

%------------------------------------ save HFUN data to file

    cell_spac = ...
        cell_spac * 2./sqrt(3.);        %%!! it seems that we need this?

    hmat.mshID = 'ELLIPSOID-GRID';
    hmat.point.coord{1} = alon ;
    hmat.point.coord{2} = alat ;
    hmat.value = cell_spac ;

    savemsh(opts.hfun_file,hmat) ;

%------------------------------------ define JIGSAW geometry

    geom.mshID = 'ELLIPSOID-MESH';
    geom.radii = 6371.*ones(3,1) ;

    savemsh(opts.geom_file,geom) ;

%------------------------------------ build mesh via JIGSAW!

    opts.hfun_scal = 'absolute';
    opts.hfun_hmax = +inf ;
    opts.hfun_hmin = +0.0 ;

    opts.mesh_dims = +2 ;               % 2-dim. simplexes

    opts.optm_qlim = 0.9375 ;

    opts.verbosity = +1 ;

    mesh = jigsaw  (opts) ;

%------------------------------------ draw mesh/cost outputs

    ang2 = triang2( ...                 % calc. tri-angles
        mesh.point.coord(:,1:3), ...
        mesh.tria3.index(:,1:3)) ;

    t_90 = max(ang2,[],2) > 90.0 ;
    t_95 = max(ang2,[],2) > 95.0 ;

    figure;
    patch ('faces',mesh.tria3.index(:,1:3), ...
        'vertices',mesh.point.coord(:,1:3), ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch ('faces',mesh.tria3.index(t_90,1:3), ...
        'vertices',mesh.point.coord(:,1:3), ...
        'facecolor','y', ...
        'edgecolor',[.2,.2,.2]) ;
    patch ('faces',mesh.tria3.index(t_95,1:3), ...
        'vertices',mesh.point.coord(:,1:3), ...
        'facecolor','r', ...
        'edgecolor',[.2,.2,.2]) ;
    set(gca,'clipping','off') ;
    title('JIGSAW TRIA mesh') ;

    drawscr(mesh.point.coord (:,1:3), ...
            [], ...
            mesh.tria3.index (:,1:3)  ...
            ) ;

    drawnow ;
    set(figure(1),'units','normalized', ...
        'position',[.35,.55,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.05,.55,.30,.35]) ;
    set(figure(3),'units','normalized', ...
        'position',[.05,.10,.30,.35]) ;
    drawnow ;

end


