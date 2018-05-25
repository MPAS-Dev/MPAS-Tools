function antarctic_coupled_vis
%ANTARCTIC-COUPLED a coupled ocean/land-ice grid designed to
%resolve dynamics in the Antarctic shelf region. 
% Visualization version

jigsaw_path = jigsaw_path_locations;

%-----------------------------------------------------------
%   Phillip Wolfram (pwolfram@lanl.gov)
%   Matthew Hoffman (mhoffman@lanl.gov)
%   Darren Engwirda (de2363@columbia.edu)
%   01/19/2018
%-----------------------------------------------------------
%------------------------------------ setup files for JIGSAW

    name = mfilename('fullpath');

    opts.geom_file = ...                % GEOM file
        [jigsaw_path, 'jigsaw/geo/EARTH-ANTARCTIC-GEOM.msh'];

    opts.jcfg_file = ...                % JCFG file
        [jigsaw_path, 'jigsaw/out/EARTH-ANTARCTIC.jig'];

    opts.mesh_file = ...                % MESH file
        [name, '-MESH.msh'];

    opts.hfun_file = ...                % HFUN file
        [name, '-HFUN.msh'];

%------------------------------------ define JIGSAW geometry

    geom.mshID = 'ELLIPSOID-MESH';
    geom.radii = 6371 * ones(3,1);

    savemsh (opts.geom_file,geom);

%------------------------------------ compute HFUN over GEOM

    topo = loadmsh([jigsaw_path, 'jigsaw/geo/topo.msh']);

    xpos = topo.point.coord{1};
    ypos = topo.point.coord{2};
    zlev = reshape( ...
    topo.value,length(ypos),length(xpos));

   [nlat,nlon] = size(zlev);

    radE = geom.radii(1);

   [XPOS,YPOS] = meshgrid (xpos,ypos) ;

    hfn0 = +240. ;                      % global spacing (km)
    hfn2 = +40.;                        % min. spacing
    hfn3 = +90.;                        % halo spacing

    dhdx = +.05;                        % max. gradients

    hfun = +hfn0*ones(nlat,nlon) ;

    htop = max(-hfn2*zlev/1000.,hfn2) ;
    htop = min(htop,hfn3);

    htop(zlev>+0.0) = hfn0;
    hfun(YPOS<-35.) = hfn3;
    hfun(YPOS<-60.) = htop(YPOS<-60.) ;

    hnew = limhfun( ...
       xpos,ypos,radE,true,hfun,dhdx) ;

    hmat.mshID = 'ELLIPSOID-GRID';
    hmat.point.coord{1} = xpos*pi/180 ;
    hmat.point.coord{2} = ypos*pi/180 ;
    hmat.value = hnew ;

    savemsh(opts.hfun_file,hmat) ;

%------------------------------------ build mesh via JIGSAW!

    opts.hfun_scal = 'absolute';
    opts.hfun_hmax = +inf ;
    opts.hfun_hmin = +0.0 ;

    opts.mesh_dims = +2 ;               % 2-dim. simplexes

    opts.optm_qlim = 0.9375 ;

    opts.verbosity = +1 ;


    mesh = jigsaw  (opts) ;

%------------------------------------ draw mesh/cost outputs

    xrad = mesh.point.coord(:,1) .^ 2 ...
         + mesh.point.coord(:,2) .^ 2 ...
         + mesh.point.coord(:,3) .^ 2 ;
    xrad = max(sqrt(xrad),eps) ;

    xlat = asin (mesh.point.coord(:,3)./xrad);
    xlon = atan2(mesh.point.coord(:,2), ...
                 mesh.point.coord(:,1)) ;

    xlat = xlat * 180 / pi;
    xlon = xlon * 180 / pi;

    xlev = interp2 (xpos,ypos,zlev,xlon,xlat);

    xmsk = xlev <= +0. | xlat < -60. ;  % point "in" mask

    tlev = xlev(mesh.tria3.index(:,1)) ...
         + xlev(mesh.tria3.index(:,2)) ...
         + xlev(mesh.tria3.index(:,3)) ;
    tlev = tlev / +3. ;                 % tria. altitudes

    tmsk = xmsk(mesh.tria3.index(:,1)) ...
         + xmsk(mesh.tria3.index(:,2)) ...
         + xmsk(mesh.tria3.index(:,3)) ;
    tmsk = tmsk > +1  ;                 % tria. "in" mask

    xrad = xrad + 15. * xlev / 1000. ;  % add fake altitude

    xlat = xlat * pi / 180. ;
    xlon = xlon * pi / 180. ;

    mesh.point.coord(:,1) = ...         % new mesh coord.'s
        xrad.*cos(xlon).*cos(xlat) ;
    mesh.point.coord(:,2) = ...
        xrad.*sin(xlon).*cos(xlat) ;
    mesh.point.coord(:,3) = ...
        xrad.*           sin(xlat) ;


    figure('color','w') ;               % HFUNCTION parts
    hnew(zlev>+0.) = inf;
    surf(XPOS,YPOS,hnew);
    view(2); axis image; hold on ;
    shading interp;
    title('JIGSAW HFUN data') ;


    figure('color','w') ;               % TRIA mesh parts
    patch ( ...
        'faces',mesh.tria3.index( tmsk,1:3), ...
        'vertices', mesh.point.coord(:,1:3), ...
        'facevertexcdata',tlev(tmsk,:), ...
        'facecolor','flat', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch ( ...
        'faces',mesh.tria3.index(~tmsk,1:3), ...
        'vertices', mesh.point.coord(:,1:3), ...
        'facecolor','w', ...
        'edgecolor','none');
    set(gca,'clipping','off') ;
    caxis([min(zlev(:))*4./3., +0.]);
    colormap('hot'); brighten(+0.75);
    title('JIGSAW TRIA mesh') ;


   [CP,CE,PV,EV] = makedual2 (  ...     % DUAL mesh parts
       mesh.point.coord(:,1:3), ...
       mesh.tria3.index(:,1:3)) ;

    cmsk = xmsk(CP(:,1)) ;
    clev = xlev(CP(:,1)) ;

    figure('color','w') ;
    drawdual2(CP( cmsk,:),CE,PV,EV,clev(cmsk));
    hold on; axis image off;
    drawdual2(CP(~cmsk,:),CE,PV,EV,'w','none');
    set(gca,'clipping','off') ;
    caxis([min(zlev(:))*4./3., +0.]);
    colormap('hot'); brighten(+0.75);
    title('JIGSAW DUAL mesh') ;

end


