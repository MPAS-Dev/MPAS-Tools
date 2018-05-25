function coast_regional
%ANTARCTIC-COUPLED a coupled ocean/land-ice grid designed to
%resolve dynamics in the Antarctic shelf region.

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
        [name, '.msh'];

    opts.jcfg_file = ...                % JCFG file
        [name, '.jig'];

    opts.mesh_file = ...                % MESH file
        [name, '-MESH.msh'];

    opts.hfun_file = ...                % HFUN file
        [name, '-HFUN.msh'];

%------------------------------------ define JIGSAW geometry

    geom.mshID = 'ELLIPSOID-MESH';
    geom.radii = 6371 * ones(3,1);

    savemsh (opts.geom_file,geom);

%------------------------------------ compute HFUN over GEOM


    %%%%%%%%%%%%%%%%%%%%%%%%%% proposed MPAS Eddy Closure mesh

    % Here is the key to variable names:
    % LL low latitude
    % ML mid latitude
    % HL high latitude

    cellWidthLL = 30.0;
    cellWidthML = 70.0;
    cellWidthHL = 35.0;
    minCellWidth = min([cellWidthLL cellWidthML cellWidthHL]);

    dtr = pi/180.0;
    dLat = 0.1;
    lat = [-90+dLat:dLat:90-dLat]*dtr; % Latitude of point to compute density for
    dLatRadians = dLat*dtr;

    latStepFunction = 40 * dtr; % Latitude to change from LL to HL function

    lat_centLL = 15.0 * dtr; % In Radians - Position of center of transition region
    lat_centHL = 73.0 * dtr; % In Radians - Position of center of transition region

    widthLL = 6.0 * dtr; % In Radians - Width of transition region
    widthHL = 9.0 * dtr; % In Radians - Width of transition region

    densityLL = (minCellWidth/cellWidthLL)^4;
    densityML = (minCellWidth/cellWidthML)^4;
    densityHL = (minCellWidth/cellWidthHL)^4;

    density = zeros(size(lat));
    for j=1:length(lat)
      if abs(lat(j))<latStepFunction
        density(j) = ((densityLL-densityML) * (1.0 + tanh( (lat_centLL - abs(lat(j)))/ widthLL)) / 2.0) + densityML;
      else
        density(j) = ((densityML-densityHL) * (1.0 + tanh( (lat_centHL - abs(lat(j)))/ widthHL)) / 2.0) + densityHL;
      end

    end

    ECcellWidth = minCellWidth./density.^0.25;

    topo = loadmsh([jigsaw_path, 'jigsaw/geo/topo.msh']);

    xpos = topo.point.coord{1};
    ypos = topo.point.coord{2};
    zlev = reshape( ...
    topo.value,length(ypos),length(xpos));

   [nlat,nlon] = size(zlev);

    radE = geom.radii(1);

   [XPOS,YPOS] = meshgrid (xpos,ypos) ;

    hfn0 = +120. ;                      % global spacing (km)
    hfn2 = +10. ;                        % min. spacing
    hfn3 = +40. ;                        % halo spacing

    % coast parameters
    coastWlon = -125; % deg
    coastElon = -110; % deg
    coastSlat = 20; % deg
    coastNlat = 50; % deg

    dhdx = +.10;                        % max. gradients

    hfun = +hfn0*ones(nlat,nlon) ;

    % interp in ECcellWidth
    hfun = interp1(lat/dtr, ECcellWidth, YPOS, 'linear', 'extrap');

    htop = max(-hfn2*zlev/1000.,hfn2) ;
    htop = min(htop,hfn3);

    htop(zlev>+0.0) = hfn0;
    region = (YPOS > coastSlat & YPOS < coastNlat) & ...
             (XPOS > coastWlon & XPOS < coastElon);
    hfun(region) = htop(region);

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
    colorbar;
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


