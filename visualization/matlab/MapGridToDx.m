

% This script open a grid.nc file and writes out the grid description
% in ascii format to be read in by OpenDX

clear all

% begin periodic parameters
doPeriodic = 0;
dc = 1000.0;
nx = 200;
ny = 200;
% end periodic parameters

doWrite = 1
doVor = 1
doTri = 1
doVector = 1

ncid = netcdf.open('grid.nc','nc_nowrite');

if (doVor == 1)
    
    xV_id = netcdf.inqVarID(ncid,'xVertex');
    yV_id = netcdf.inqVarID(ncid,'yVertex');
    zV_id = netcdf.inqVarID(ncid,'zVertex');
    nEdgesOnCell_id = netcdf.inqVarID(ncid,'nEdgesOnCell');
    verticesOnCell_id = netcdf.inqVarID(ncid,'verticesOnCell');
    areaCell_id = netcdf.inqVarID(ncid,'areaCell');

    xV=netcdf.getVar(ncid, xV_id);
    yV=netcdf.getVar(ncid, yV_id);
    zV=netcdf.getVar(ncid, zV_id);
    nEdgesOnCell=netcdf.getVar(ncid, nEdgesOnCell_id);
    verticesOnCell=netcdf.getVar(ncid, verticesOnCell_id);
    areaCell = netcdf.getVar(ncid, areaCell_id);

    xC_id = netcdf.inqVarID(ncid,'xCell');
    yC_id = netcdf.inqVarID(ncid,'yCell');
    zC_id = netcdf.inqVarID(ncid,'zCell');

    xC=netcdf.getVar(ncid, xC_id);
    yC=netcdf.getVar(ncid, yC_id);
    zC=netcdf.getVar(ncid, zC_id);
    
    work=size(nEdgesOnCell(:,1));
    nCells=work(1)

    if (doWrite == 1)
    system('rm -f ./dx/vor.position.data');
    system('rm -f ./dx/vor.edge.data');
    system('rm -f ./dx/vor.loop.data');
    system('rm -f ./dx/vor.face.data');
    system('rm -f ./dx/vor.area.data');

    iloop=0;
    iedge=0;
    for i=1:nCells
     dlmwrite('./dx/vor.face.data', i-1, '-append');
     dlmwrite('./dx/vor.area.data', areaCell(i), ...
        'precision', '%18.10e', '-append');
       dlmwrite('./dx/vor.loop.data', iloop, ...
        'precision', '%10i', '-append');
       edge(1:nEdgesOnCell(i)) = iedge;
     
     for j=1:nEdgesOnCell(i)
         x(1,j) = xV(verticesOnCell(j,i));
         x(2,j) = yV(verticesOnCell(j,i));
         x(3,j) = zV(verticesOnCell(j,i));
     end;
     
     if (doPeriodic == 1);
         for j=1:nEdgesOnCell(i);
             dx = x(1,j)-xC(i);
             dy = x(2,j)-yC(i);
             if(abs(dx) > 0.1*nx*dc);
                 if(dx > 0);, x(1,j) = x(1,j) - nx*dc;, end;
                 if(dx < 0);, x(1,j) = x(1,j) + nx*dc;, end;
             end;
             if(abs(dy) > 0.1*ny*dc*sqrt(3)/2);
                 if(dy > 0);, x(2,j) = x(2,j) - sqrt(3)*nx*dc/2;, end;
                 if(dy < 0);, x(2,j) = x(2,j) + sqrt(3)*nx*dc/2;, end;
             end;
         end;
     end;
     
     for j=1:nEdgesOnCell(i)
         dlmwrite('./dx/vor.position.data', x(:,j), 'delimiter', '\t', ...
             'precision', '%18.10e', '-append');
         edge(j) = iedge + j - 1;
       end;
       dlmwrite('./dx/vor.edge.data', edge(1:nEdgesOnCell(i)), ...
        'delimiter', '\t', 'precision', '%10i', '-append')
       iloop = iloop + nEdgesOnCell(i);
      iedge = iedge + nEdgesOnCell(i);
    end;
    
    end;

end;

if (doTri == 1)

    xC_id = netcdf.inqVarID(ncid,'xCell');
    yC_id = netcdf.inqVarID(ncid,'yCell');
    zC_id = netcdf.inqVarID(ncid,'zCell');
    nCellsOnVertex = 3;
    cellsOnVertex_id = netcdf.inqVarID(ncid, 'cellsOnVertex');
    areaTriangle_id = netcdf.inqVarID(ncid,'areaTriangle');

    xC=netcdf.getVar(ncid, xC_id);
    yC=netcdf.getVar(ncid, yC_id);
    zC=netcdf.getVar(ncid, zC_id);
    cellsOnVertex=netcdf.getVar(ncid, cellsOnVertex_id);
    areaTriangle = netcdf.getVar(ncid, areaTriangle_id);
    
    xV_id = netcdf.inqVarID(ncid,'xVertex');
    yV_id = netcdf.inqVarID(ncid,'yVertex');
    zV_id = netcdf.inqVarID(ncid,'zVertex');

    xV=netcdf.getVar(ncid, xV_id);
    yV=netcdf.getVar(ncid, yV_id);
    zV=netcdf.getVar(ncid, zV_id);

    work=size(cellsOnVertex);
    nVertices = work(:,2)

    if (doWrite == 1)
    system('rm -f ./dx/tri.position.data');
    system('rm -f ./dx/tri.edge.data');
    system('rm -f ./dx/tri.loop.data');
    system('rm -f ./dx/tri.face.data');
    system('rm -f ./dx/tri.area.data');
    
    iloop=0;
    iedge=0;
    for i=1:nVertices
     dlmwrite('./dx/tri.face.data', i-1, '-append');
     dlmwrite('./dx/tri.area.data', areaTriangle(i), ...
        'precision', '%18.10e', '-append');
     dlmwrite('./dx/tri.loop.data', iloop, ...
        'precision', '%10i', '-append');
     edge(1:3) = iedge;
     for j=1:nCellsOnVertex
         x(1,j) = xC(cellsOnVertex(j,i));
         x(2,j) = yC(cellsOnVertex(j,i));
         x(3,j) = zC(cellsOnVertex(j,i));
     end;
     
     if (doPeriodic == 1);
         for j=1:nCellsOnVertex;
             dx = x(1,j)-xV(i);
             dy = x(2,j)-yV(i);
             if(abs(dx) > 0.1*nx*dc);
                 if(dx > 0);, x(1,j) = x(1,j) - nx*dc;, end;
                 if(dx < 0);, x(1,j) = x(1,j) + nx*dc;, end;
             end;
             if(abs(dy) > 0.1*ny*dc*sqrt(3)/2);
                 if(dy > 0);, x(2,j) = x(2,j) - sqrt(3)*nx*dc/2;, end;
                 if(dy < 0);, x(2,j) = x(2,j) + sqrt(3)*nx*dc/2;, end;
             end;
         end;
     end;
     
     for j=1:nCellsOnVertex;
         dlmwrite('./dx/tri.position.data', x(:,j), 'delimiter', '\t', ...
             'precision', '%18.10e', '-append')
         edge(j) = iedge + j - 1;
     end;
     dlmwrite('./dx/tri.edge.data', edge(1:3), ...
         'delimiter', '\t', 'precision', '%10i', '-append')
     iloop = iloop + 3;
     iedge = iedge + 3;
    end;
    
    end;

end;

if (doVector == 1)
    
    if (doWrite == 1)
        system('rm -f ./dx/vector.position.data');
        system('rm -f ./dx/vector.data');
    end;
    
    nEdgesOnCell_id = netcdf.inqVarID(ncid,'nEdgesOnCell');
    nEdgesOnCell=netcdf.getVar(ncid, nEdgesOnCell_id);
    work=size(nEdgesOnCell(:,1));
    nCells=work(1)
    
    xC_id = netcdf.inqVarID(ncid,'xCell');
    yC_id = netcdf.inqVarID(ncid,'yCell');
    zC_id = netcdf.inqVarID(ncid,'zCell');

    xC=netcdf.getVar(ncid, xC_id);
    yC=netcdf.getVar(ncid, yC_id);
    zC=netcdf.getVar(ncid, zC_id);
    
    xP = 0.0;
    yP = 0.0;
    zP = 1.0;
    
    for i=1:nCells
        
        a(1) = xC(i); 
        a(2) = yC(i);
        a(3) = zC(i);
        
        b(1) = xP;
        b(2) = yP;
        b(3) = zP;
        
        c(1) = a(2)*b(3) - a(3)*b(2);
        c(2) = a(3)*b(1) - a(1)*b(3);
        c(3) = a(1)*b(2) - a(2)*b(1);

    
        if (doWrite == 1)
        
            dlmwrite('./dx/vector.position.data', xC(i), ...
             'precision', '%18.10e', '-append')
       
            dlmwrite('./dx/vector.position.data', yC(i), ...
             'precision', '%18.10e', '-append')
       
            dlmwrite('./dx/vector.position.data', zC(i), ...
             'precision', '%18.10e', '-append')
         
            dlmwrite('./dx/vector.data', c(1), ...
             'precision', '%18.10e', '-append')
       
            dlmwrite('./dx/vector.data', c(2), ...
             'precision', '%18.10e', '-append')
       
            dlmwrite('./dx/vector.data', c(3), ...
             'precision', '%18.10e', '-append')
       
    
    end;
    
end;    
    
end;

netcdf.close(ncid)