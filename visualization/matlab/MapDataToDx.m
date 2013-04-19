
% This script open a output.nc file and writes out the output
% in ascii format to be read in by OpenDX

clear all

eps = 1.0e-12

ncid = netcdf.open('output.nc','nc_nowrite')

doThickness = 0;
doKE = 0;
doVorticity = 0;
doVelocity = 1;

%%%%%
%  CHECK THAT THE DIMENSION ORDER (AND NUMBER) AGREES WITH OUR OUTPUT.NC
%%%%%

[TimeName, TimeLength] = netcdf.inqDim(ncid,0);
[nCellsName, nCellsLength] = netcdf.inqDim(ncid,1);
[nEdgesName, nEdgesLength] = netcdf.inqDim(ncid,2);
[nVerticesName, nVerticesLength] = netcdf.inqDim(ncid,5);
[nVertLevelsName, nVertLevelsLength] = netcdf.inqDim(ncid,9);
[nTracersName, nTracersLength] = netcdf.inqDim(ncid,10);

TimeLength
nCellsLength
nEdgesLength
nVerticesLength
nVertLevelsLength
nTracersLength


if (doThickness == 1)
    
thicknessID = netcdf.inqVarID(ncid,'h');
work =  netcdf.getVar(ncid,thicknessID);
[thicknessName,xtype,dimids,natts] = netcdf.inqVar(ncid,thicknessID);
thickness=work;

system('rm -f ./dx/h.*.*.data')
for iLevel=1:nVertLevelsLength
for iTime=0:TimeLength-1
    iTime
    stringTime = int2str(iTime);
    stringVert = int2str(iLevel);
    FileName = strcat('./dx/', thicknessName, '.', ...
        stringVert, '.', stringTime, '.', 'data')
    for iCell=1:nCellsLength
      data = thickness(iLevel,iCell,iTime+1);
      if abs(data) < eps, data=0, end
      dlmwrite(FileName, data, ...
         'precision', '%18.10e', '-append')
    end
end
end
end

if (doKE == 1)

keID = netcdf.inqVarID(ncid,'ke');
work =  netcdf.getVar(ncid,keID);
[keName,xtype,dimids,natts] = netcdf.inqVar(ncid,keID);
ke=work;

system('rm -f ./dx/ke.*.*.data')
for iLevel=1:nVertLevelsLength
for iTime=0:TimeLength-1
    stringTime = int2str(iTime);
    stringVert = int2str(iLevel);
    FileName = strcat('./dx/', keName, '.', ...
        stringVert, '.', stringTime, '.', 'data')
    for iCell=1:nCellsLength
      data = ke(iLevel,iCell,iTime+1);
      if abs(data) < eps, data=0;, end
      dlmwrite(FileName, data, ...
         'precision', '%18.10e', '-append')
    end
end
end
end

if (doVorticity == 1)
    
vorticityID = netcdf.inqVarID(ncid,'vorticity');
work =  netcdf.getVar(ncid,vorticityID);
[vorticityName,xtype,dimids,natts] = netcdf.inqVar(ncid,vorticityID);
vorticity=work;

system('rm -f ./dx/vorticity.*.*.data')
for iLevel=1:nVertLevelsLength
for iTime=0:TimeLength-1
    stringTime = int2str(iTime);
    stringVert = int2str(iLevel);
    FileName = strcat('./dx/', vorticityName, '.', ...
        stringVert, '.', stringTime, '.', 'data')
    for iVertex=1:nVerticesLength
      data = vorticity(iLevel,iVertex,iTime+1);
      if abs(data) < eps, data=0;, end
      dlmwrite(FileName, data, ...
         'precision', '%18.10e', '-append')
    end
end
end
end

if (doVelocity == 1)

uID = netcdf.inqVarID(ncid,'uReconstructX');
work =  netcdf.getVar(ncid,uID);
[uName,xtype,dimids,natts] = netcdf.inqVar(ncid,uID);
u=work;

vID = netcdf.inqVarID(ncid,'uReconstructY');
work =  netcdf.getVar(ncid,vID);
[vName,xtype,dimids,natts] = netcdf.inqVar(ncid,vID);
v=work;

wID = netcdf.inqVarID(ncid,'uReconstructZ');
work =  netcdf.getVar(ncid,wID);
[wName,xtype,dimids,natts] = netcdf.inqVar(ncid,wID);
w=work;

system('rm -f ./dx/velocity.*.*.data')
for iLevel=1:nVertLevelsLength
for iTime=0:TimeLength-1
    stringTime = int2str(iTime);
    stringVert = int2str(iLevel);
    FileName = strcat('./dx/', 'velocity', '.', ...
        stringVert, '.', stringTime, '.', 'data')
    for iCell=1:nCellsLength
      r(1) = u(iLevel,iCell,iTime+1);
      r(2) = v(iLevel,iCell,iTime+1);
      r(3) = w(iLevel,iCell,iTime+1);
      dlmwrite(FileName, r(1), ...
         'precision', '%18.10e', '-append')
     dlmwrite(FileName, r(2), ...
         'precision', '%18.10e', '-append')
     dlmwrite(FileName, r(3), ...
         'precision', '%18.10e', '-append')
    end
end
end
end