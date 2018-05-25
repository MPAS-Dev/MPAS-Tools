load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  plotfield = "ke"

  maxedges = 6 

  wks = gsn_open_wks("pdf","cells")
  gsn_define_colormap(wks,"wh-bl-gr-ye-re")

  f = addfile("output.nc","r")

  xCell   = f->xCell(:)
  yCell   = f->yCell(:)
  zCell   = f->zCell(:)
  xEdge   = f->xEdge(:)
  yEdge   = f->yEdge(:)
  zEdge   = f->zEdge(:)
  xVertex = f->xVertex(:)
  yVertex = f->yVertex(:)
  zVertex = f->zVertex(:)
  verticesOnCell = f->verticesOnCell(:,:)
  edgesOnCell = f->edgesOnCell(:,:)
  edgesOnEdge = f->edgesOnEdge(:,:)
  verticesOnEdge = f->verticesOnEdge(:,:)
  cellsOnEdge = f->cellsOnEdge(:,:)
  cellsOnVertex = f->cellsOnVertex(:,:)
  edgesOnVertex = f->edgesOnVertex(:,:)

  res                      = True

  t = stringtointeger(getenv("T"))

  xpoly = new((/6/), "double")
  ypoly = new((/6/), "double")

;  cmap = RGBtoCmap("greymap.txt")
;  gsn_define_colormap(wks,cmap)

  pres = True
  pres@gsnFrame = False
  pres@trXMinF = 0.0
  pres@trYMinF = 0.0
  pres@trXMaxF = 201000.0
  pres@trYMaxF = 201000.0
  plot = gsn_xy(wks,xCell,yCell,pres)

  if (plotfield .eq. "tracer") then
     fld   = f->tracers(t,:,0,0)
     minfld = min(fld)
     maxfld = max(fld)
  end if
  if (plotfield .eq. "h") then
     fld   = f->h(t,:,0)
     minfld = min(fld)
     maxfld = max(fld)
  end if
  if (plotfield .eq. "ke") then
     fld   = f->ke(t,:,0)
     minfld = min(fld)
     maxfld = max(fld)
  end if
  if (plotfield .eq. "vorticity") then
     fld   = f->vorticity(t,:,0)
     minfld = min(fld)
     maxfld = max(fld)
  end if
  scalefac = 198.0/(maxfld - minfld)

  if (plotfield .eq. "vorticity") then
     do iRow=1,198
     do iCol=1,198
     iCell = iRow*200+iCol
     do iVertex=2*iCell,2*iCell+1
     do i=0,2
        xpoly(i) = xCell(cellsOnVertex(iVertex,i)-1)
        ypoly(i) = yCell(cellsOnVertex(iVertex,i)-1)
        res@gsFillColor = doubletointeger((fld(iVertex)-minfld)*scalefac)+2
     end do
     gsn_polygon(wks,plot,xpoly,ypoly,res);
     end do
     end do
     end do
  end if

  if (plotfield .eq. "h" .or. plotfield .eq. "ke") then
     do iRow=1,198
     do iCol=1,198
     iCell = iRow*200+iCol
     do i=0,5
        xpoly(i) = xVertex(verticesOnCell(iCell,i)-1)
        ypoly(i) = yVertex(verticesOnCell(iCell,i)-1)
        res@gsFillColor = doubletointeger((fld(iCell-1)-minfld)*scalefac)+2
     end do
     gsn_polygon(wks,plot,xpoly,ypoly,res);
     end do
     end do
  end if

;  u   = f->u(t,:,0)
;  v   = f->v(t,:,0)
;  alpha = f->angleEdge(:)
;  esizes = dimsizes(u)
;  u_earth = new(dimsizes(u),float)
;  v_earth = new(dimsizes(u),float)
;  xwind = new(dimsizes(u),float)
;  ywind = new(dimsizes(u),float)
;  do i=0,esizes(0)-1
;     u_earth(i) = doubletofloat(u(i)*cos(alpha(i)) - v(i)*sin(alpha(i)))
;     v_earth(i) = doubletofloat(u(i)*sin(alpha(i)) + v(i)*cos(alpha(i)))
;     xwind(i) = doubletofloat(xEdge(i))
;     ywind(i) = doubletofloat(yEdge(i))
;  end do
;
;  wmsetp("VCH",0.0010)
;  wmsetp("VRN",0.010)
;  wmsetp("VRS",100.0)
;  wmsetp("VCW",0.10)
;
;  wmvect(wks, xwind, ywind, u_earth, v_earth)

  frame(wks)

end

