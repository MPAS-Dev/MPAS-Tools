load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  plotfield = "area"
  level = 5
  winds = False
  nrows = 100
  ncols = 100
  maxedges = 9

  wks = gsn_open_wks("pdf","cells")
  gsn_define_colormap(wks,"wh-bl-gr-ye-re")

  f = addfile("grid.nc","r")

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
  nEdgesOnCell = f->nEdgesOnCell(:)

  res                      = True

  t = stringtointeger(getenv("T"))

  xpoly = new((/maxedges/), "double")
  ypoly = new((/maxedges/), "double")

  xcb = new((/4/), "float")
  ycb = new((/4/), "float")

  pres = True
  pres@gsnFrame = False
  pres@xyLineColor = "Background"
  plot = gsn_xy(wks,xCell,yCell,pres)

  if (plotfield .eq. "area") then
     fld   = f->areaCell(:)
     minfld = min(fld)
     maxfld = max(fld)
  end if
  if (plotfield .eq. "w") then
     fld   = f->w(t,:,level)
     minfld = min(fld)
     maxfld = max(fld)
  end if
  if (plotfield .eq. "t") then
     fld   = f->theta(t,:,level)
     minfld = min(fld)
     maxfld = max(fld)
  end if
  if (plotfield .eq. "qr") then
     fld   = f->qr(t,:,level)
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
     do iRow=1,nrows-2
     do iCol=1,ncols-2
     iCell = iRow*ncols+iCol
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

  if (plotfield .eq. "area" .or. plotfield .eq. "ke" .or. plotfield .eq. "t" .or. plotfield .eq. "w" .or. plotfield .eq. "qr") then
     do iRow=1,nrows-2
     do iCol=1,ncols-2
        iCell = iRow*ncols+iCol
        do i=0,nEdgesOnCell(iCell)-1
           xpoly(i) = xVertex(verticesOnCell(iCell,i)-1)
           ypoly(i) = yVertex(verticesOnCell(iCell,i)-1)
        end do
        do i=nEdgesOnCell(iCell),maxedges-1
           xpoly(i) = xpoly(0)
           ypoly(i) = xpoly(0)
        end do
        gsn_polyline(wks,plot,xpoly(0:nEdgesOnCell(iCell)-1),ypoly(0:nEdgesOnCell(iCell)-1),res);
     end do
     end do
     do iCell=10000,16083
        do i=0,nEdgesOnCell(iCell)-1
           xpoly(i) = xVertex(verticesOnCell(iCell,i)-1)
           ypoly(i) = yVertex(verticesOnCell(iCell,i)-1)
        end do
        do i=nEdgesOnCell(iCell),maxedges-1
           xpoly(i) = xpoly(0)
           ypoly(i) = xpoly(0)
        end do
        gsn_polyline(wks,plot,xpoly(0:nEdgesOnCell(iCell)-1),ypoly(0:nEdgesOnCell(iCell)-1),res);
     end do
  end if

  if (winds) then
     u   = 2.*f->u(t,:,level)
     v   = 2.*f->v(t,:,level)
     alpha = f->angleEdge(:)
     esizes = dimsizes(u)
     u_earth = new(dimsizes(u),float)
     v_earth = new(dimsizes(u),float)
     xwind = new(dimsizes(u),float)
     ywind = new(dimsizes(u),float)
     do i=0,esizes(0)-1
        u_earth(i) = doubletofloat(u(i)*cos(alpha(i)) - v(i)*sin(alpha(i)))
        v_earth(i) = doubletofloat(u(i)*sin(alpha(i)) + v(i)*cos(alpha(i)))
        xwind(i) = doubletofloat(xEdge(i))
        ywind(i) = doubletofloat(yEdge(i))
     end do
   
     wmsetp("VCH",0.0010)
     wmsetp("VRN",0.010)
     wmsetp("VRS",100.0)
     wmsetp("VCW",0.10)

     wmvect(wks, xwind, ywind, u_earth, v_earth)
  end if

  ;
  ; Draw label bar
  ;
  tres = True
  tres@txAngleF = 90.0
  tres@txFontHeightF = 0.015
  do i=2,200
     xcb(0) = 0.1 + i*0.8/198
     ycb(0) = 0.1

     xcb(1) = 0.1 + (i+1)*0.8/198
     ycb(1) = 0.1

     xcb(2) = 0.1 + (i+1)*0.8/198
     ycb(2) = 0.15

     xcb(3) = 0.1 + i*0.8/198
     ycb(3) = 0.15

     res@gsFillColor = i

     gsn_polygon_ndc(wks,xcb,ycb,res);

     j = (i-2) % 20
     if ((j .eq. 0) .or. (i .eq. 200)) then
        ff = minfld + (i-2) / scalefac
        label = sprintf("%7.3g", ff)
        gsn_text_ndc(wks, label, xcb(0), 0.05, tres)
     end if

  end do

  frame(wks)

end

