load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  r2d = 57.2957795             ; radians to degrees
  pi  = 3.14159265

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;
  ; Which field to plot
  ;
  plotfield = "h"
;  plotfield = "ke"
;  plotfield = "vorticity"


  ;
  ; Whether to plot horizontal wind vectors
  ;
;  horiz_winds    = True
  horiz_winds    = False

  ;
  ; Whether to do color-filled plot (filled=True) or
  ;   to plot contours of height field (filled=False)
  ;
  filled   = True
;  filled   = False

  ;
  ; The longitude of the pole-to-pole cross section
  ;
  xsec_longitude = -1.0 * pi / 6.0

  ;
  ; The number of points along the cross section
  ;
  nsec = 200

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  wks = gsn_open_wks("pdf","xsec")
  gsn_define_colormap(wks,"gui_default")

  f = addfile("output.nc","r")

  lonCell   = f->lonCell(:) * r2d
  latCell   = f->latCell(:) * r2d
  xCell     = f->xCell(:)
  yCell     = f->yCell(:)
  zCell     = f->zCell(:)
  lonVertex = f->lonVertex(:) * r2d
  latVertex = f->latVertex(:) * r2d
  xVertex = f->xVertex(:)
  yVertex = f->yVertex(:)
  zVertex = f->zVertex(:)
  lonEdge = f->lonEdge(:) * r2d
  latEdge = f->latEdge(:) * r2d
  xEdge = f->xEdge(:)
  yEdge = f->yEdge(:)
  zEdge = f->zEdge(:)
  verticesOnCell = f->verticesOnCell(:,:)
  edgesOnCell = f->edgesOnCell(:,:)
  nCellsOnCell = f->nEdgesOnCell(:)
  cellsOnCell = f->cellsOnCell(:,:)
  alpha = f->angleEdge(:)

  dims = dimsizes(latCell)
  nCells = dims(0)

  radius = 6371220.0
  xsec_latitude = 3.141592653 / 2.0
  xsec_lat_inc = 3.141592653 / (int2flt(nsec) - 1.0)

  xsecx = new((/nsec/),float)
  xsecy = new((/nsec/),float)
  xsecz = new((/nsec/),float)
  xsec_id = new((/nsec/),integer)
  xsec_edge_id = new((/nsec/),integer)
  xsec_vtx_id = new((/nsec/),integer)

  ; Compute (x,y,z) coordinates for points on cross section
  do i=0,nsec-1
     xsecx(i) = radius * cos(xsec_longitude) * cos(xsec_latitude)
     xsecy(i) = radius * sin(xsec_longitude) * cos(xsec_latitude)
     xsecz(i) = radius * sin(xsec_latitude)
     xsec_latitude = xsec_latitude - xsec_lat_inc
  end do

  ; Find cell containing first cross section point
  dmin = 2.0 * radius
  cellmin = -1
  do i=0,nCells-1
     d = sqrt((xCell(i) - xsecx(0))^2.0 + (yCell(i) - xsecy(0))^2.0 + (zCell(i) - xsecz(0))^2.0)
     if (d .lt. dmin) then
        cellmin = i
        dmin = doubletofloat(d)
     end if
  end do
  xsec_id(0) = cellmin

  ; For all other cross section points, find the grid cell containing them
  do j=1,nsec-1
     moved = 1
     do while (moved .ne. 0)
        moved = 0
        d = sqrt((xCell(cellmin) - xsecx(j))^2.0 + (yCell(cellmin) - xsecy(j))^2.0 + (zCell(cellmin) - xsecz(j))^2.0)
        do k=0,nCellsOnCell(cellmin)-1
           dn = sqrt((xCell(cellsOnCell(cellmin,k)-1) - xsecx(j))^2.0 + (yCell(cellsOnCell(cellmin,k)-1) - xsecy(j))^2.0 + (zCell(cellsOnCell(cellmin,k)-1) - xsecz(j))^2.0)
           if (dn .lt. d) then
              d = dn
              nearest = (/cellsOnCell(cellmin,k)/)-1
              moved = 1
           end if
        end do
        if (moved .eq. 1) then
           cellmin = nearest
        end if
     end do
     xsec_id(j) = cellmin
  end do

  ; For all cross section points, find the nearest vertex and edge
  do i=0,nsec-1
     iVtx = verticesOnCell(xsec_id(i),0) - 1
     iEdge = edgesOnCell(xsec_id(i),0) - 1
     xsec_edge_id(i) = iEdge
     xsec_vtx_id(i) = iVtx
     de = sqrt((xEdge(iEdge) - xsecx(i))^2.0 + (yEdge(iEdge) - xsecy(i))^2.0 + (zEdge(iEdge) - xsecz(i))^2.0)
     dv = sqrt((xVertex(iVtx) - xsecx(i))^2.0 + (yVertex(iVtx) - xsecy(i))^2.0 + (zVertex(iVtx) - xsecz(i))^2.0)
     do j=1,nCellsOnCell(xsec_id(i))-1
        iVtx = verticesOnCell(xsec_id(i),j) - 1
        iEdge = edgesOnCell(xsec_id(i),j) - 1
        de_test = sqrt((xEdge(iEdge) - xsecx(i))^2.0 + (yEdge(iEdge) - xsecy(i))^2.0 + (zEdge(iEdge) - xsecz(i))^2.0)
        dv_test = sqrt((xVertex(iVtx) - xsecx(i))^2.0 + (yVertex(iVtx) - xsecy(i))^2.0 + (zVertex(iVtx) - xsecz(i))^2.0)
        if (de_test .lt. de) then
           de = de_test
           xsec_edge_id(i) = iEdge
        end if     
        if (dv_test .lt. dv) then
           dv = dv_test
           xsec_vtx_id(i) = iVtx
        end if     
     end do
  end do

  res                      = True
  res@gsnMaximize          = True
  res@gsnSpreadColors      = True

  res@cnFillMode           = "AreaFill"

  if (filled) then
     res@cnFillOn             = True
     res@cnLinesOn            = False
     res@cnLineLabelsOn       = False
  else
     res@cnFillOn             = False
     res@cnLinesOn            = True
     res@cnLineLabelsOn       = True
  end if

  res@cnLevelSpacingF      =  50.0
  res@cnInfoLabelOn        = True

  res@lbLabelAutoStride    = True
  res@lbBoxLinesOn         = False

  res@gsnFrame          = False

  t = stringtointeger(getenv("T"))
  if (plotfield .eq. "h") then
     fld = f->h(t,:,:)
     hs  = f->h_s(:)
     ldims = dimsizes(fld)
     nVertLevels = ldims(1)
     do i=0,nVertLevels-1
        fld(:,i) = fld(:,i) + hs(:)
     end do
  end if
  if (plotfield .eq. "ke") then
     fld = f->ke(t,:,:)
     ldims = dimsizes(fld)
     nVertLevels = ldims(1)
  end if
  if (plotfield .eq. "vorticity") then
     fld = f->vorticity(t,:,:)
     ldims = dimsizes(fld)
     nVertLevels = ldims(1)
     xsec_id(:) = xsec_vtx_id(:)
  end if
  res@cnLineDashPattern = 0

  ; Extract field from along cross section into plotting array
  arr = new((/nVertLevels,nsec/),float) 
  do i=0,nsec-1
     do j=0,nVertLevels-1
        arr(j,i) = doubletofloat(fld(xsec_id(i),j))
     end do
  end do
  
  map = gsn_csm_contour(wks,arr,res)

  if (horiz_winds) then
     u   = f->u(t,:,:)
     v   = f->v(t,:,:)
     esizes = dimsizes(u)
     u_earth = new((/nVertLevels,nsec/),float)
     v_earth = new((/nVertLevels,nsec/),float)
     x_edge = new((/nVertLevels,nsec/),float)
     y_edge = new((/nVertLevels,nsec/),float)
     do i=0,nsec-1
     do j=0,nVertLevels-1
        u_earth(j,i) = doubletofloat(u(xsec_edge_id(i),j)*cos(alpha(xsec_edge_id(i))) - v(xsec_edge_id(i),j)*sin(alpha(xsec_edge_id(i))))
        v_earth(j,i) = doubletofloat(u(xsec_edge_id(i),j)*sin(alpha(xsec_edge_id(i))) + v(xsec_edge_id(i),j)*cos(alpha(xsec_edge_id(i))))
        x_edge(j,i) = i
        y_edge(j,i) = j
     end do
     end do

     wmsetp("VCH",0.0010)
     wmsetp("VRN",0.010)
     wmsetp("VRS",100.0)
     wmsetp("VCW",0.10)

     wmvect(wks, x_edge, y_edge, u_earth, v_earth)
  end if

  frame(wks)

end

