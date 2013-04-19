load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  r2d = 57.2957795             ; radians to degrees

  maxedges = 7 

  wks = gsn_open_wks("pdf","cells")
  gsn_define_colormap(wks,"gui_default")

  f = addfile("output.nc","r")

  x   = f->lonCell(:) * r2d
  y   = f->latCell(:) * r2d
  lonVertex = f->lonVertex(:) * r2d
  latVertex = f->latVertex(:) * r2d
  lonEdge = f->lonEdge(:) * r2d
  latEdge = f->latEdge(:) * r2d
  verticesOnCell = f->verticesOnCell(:,:)
  alpha = f->angleEdge(:)

  sizes = dimsizes(verticesOnCell)

  cx = new((/sizes(0),maxedges/),double)
  cy = new((/sizes(0),maxedges/),double)

  do j=0,sizes(0)-1
  do i=0,maxedges-2
     cy(j,i) = latVertex(verticesOnCell(j,i)-1)
     cx(j,i) = lonVertex(verticesOnCell(j,i)-1)
  end do
  cx(j,maxedges-1) = cx(j,0)
  cy(j,maxedges-1) = cy(j,0)
  end do

  res                      = True
  res@gsnMaximize          = True
  res@gsnSpreadColors      = True

  res@sfXArray             = x
  res@sfYArray             = y
  res@sfXCellBounds        = cx
  res@sfYCellBounds        = cy

  res@cnFillOn             = True
  res@cnFillMode           = "RasterFill"
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
;  res@cnMaxLevelCount      = 22
;  res@cnLevelSelectionMode = "ManualLevels"
;  res@cnMinLevelValF       = -100.
;  res@cnMaxLevelValF       = 1000.
  res@cnLevelSpacingF      =  50.
  res@cnInfoLabelOn        = True

  res@lbLabelAutoStride    = True
  res@lbBoxLinesOn         = False

;  res@mpProjection      = "Satellite"
  res@mpDataBaseVersion = "MediumRes"
  res@mpCenterLatF      = 0.
  res@mpCenterLonF      = 0.
  res@mpGridAndLimbOn   = False
  res@mpOutlineOn       = False
  res@mpFillOn          = False
  res@mpPerimOn         = False
  res@gsnFrame          = False

  t = stringtointeger(getenv("T"))
  u   = f->u(t,:,0)
  v   = f->v(t,:,0)
  h   = f->tracers(t,:,0,1)
  map = gsn_csm_contour_map(wks,h,res)

  esizes = dimsizes(u)
  u_earth = new(dimsizes(u),float)
  v_earth = new(dimsizes(u),float)
  lat_edge = new(dimsizes(u),float)
  lon_edge = new(dimsizes(u),float)
  do i=0,esizes(0)-1
     j = i % 20
     if (j .eq. 0) then
        u_earth(i) = doubletofloat(u(i)*cos(alpha(i)) - v(i)*sin(alpha(i)))
        v_earth(i) = doubletofloat(u(i)*sin(alpha(i)) + v(i)*cos(alpha(i)))
     else
        u_earth(i) = -999.
        v_earth(i) = -999.
     end if
     lat_edge(i) = doubletofloat(latEdge(i))
     lon_edge(i) = doubletofloat(lonEdge(i))
  end do

  wmsetp("VCH",0.003)
  wmsetp("VRN",0.010)
  wmsetp("VRS",20.0)
;  wmvectmap(wks, lat_edge, lon_edge, u_earth, v_earth)

  frame(wks)
end

