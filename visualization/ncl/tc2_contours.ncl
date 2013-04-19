load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  r2d = 57.2957795             ; radians to degrees

  maxedges = 8 

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
  cellArea = f->areaCell(:)

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
;  res@gsnMaximize          = True
  res@gsnSpreadColors      = True

  res@sfXArray             = x
  res@sfYArray             = y
  res@sfXCellBounds        = cx
  res@sfYCellBounds        = cy

  res@cnLinesOn            = True
  res@cnLineLabelsOn       = False
  res@cnLevelSelectionMode = "AutomaticLevels"
  res@cnLevelSpacingF      = 100.
  res@cnMonoLineDashPattern = True
  res@cnLineDashPattern    = 2
  res@cnInfoLabelOn        = False

  res@lbBoxLinesOn         = False

  res@mpProjection      = "Stereographic"
  res@mpDataBaseVersion = "MediumRes"
  res@mpCenterLatF      = 90.
  res@mpCenterLonF      = 270.
  res@mpGridAndLimbOn   = False
  res@mpOutlineOn       = False
  res@mpPerimOn         = True
  res@mpAreaMaskingOn   = True
  res@mpMaskAreaSpecifiers = "Land"
  res@gsnFrame          = False

  h0   = f->h(0,:,0)
  print("Min/max of initial height field: "+min(h0)+", "+max(h0))
  map0 = gsn_csm_contour_map(wks,h0,res)

  t = stringtointeger(getenv("T"))
  res@cnLineDashPattern    = 0
  h   = f->h(t,:,0)
  print("Min/max of height field after one rotation: "+min(h)+", "+max(h))
  map1 = gsn_csm_contour_map(wks,h,res)
  tres = False

  frame(wks)

  hdiff = h - h0
  print("Min/max error after one rotation: "+min(hdiff)+", "+max(hdiff))
  map1 = gsn_csm_contour_map(wks,hdiff,res)

  frame(wks)

  cellsizes = dimsizes(h)
  intgrl = 0.0D
  do i=0,cellsizes(0)-1
  intgrl = intgrl + h0(i)*cellArea(i)
  end do
  print("Original integrated volume: "+intgrl)
  intgrl = 0.0D
  do i=0,cellsizes(0)-1
  intgrl = intgrl + h(i)*cellArea(i)
  end do
  print("Integrated volume after one rotation: "+intgrl)
  
end

