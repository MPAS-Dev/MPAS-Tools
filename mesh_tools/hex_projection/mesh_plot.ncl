load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  ;******************************************************************************

  wks = gsn_open_wks("pdf","mesh_ncl_plot")

  f_in = addfile("mpas_hex_mesh.nc","r")

  ;******************************************************************************

  res                      = True
;  res@gsnMaximize          = True
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@gsnSpreadColors      = True

;  res@mpProjection         = "CylindricalEquidistant"
  res@mpProjection         = "Orthographic"
  res@mpDataBaseVersion    = "MediumRes"
  res@mpCenterLatF         = 38.5
;  res@mpCenterLatF         = 0.0
  res@mpCenterLonF         = -97.5
 
  res@mpMinLonF            = -157.5
  res@mpMaxLonF            = -37.5
  res@mpMinLatF            =  10.
  res@mpMaxLatF            =  70.

  res@cnFillOn             = True
  res@cnFillMode           = "AreaFill"
  res@cnLinesOn            = False
  ;res@cnLineLabelsOn       = True
  ;res@cnLineThicknessF     = 2.

  res@cnLevelSelectionMode = "ManualLevels"

  res@tiMainFontHeightF    = 0.016
  res@tiMainString         = "MPAS regional mesh"
  res@gsnLeftStringFontHeightF = 0.015
  res@gsnRightStringFontHeightF = 0.015

  ;res@lbLabelStride = 2

  ;-----1. contour for global z500 -------------------------------------------

  ; gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
  ; gsn_define_colormap(wks,"BlAqGrYeOrReVi200")
  gsn_define_colormap(wks,"perc2_9lev")

  ;res@gsnSpreadColors      = False
  ;res@gsnSpreadColorStart  = 2
  ;res@gsnSpreadColorEnd    = 201

  r2d = 180.D/(4.*atan(1.))
  lat1d = f_in->latCell(:)*r2d
  lon1d = f_in->lonCell(:)*r2d

  res@sfXArray = lon1d
  res@sfYArray = lat1d


  ibdymask = f_in->bdyMaskCell
  bdymask = int2dble(ibdymask)
  nCells = dimsizes(lat1d)

    res@cnMinLevelValF       =  0.5
    res@cnMaxLevelValF       =  7.5
    res@cnLevelSpacingF      =  1.0


;    plot = gsn_csm_contour_map_ce(wks,bdymask,res)
    plot = gsn_csm_contour_map(wks,bdymask,res)
    draw(plot)
    frame(wks)


  end 

