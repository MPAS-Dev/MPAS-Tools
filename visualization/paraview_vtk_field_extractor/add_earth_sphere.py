"""
A macro for adding an opaque sphere just under MPAS paraview data on the
sphere.

Add to ParaView via Macros -> Add new macro..

Apply by selecting the imported pdv file in the pipeline browser, then running
Macros -> add_earth_sphere

Xylar Asay-Davis
24-Aug-2017
"""
import paraview.simple


# disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find view
renderView1 = paraview.simple.GetActiveView()

# create a new 'Sphere'
sphere1 = paraview.simple.Sphere()
sphere1.Radius = 6369000.0
sphere1.ThetaResolution = 1000
sphere1.PhiResolution = 1000

# show data from sphere1
sphere1Display = paraview.simple.Show(sphere1, renderView1)
