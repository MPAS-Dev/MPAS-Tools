import netCDF4
import matplotlib.pyplot as plt
import numpy as np


f=netCDF4.Dataset("grid.nc",'r')
x=f.variables['xCell'][:]
y=f.variables['yCell'][:]
xv=f.variables['xVertex'][:]
yv=f.variables['yVertex'][:]
d=f.variables['meshDensity']
cONv=f.variables['cellsOnVertex']

print "min/max xCell", x.min(), x.max()
print "min/max yCell", y.min(), y.max()
print "min/max xVertex", xv.min(), xv.max()
print "min/max yVertex", yv.min(), yv.max()

fig = plt.figure(1)
#plt.plot(x,y,'o')
plt.scatter(x, y, 12, d, linewidth='0')
#plt.plot(xv, yv, 'x')  # optional: also plot vertex locations.
plt.axis('equal')
plt.colorbar()
plt.title("Cell center locations colored by density")


#fig = plt.figure(2)
#for v in range(len(cONv)):
#    plt.plot( np.hstack( (x[cONv[v,:]-1], x[cONv[v,0]-1]) ), np.hstack( (y[cONv[v,:]-1], y[cONv[v,0]-1]) ), '-k')
#plt.plot(xv, yv, 'xr')
#plt.title("Triangulation")

plt.show()
