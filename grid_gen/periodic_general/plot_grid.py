import netCDF4
import matplotlib.pyplot as plt


f=netCDF4.Dataset("grid.nc",'r')
x=f.variables['xCell']
y=f.variables['yCell']
d=f.variables['meshDensity']

#plt.plot(x,y,'o')
plt.scatter(x, y, 12, d, linewidth='0')
plt.axis('equal')
plt.colorbar()
plt.title("Cell center locations colored by density")

plt.show()
