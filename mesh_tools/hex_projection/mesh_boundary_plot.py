import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import math
try:
    from numba import njit
except ImportError:
    njit = lambda f: f

@njit
def boundary_vertices(lonVertex, latVertex, bdyMaskVertex, edgesOnVertex, verticesOnEdge, bdyMask):
    bdyLats = []
    bdyLons = []

    for startVertex in range(bdyMaskVertex.size):
        if bdyMaskVertex[startVertex] == bdyMask:
            break

    for edge in edgesOnVertex[startVertex][:]:
        if edge == -1:
            continue

        if verticesOnEdge[edge][0] != startVertex:
            neighbor = verticesOnEdge[edge][0]
        else:
            neighbor = verticesOnEdge[edge][1]

        if bdyMaskVertex[neighbor] == bdyMask:
            nextVertex = neighbor
            bdyLons.append(lonVertex[nextVertex])
            bdyLats.append(latVertex[nextVertex])
            break

    prevVertex = startVertex

    while nextVertex != startVertex:
        for edge in edgesOnVertex[nextVertex][:]:
            if edge == -1:
                continue

            if verticesOnEdge[edge][0] != nextVertex:
                neighbor = verticesOnEdge[edge][0]
            else:
                neighbor = verticesOnEdge[edge][1]

            if bdyMaskVertex[neighbor] == bdyMask and neighbor != prevVertex:
                prevVertex = nextVertex
                nextVertex = neighbor
                bdyLons.append(lonVertex[nextVertex])
                bdyLats.append(latVertex[nextVertex])
                break

    return np.asarray(bdyLons), np.asarray(bdyLats)


def mesh_extent(latField, lonField):
    earthRadius = 6371229.0
    extentLat = (max(latField) - min(latField)) * earthRadius * math.pi / 180.0
    extentLon = (max(lonField) - min(lonField)) * math.cos(cenLat * math.pi / 180.0) * earthRadius * math.pi / 180.0
    return max(extentLat, extentLon)


def map_background(cenLat, cenLon, extent):
    import cartopy.feature as cfeature
    import matplotlib.ticker as mticker

    proj = ccrs.Stereographic(cenLat, cenLon)
    ax = plt.axes(projection=proj)

    #### MGD
    # print(proj.transform_point(-110.0, 40.0, ccrs.PlateCarree()))

    ax.set_title('Limited-area mesh')

    scaling = 1.5
    ax.set_extent([-scaling * extent, scaling * extent, -scaling * extent, scaling * extent], crs=proj)

    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.1)
    ax.add_feature(cfeature.BORDERS, linewidth=0.1)
    ax.add_feature(cfeature.LAKES, linewidth=0.1)
    ax.add_feature(cfeature.RIVERS, linewidth=0.1)
    ax.add_feature(cfeature.STATES, linewidth=0.1)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
		      linewidth=0.1, color='black', alpha=1.0, linestyle='--')

    xticks = np.arange(-180, 180, 15)
    yticks = np.arange(-90, 90, 15)

    gl.ylocator = mticker.FixedLocator(yticks)
    gl.xlocator = mticker.FixedLocator(xticks)

    return proj, ax


def draw_domain(ax, bdyLats, bdyLons, bdyLatsSpec, bdyLonsSpec):
    import matplotlib.patches as mpatches

    poly_corners = np.zeros((len(bdyLons), 2), np.float64)
    poly_corners[:,0] = np.asarray(bdyLons)
    poly_corners[:,1] = np.asarray(bdyLats)

    poly = mpatches.Polygon(poly_corners, closed=True, ec='black', fill=True, lw=0.1, fc='black', alpha=0.3, transform=ccrs.Geodetic())
    ax.add_patch(poly)

    poly_corners = np.zeros((len(bdyLonsSpec), 2), np.float64)
    poly_corners[:,0] = np.asarray(bdyLonsSpec)
    poly_corners[:,1] = np.asarray(bdyLatsSpec)

    poly = mpatches.Polygon(poly_corners, closed=True, ec='black', fill=True, lw=0.1, fc='black', alpha=0.3, transform=ccrs.Geodetic())
    ax.add_patch(poly)


if __name__ == '__main__':
    from time import time
    from netCDF4 import Dataset

    t1 = time()
    f = Dataset('mpas_hex_mesh.nc')

    latVertex = np.ma.getdata(f.variables['latVertex'][:]) * 180.0 / math.pi
    lonVertex = np.ma.getdata(f.variables['lonVertex'][:]) * 180.0 / math.pi
    bdyMaskVertex = np.ma.getdata(f.variables['bdyMaskVertex'][:])
    edgesOnVertex = np.ma.getdata(f.variables['edgesOnVertex'][:]) - 1
    verticesOnEdge = np.ma.getdata(f.variables['verticesOnEdge'][:]) - 1

    f.close()
    t2 = time()
    print('Time to reading input fields: ', (t2 - t1))

    t1 = time()
    cenLat = np.sum(latVertex) / latVertex.size
    cenLon = np.sum(lonVertex) / lonVertex.size
    extent = mesh_extent(latVertex, lonVertex)
    t2 = time()
    print('Time to calculate domain extents ', (t2 - t1))

    t1 = time()
    proj, ax = map_background(cenLat, cenLon, extent)
    t2 = time()
    print('Time to set up map background: ', (t2 - t1))

    t1 = time()
    bdyLons, bdyLats = boundary_vertices(lonVertex, latVertex, bdyMaskVertex, edgesOnVertex, verticesOnEdge, 1)
    bdyLonsSpec, bdyLatsSpec = boundary_vertices(lonVertex, latVertex, bdyMaskVertex, edgesOnVertex, verticesOnEdge, 7)
    t2 = time()
    print('Time to find domain boundaries: ', (t2 - t1))

    t1 = time()
    draw_domain(ax, bdyLats, bdyLons, bdyLatsSpec, bdyLonsSpec)
    t2 = time()
    print('Time to draw boundaries: ', (t2 - t1))

    t1 = time()
    plt.savefig('mesh.pdf')
    t2 = time()
    print('Time to save figure: ', (t2 - t1))
