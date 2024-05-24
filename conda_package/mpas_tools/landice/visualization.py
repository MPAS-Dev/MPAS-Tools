import numpy as np
import csv
from netCDF4 import Dataset
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm


def plot_transect(data_path, variable, times=[0], ax=None,
                  coords_file=None, x=None, y=None,
                  time_cmap_name="plasma",
                  temperature_cmap_name='YlGnBu_r'):
    """
    Plot transects of desired variable from MALI output.

    Parameters
    ----------
    data_path : str
        Path to .nc file containing MALI mesh and variables to plot

    variable : str
        MALI variable to plot. Can also be "geometry", which will
        calculate upper and lower surfaces from thickness and bed topography.

    times : list of ints, optional
        Time indices at which to plot variable.

    ax : matplotlib.axes._axes.Axes
        Axes on which to plot variable

    coords_file : str, optional
        Path to file containing coordinates for transect

    x : list of floats, optional
        x coordinates defining transect if not using coords_file

    y : list of floats, optional
        y coordinates defining transect if not using coords_file

    time_cmap_name : str, optional
        Name of matplotlib colormap for multiple time levels of variable

    temperature_cmap_nam : str, optional
        Name of matplotlib colormap for temperature transect
    """

    dataset = Dataset(data_path)
    xCell = dataset.variables["xCell"][:]
    yCell = dataset.variables["yCell"][:]
    nCells = dataset.dimensions['nCells'].size
    areaCell = dataset.variables["areaCell"][:]
    layerThicknessFractions = dataset.variables["layerThicknessFractions"][:]
    nVertLevels = dataset.dimensions['nVertLevels'].size
    thk = dataset.variables["thickness"][:]
    bed = dataset.variables["bedTopography"][:]
    # ensure that times is a list, otherwise indexing will cause errors
    times = list(times)
    if "daysSinceStart" in dataset.variables.keys():
        use_yrs = True
        yrs = dataset.variables["daysSinceStart"][:] / 365.
        times_list = [f'{yrs[i]:.1f}' for i in times]
    else:
        use_yrs = False
        times_list =  [str(i) for i in times]

    # Replace -1 time index with last forward index of time array.
    # It is unclear why these need to be in separate if-statements, but they do.
    if -1 in times:
        times[times.index(-1)] = int(dataset.dimensions['Time'].size - 1)
    if '-1' in times_list:
        times_list[times_list.index('-1')] = str(dataset.dimensions['Time'].size - 1)

    # Use coordinates from CSV file or x,y options, but not both.
    # CSV file takes precedent if both are present.
    if coords_file is not None:
        x = []
        y = []
        with open(options.coords_file, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')

            for row in reader:
                x.append(float(row[0]))
                y.append(float(row[1]))
        if [options.x_coords, options.y_coords] != [None, None]:
            print('-c and -x/-y options were both provided. Reading from ',
                  f'{options.coords_file} and ignoring -x and -y settings.')
        x = np.asarray(x)
        y = np.asarray(y)

    # increase sampling to match highest mesh resolution
    total_distance, = np.cumsum( np.sqrt( np.diff(x)**2. + np.diff(y)**2. ) )
    n_samples = int(round(total_distance / np.min(dataset.variables["dcEdge"][:])))
    x_interp = np.interp(np.linspace(0, len(x)-1, n_samples),
                         np.linspace(0, len(x)-1, len(x)), x)
    y_interp = np.interp(np.linspace(0, len(y)-1, n_samples),
                         np.linspace(0, len(y)-1, len(y)), y)

    d_distance = np.zeros(len(x_interp))
    for ii in np.arange(1, len(x_interp)):
        d_distance[ii] = np.sqrt( (x_interp[ii] - x_interp[ii-1])**2 +
                                  (y_interp[ii] - y_interp[ii-1])**2 )

    distance = np.cumsum(d_distance) / 1000.  # in km for plotting

    # Handle colormaps
    time_cmap = plt.get_cmap(time_cmap_name)
    time_colors = time_cmap(np.linspace(0,1,len(times)))

    if variable in ["temperature", "geometry"]:
        # Assume constant bed
        bed_interpolator = LinearNDInterpolator(
                     np.vstack((xCell, yCell)).T, bed[0,:])
        bed_transect = bed_interpolator(np.vstack((x_interp, y_interp)).T)
        ax.plot(distance, bed_transect, color='black')
        for i, time in enumerate(times):
            thk_interpolator = LinearNDInterpolator(
                np.vstack((xCell, yCell)).T, thk[time,:])
            thk_transect = thk_interpolator(np.vstack((x_interp, y_interp)).T)
            lower_surf = np.maximum( -910. / 1028. * thk_transect, bed_transect)
            lower_surf_nan = lower_surf.copy()  # for plotting
            lower_surf_nan[thk_transect==0.] = np.nan
            upper_surf = lower_surf + thk_transect
            upper_surf_nan = upper_surf.copy()  # for plotting
            upper_surf_nan[thk_transect==0.] = np.nan
            ax.plot(distance, lower_surf_nan, color=time_colors[i])
            ax.plot(distance, upper_surf_nan, color=time_colors[i])

        if variable == "temperature":
            time = times[-1]
            if len(times) > 1:
                print("Cannot plot temperature at more than one time."
                      " Only plotting temperature for final time.")
            temperature = dataset.variables["temperature"][:]
            layer_thk = np.zeros((len(thk_transect), nVertLevels + 1))
            layer_midpoints = np.zeros((len(thk_transect), nVertLevels))
            layer_interfaces = np.zeros((len(thk_transect), nVertLevels + 1))
            layer_thk[:,0] = 0.
            for ii in range(len(thk_transect)):
                layer_thk[ii,1:] = np.cumsum(layerThicknessFractions *
                                             thk_transect[ii])
                layer_midpoints[ii,:] = upper_surf[ii] - (layer_thk[ii,1:] +
                                                         layer_thk[ii,0:-1]) / 2.
                layer_interfaces[ii,:] = upper_surf[ii] - layer_thk[ii,:]

            temp_transect = np.zeros((len(x_interp), nVertLevels))
            for lev in range(nVertLevels):
                print(f'Interpolating temperature for level {lev} of {nVertLevels}. ', end="")
                temp_interpolant = LinearNDInterpolator(
                                        np.vstack((xCell, yCell)).T,
                                        temperature[time,:,lev])
                temp_transect[:, lev] = temp_interpolant(
                                            np.vstack((x_interp, y_interp)).T)
                temp_transect[temp_transect == 0.] = np.nan
            temp_plot = ax.pcolormesh(np.tile(distance, (nVertLevels+1,1)).T,
                                layer_interfaces[:,:], temp_transect[1:,:],
                                cmap=temperature_cmap_name,
                                vmin=240., vmax=273.15)
            temp_cbar = plt.colorbar(temp_plot)
            temp_cbar.set_label('Temperature (K)')
            temp_cbar.ax.tick_params(labelsize=12)

    else:
        var = dataset.variables[variable][:]
        for i, time in enumerate(times):
            var_interpolator = LinearNDInterpolator(
                                  np.vstack((xCell, yCell)).T, var[i,:])
            var_transect = var_interpolator(np.vstack((x_interp, y_interp)).T)
            ax.plot(distance, var_transect, color=time_colors[i])
    if (len(times) > 1):
        time_cbar = plt.colorbar(cm.ScalarMappable(cmap=time_cmap_name), ax=ax)
        time_cbar.ax.tick_params(labelsize=12)
        if use_yrs:
            time_cbar.set_label('Year')
        else:
            time_cbar.set_label('time index')
        time_cbar.set_ticks(times / np.max(times))
        time_cbar.set_ticklabels(times_list)

