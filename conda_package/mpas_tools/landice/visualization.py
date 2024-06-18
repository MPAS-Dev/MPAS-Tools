import numpy as np
import csv
from netCDF4 import Dataset
from scipy.interpolate import LinearNDInterpolator
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colorbar import Colorbar
from matplotlib.colors import Normalize, TwoSlopeNorm


def plot_transect(data_path, variable, ax, times=[0],
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

    ax : matplotlib.axes._axes.Axes
        Axes on which to plot variable

    times : list of ints, optional
        Time indices at which to plot variable.

    coords_file : str, optional
        Path to file containing coordinates for transect

    x : list of floats, optional
        x coordinates defining transect if not using coords_file

    y : list of floats, optional
        y coordinates defining transect if not using coords_file

    time_cmap_name : str, optional
        Name of matplotlib colormap for multiple time levels odataset.variable

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


def plot_map(data_path, variable, ax, time=0, cmap=None,
             vmin=None, vmax=None, log_plot=False, mesh_file=None,
             triangles=None, plot_grounding_line=False,
             variable_name=None):
    """
    Plot map of MALI output either by specifying a variable name or
    a pre-computed field.

    Parameters
    ----------
    data_path : str
        Path to .nc file containing variables to plot. May contain
        MALI mesh fields, or you can use the mesh_file argument.

    variable : str or numpy.ndarray
        MALI variable to plot. If a string is specified, the variable with that
        name will be read from the .nc file at data_path. If a numpy array is
        given, that array will be plotted directly.

    ax : matplotlib.axes._axes.Axes
        Axes on which to plot variable

    time : int, optional
        Time index at which to plot variable.

    cmap : str, optional
        Name of matplotlib colormap for multiple time levels of variable

    vmin : float, optional
        Minimum value to use for plotting. If not specified, the 1st
        percentile value will be used (not weighted by cell area).

    vmax : float, optional
        Maximum value to use for plotting. If not specified, the 99th
        percentile value will be used (not weighted by cell area).

    log_plot : boolean, optional
        Whether to plot log10(variable)

    mesh_file : str, optional
        Optional file used to specify mesh variables. If not provided, mesh
        variables will be read from data_path

    triangles : matplotlib.tri._triangulation.Triangulation, optional
        Triangles to use for plotting. If not specified,
        they will be calculated.

    plot_grounding_line : boolean, optional
        Whether to plot the grounding line along with variable.

    variable_name : str
        Name to use for colorbar if specifying `variable` as a numpy array.

    Returns
    -------
    var_plot : matplotlib.collections.PolyCollection
        Plot of variable

    cbar : matplotlib.colorbar.Colorbar
        Colorbar object corresponded to var_plot

    gl_plot : matplotlib.tri._tricontour.TriContourSet or None
        Contour object of the grounding line, if plot_grounding_line=True
    """

    sec_per_year = 60. * 60. * 24. * 365.

    dataset = Dataset(data_path)
    dataset.set_auto_mask(False)

    if triangles is None:
        if mesh_file is None:
           mesh = Dataset(data_path)
        else:
           mesh = Dataset(mesh_file)

        triangles, tri_mask = _get_triangles(mesh)
        mesh.close()

    if type(variable) is str:
        variable_name = variable
        if variable == 'observedSpeed':
            var_to_plot = np.sqrt(dataset.variables['observedSurfaceVelocityX'][:]**2 +
                                  dataset.variables['observedSurfaceVelocityY'][:]**2)
        else:
            var_to_plot = dataset.variables[variable][:]
    else:
        var_to_plot = variable

    if len(np.shape(var_to_plot)) == 1:
       var_to_plot = var_to_plot.reshape((1, np.shape(var_to_plot)[0]))

    if 'Speed' in variable:
        units = 'm yr^{-1}'
        var_to_plot *= sec_per_year
    else:
        try:
            units = f'({dataset.variables[variable].units})'
        except AttributeError:
            units = "{}"  # This leaves out units on the colorbar
        except TypeError:
            units = "{}"

    default_colors = {'thickness' : 'Blues',
                     'surfaceSpeed' : 'plasma',
                     'basalSpeed' : 'plasma',
                     'bedTopography' : 'BrBG',
                     'floatingBasalMassBalApplied' : 'cividis'
                    }
    # List of diverging colormaps for use in plotting bedTopography.
    # I don't see a way around hard-coding this.
    div_color_maps = ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
                          'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']

    if cmap is None:
        if variable_name in default_colors.keys():
            cmap = default_colors[variable]
        else:
            cmap = "viridis"

    if log_plot:
        var_to_plot = np.log10(var_to_plot)
        # Get rid of +/- inf values that ruin vmin and vmax
        # calculations below.
        var_to_plot[np.isinf(var_to_plot)] = np.nan
        colorbar_label_prefix = 'log10 '
    else:
        colorbar_label_prefix = ''

    # Set lower and upper bounds for plotting
    if vmin is None:
        # 0.1 m/yr is a pretty good lower bound for speed
        first_quant = np.nanquantile(var_to_plot[time, :], 0.01)
        if 'Speed' in variable and log_plot:
            vmin = max(first_quant, -1.)
        else:
            vmin = first_quant
    if vmax is None:
        vmax = np.nanquantile(var_to_plot[time, :], 0.99)
    # Plot bedTopography on an asymmetric colorbar if appropriate
    if ( (variable_name == 'bedTopography') and
         (np.nanquantile(var_to_plot[time, :], 0.99) > 0.) and
         (cmap in div_color_maps) ):
        norm = TwoSlopeNorm(vmin=vmin, vmax=vmax, vcenter=0.)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    var_plot = ax.tripcolor(
        triangles, var_to_plot[time, :], cmap=cmap,
        shading='flat', norm=norm)
    ax.set_aspect('equal')

    cbar = plt.colorbar(ax=ax, mappable=var_plot,
                        orientation='vertical',
                        label=f'{colorbar_label_prefix}{variable_name} ${units}$')
    if plot_grounding_line:
        valid_masks, grounding_line_mask, _, _, _ = _calculate_masks(dataset)
        if valid_masks:
            gl_plot = ax.tricontour(triangles, grounding_line_mask[time, :],
                                    levels=[0.9999], colors='white',
                                    linestyles='solid')
        else:
            gl_plot = None
    else:
        gl_plot = None

    return var_plot, cbar, gl_plot


def plot_grounding_lines(data_paths, ax, times=[0],
                         cmap="plasma_r", mesh_file=None,
                         triangles=None):
    """
    Plot MALI grounding line at arbitrary number of times.

    Parameters
    ----------
    data_paths : str or list of str
        Path(s) to MALI file. May contain MALI mesh fields,
        or you can use the mesh_file argument.

    ax : matplotlib.axes._axes.Axes
        Axes on which to plot variable

    time : list of ints, optional
        Time indices at which to plot variable.

    cmap : str, optional
        Name of matplotlib colormap for multiple time levels of variable

    mesh_file : str, optional
        Optional file used to specify mesh variables. If not provided, mesh
        variables will be read from data_path

    triangles : matplotlib.tri._triangulation.Triangulation, optional
        Triangles to use for plotting. If not specified,
        they will be calculated.

    Returns
    -------
    gl_plots : list of matplotlib.tri._tricontour.TriContourSet
        List of grounding line contour objects
    """

    # Ensure data_path is a list for flexibility in loops
    if type(data_paths) != list:
        data_paths = [data_paths]

    # If triangles are not specified, use first
    # data file to define mesh, or use mesh_file
    if triangles is None:
        if mesh_file is None:
            mesh = Dataset(data_paths[0])
        else:
            mesh = Dataset(mesh_file)

        triangles, tri_mask = _get_triangles(mesh)
        mesh.close()

    # Loop over all files and time levels to create
    # lists of the grounding lines, and their associated years
    plot_times = []
    grounding_line_masks = []
    gl_plots = []
    for file in data_paths:
        f = Dataset(file)
        f.set_auto_mask(False)
        if 'daysSinceStart' in f.variables.keys():
            yr = f.variables['daysSinceStart'][times] / 365.
        else:
            yr = times
        plot_times.append(yr)
        valid_masks, grounding_line_mask, _, _, _ = _calculate_masks(f)
        if valid_masks:
            for time in times:
                grounding_line_masks.append(grounding_line_mask[time, :])
        f.close()

    # Determine mapping between plot time and colormap.
    # If just plotting one time, then use maximum value
    # of the specified colormap.
    n_times = len(times) * len(data_paths)
    gl_cmap = plt.get_cmap(cmap)
    if n_times > 1:
        plot_times = np.squeeze(np.ravel(plot_times))
        plot_times_norm = ( (plot_times - np.min(plot_times)) /
                           np.max(plot_times - np.min(plot_times)) )
    else:
        plot_times_norm = np.ones_like(plot_times)

    time_colors = gl_cmap(plot_times_norm)

    for ii, mask in enumerate(grounding_line_masks):
        gl_plots.append(ax.tricontour(triangles, mask,
                                levels=[0.9999], linestyles='solid',
                                colors=time_colors[ii, None]))

    if len(plot_times) > 1:
        time_cbar = plt.colorbar(cm.ScalarMappable(cmap=cmap), ax=ax,
                                 location='bottom', label="Grounding line year")
        time_cbar.ax.tick_params(labelsize=12)
        time_cbar.set_ticks(plot_times_norm)
        time_cbar.set_ticklabels(str(i) for i in plot_times)

    return gl_plots, time_cbar


def _dist(i1, i2, xCell, yCell):

    dist = ((xCell[i1]-xCell[i2])**2 + (yCell[i1]-yCell[i2])**2)**0.5
    return dist


def _calculate_masks(dataset):

    # Set bitmask values
    initial_extent_value = 1
    dynamic_value = 2
    float_value = 4
    grounding_line_value = 256
    rhoi = 910.
    rhosw = 1028.

    if 'cellMask' in dataset.variables.keys():
        valid_masks = True
        cellMask = dataset.variables["cellMask"][:]
        float_mask = (cellMask & float_value) // float_value
        dynamic_mask = (cellMask & dynamic_value) // dynamic_value
        grounding_line_mask = (cellMask & grounding_line_value) // grounding_line_value
        initial_extent_mask = (cellMask & initial_extent_value) // initial_extent_value
    elif ( 'cellMask' not in dataset.variables.keys() and
         'thickness' in dataset.variables.keys() and
         'bedTopography' in dataset.variables.keys() ):
        print(f'cellMask is not present in output file {run};'
               ' calculating masks from ice thickness')
        valid_masks = True
        grounded_mask = (dataset.variables['thickness'][:] >
                        (-rhosw / rhoi *
                         dataset.variables['bedTopography'][:]))
        # This isn't technically correct, but works for plotting
        grounding_line_mask = grounded_mask.copy()
        initial_extent_mask = (dataset.variables['thickness'][:] > 0.)
    else:
        print('cellMask and thickness and/or bedTopography'
              f' not present in output file {run};'
               ' Skipping mask calculation.')
        valid_masks = False

    return valid_masks, grounding_line_mask, float_mask, dynamic_mask, initial_extent_mask


def _get_triangles(mesh):

    xCell = mesh.variables["xCell"][:]
    yCell = mesh.variables["yCell"][:]
    dcEdge = mesh.variables["dcEdge"][:]

    triang = tri.Triangulation(xCell, yCell)
    tri_mask = np.zeros(len(triang.triangles))

    # Maximum distance in m of edges between points.
    # Make twice dcEdge to be safe
    max_dist = np.max(dcEdge) * 2.0
    for t in range(len(triang.triangles)):
        thisTri = triang.triangles[t, :]
        if _dist(thisTri[0], thisTri[1], xCell, yCell) > max_dist:
            tri_mask[t] = True
        if _dist(thisTri[1], thisTri[2], xCell, yCell) > max_dist:
            tri_mask[t] = True
        if _dist(thisTri[0], thisTri[2], xCell, yCell) > max_dist:
            tri_mask[t] = True
    triang.set_mask(tri_mask)

    return triang, tri_mask
