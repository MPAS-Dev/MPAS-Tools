#!/usr/bin/env python
'''
Script to calculate deltat for the ISMIP6 melt param to best match observed melt rates
for Antarctica.

Note that gamma0 is set in the script, as well as the range of deltaT to search over.

The tuned basin-by-basin deltaT values are written to a file called
'basin_and_coeff_gamma0_DeltaT_quadratic_non_local_gammaX.nc'. You will want to rename it to
avoid clobbering it if the script is rerun.  In addition to saving deltaT, the file also
includes the basin info and gamma0 so it can dropped directly into a streams file to run with.

For high res meshes, there may be efficiency gains that can be implemented, or the code may
need to move to Fortran.

The target melt rate for each region is computed by spatially averaging the
observational melt field over floating cells in each region and converting
to total melt (Gt/yr).

Matt Hoffman, 9/8/2022
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import numpy as np
from netCDF4 import Dataset, num2date
from optparse import OptionParser
import matplotlib.pyplot as plt


# -----------------------
# --- Set gamma0 here ---
gamma0 = 14500.0 # MeanAnt
#gamma0 = 159000.0 # PIGL
# -----------------------
# Select range of deltaT values to search through.
# np.arange(-1.0, 1.5, 0.05) has been wide enough for MeanAnt with default gamma0 values,
# but range will be affected by gamma0 and a wider range may be necessary if any optimal deltaTs are outside this range.
# Confirm that the output deltaTs are more than increment away from boundary.
# Increments of 0.05 seems than fine enough to get a smooth function for interpolating, but use larger increments
# when testing for faster execution.
#dTs = np.arange(-1.5, 2.0, 0.25)  # MeanAnt - coarse spacing for rapid testing
dTs = np.arange(-1.0, 1.5, 0.05)  # MeanAnt - fine spacing for accurate calculation
#dTs = np.arange(-1.5, 0.0, 0.05)  # PIGL
# -----------------------


print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = OptionParser(description=__doc__)
parser.add_option("-g", dest="fileName", help="input filename that includes ice geometry information.", metavar="FILENAME")
parser.add_option("-n", "--region-file", dest="fileRegionNames", help="region definition filename with regionNames and regionCellMasks.", metavar="FILENAME")
parser.add_option("-o", dest="fcgFileName", help="ocean forcing filename. uses first time level", metavar="FILENAME")
parser.add_option("--obs-file", dest="obsFileName", help="observational melt filename (e.g. NSIDC-0792)", metavar="FILENAME")
parser.add_option("--start-year", dest="startYear", type="int", help="starting year for averaging observations (inclusive)")
parser.add_option("--end-year", dest="endYear", type="int", help="ending year for averaging observations (inclusive)")
options, args = parser.parse_args()

required_args = [
    ("-g", options.fileName),
    ("-n/--region-file", options.fileRegionNames),
    ("-o", options.fcgFileName),
    ("--obs-file", options.obsFileName),
    ("--start-year", options.startYear),
    ("--end-year", options.endYear)
]
missing = [name for name, value in required_args if value is None]
if missing:
    parser.error("Missing required argument(s): {}".format(", ".join(missing)))
if options.startYear > options.endYear:
    parser.error("--start-year must be <= --end-year")


def _decode_region_name(region_name_row):
    """Decode a region name row from a NetCDF char array into a clean string."""
    if hasattr(region_name_row, 'tobytes'):
        decoded = region_name_row.tobytes().decode('utf-8', errors='ignore')
        return decoded.replace('\x00', '').strip()
    return str(region_name_row).strip()


def _validate_regular_axis(axis, axis_name):
    """Validate that a coordinate axis is 1D, nonzero, and regularly spaced."""
    if axis.ndim != 1 or len(axis) < 2:
        raise ValueError("Observation {} axis must be 1D with at least 2 points".format(axis_name))
    delta = axis[1] - axis[0]
    if delta == 0.0:
        raise ValueError("Observation {} axis has zero spacing".format(axis_name))
    if not np.allclose(np.diff(axis), delta, rtol=0.0, atol=1.0e-9):
        raise ValueError("Observation {} axis must be regularly spaced for nearest-neighbor indexing".format(axis_name))
    return delta


def _compute_obs_mean_melt(obs_file_name, start_year, end_year):
    """Load NSIDC melt data and return a time-mean melt map over the requested years."""
    obs_melt_var = 'melt'
    obs_time_var = 'time'
    obs_x_var = 'x'
    obs_y_var = 'y'

    with Dataset(obs_file_name, 'r') as fobs:
        if obs_time_var not in fobs.variables:
            raise ValueError("Time variable '{}' not found in obs file".format(obs_time_var))
        if obs_melt_var not in fobs.variables:
            raise ValueError("Melt variable '{}' not found in obs file".format(obs_melt_var))
        if obs_x_var not in fobs.variables or obs_y_var not in fobs.variables:
            raise ValueError("Observation x/y variables '{}' and/or '{}' not found".format(obs_x_var, obs_y_var))

        obs_time = fobs.variables[obs_time_var]
        time_values = obs_time[:]
        time_units = getattr(obs_time, 'units', None)
        time_calendar = getattr(obs_time, 'calendar', 'standard')
        if time_units is None:
            raise ValueError("Observation time variable '{}' is missing units".format(obs_time_var))

        time_dates = num2date(time_values, units=time_units, calendar=time_calendar)
        time_years = np.array([date.year for date in time_dates], dtype=int)
        selected_time_inds = np.nonzero((time_years >= start_year) * (time_years <= end_year))[0]
        if len(selected_time_inds) == 0:
            raise ValueError(
                "No observation time records found in requested year range {}-{}. "
                "Available years are {}-{}.".format(start_year, end_year, time_years.min(), time_years.max()))

        obs_x = np.array(fobs.variables[obs_x_var][:], dtype=np.float64)
        obs_y = np.array(fobs.variables[obs_y_var][:], dtype=np.float64)
        melt_var = fobs.variables[obs_melt_var]

        if melt_var.ndim != 3:
            raise ValueError("Observation melt variable '{}' must have dimensions (time, y, x)".format(obs_melt_var))

        fill_value = getattr(melt_var, '_FillValue', None)
        melt_sum = np.zeros((len(obs_y), len(obs_x)), dtype=np.float64)
        melt_count = np.zeros((len(obs_y), len(obs_x)), dtype=np.int32)

        for t_ind in selected_time_inds:
            melt_slice = np.array(melt_var[t_ind, :, :], dtype=np.float64)
            if fill_value is not None:
                melt_slice[melt_slice == fill_value] = np.nan
            valid = np.isfinite(melt_slice)
            melt_sum[valid] += melt_slice[valid]
            melt_count[valid] += 1

        obs_melt_mean = np.full((len(obs_y), len(obs_x)), np.nan, dtype=np.float64)
        valid_count = melt_count > 0
        obs_melt_mean[valid_count] = melt_sum[valid_count] / melt_count[valid_count]

        obs_density = getattr(melt_var, 'density', np.nan)
        if not np.isfinite(obs_density):
            obs_density = np.nan

        return {
            'melt_mean': obs_melt_mean,
            'obs_x': obs_x,
            'obs_y': obs_y,
            'selected_count': len(selected_time_inds),
            'available_year_min': int(time_years.min()),
            'available_year_max': int(time_years.max()),
            'obs_density': float(obs_density)
        }


def _map_obs_to_model_cells(obs_melt_mean, obs_x, obs_y, x_cell, y_cell):
    """Map the gridded obs melt field to model cell centers using nearest neighbors."""
    dx = _validate_regular_axis(obs_x, 'x')
    dy = _validate_regular_axis(obs_y, 'y')

    ix = np.rint((x_cell - obs_x[0]) / dx).astype(int)
    iy = np.rint((y_cell - obs_y[0]) / dy).astype(int)

    in_bounds = (ix >= 0) * (ix < len(obs_x)) * (iy >= 0) * (iy < len(obs_y))
    obs_melt_cell = np.full(x_cell.shape, np.nan, dtype=np.float64)
    obs_melt_cell[in_bounds] = obs_melt_mean[iy[in_bounds], ix[in_bounds]]

    return obs_melt_cell, in_bounds

# Get region names from file
fn = Dataset(options.fileRegionNames, 'r')
rNamesIn = fn.variables['regionNames'][:]
regionCellMasks = fn.variables['regionCellMasks'][:]
nRegions = len(fn.dimensions['nRegions'])
# Process region names
rNames = list()
for reg in range(nRegions):
    this_name = _decode_region_name(rNamesIn[reg, :])
    rNames.append(this_name if this_name else "Region_{}".format(reg + 1))

ff = Dataset(options.fcgFileName, 'r')
zOcean = ff.variables['ismip6shelfMelt_zOcean'][:]
TFocean = ff.variables['ismip6shelfMelt_3dThermalForcing'][0,:,:]

rhoi = 910.0
rhosw = 1028.0
cp_seawater = 3.974e3
latent_heat_ice = 335.0e3

coef = (rhosw*cp_seawater/(rhoi*latent_heat_ice))**2  # in K^(-2)

f = Dataset(options.fileName,'r')
thickness = f.variables['thickness'][0,:]
bedTopography = f.variables['bedTopography'][0,:]
floatMask = ((thickness*910/1028+bedTopography)<0)*(thickness>0)
lowerSurface = -rhoi/rhosw*thickness  #only works for floating areas
nCells = len(f.dimensions['nCells'])
areaCell = f.variables['areaCell'][:]
xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]

obs_info = _compute_obs_mean_melt(
    options.obsFileName,
    options.startYear,
    options.endYear)

obsMeltCell, inBounds = _map_obs_to_model_cells(
    obs_info['melt_mean'], obs_info['obs_x'], obs_info['obs_y'], xCell, yCell)

if np.count_nonzero(inBounds) == 0:
    raise ValueError("No model cells fall within the observational grid extent")

obs_density = obs_info['obs_density'] if np.isfinite(obs_info['obs_density']) else rhoi
obsTargetMelt = np.full((nRegions,), np.nan, dtype=np.float64)
for reg in range(nRegions):
    ind = np.nonzero(floatMask * regionCellMasks[:, reg])[0]
    if len(ind) == 0:
        raise ValueError("Region '{}' has no floating cells".format(rNames[reg]))

    valid = np.isfinite(obsMeltCell[ind])
    if np.count_nonzero(valid) == 0:
        raise ValueError("Region '{}' has no valid mapped observational melt values".format(rNames[reg]))

    obsMeanMeltRate = np.sum(obsMeltCell[ind][valid] * areaCell[ind][valid]) / np.sum(areaCell[ind][valid])
    obsTargetMelt[reg] = obsMeanMeltRate * np.sum(areaCell[ind][valid]) * obs_density / 1.0e12

print("Observation summary:")
print("  file: {}".format(options.obsFileName))
print("  selected years: {}-{}".format(options.startYear, options.endYear))
print("  available years in file: {}-{}".format(obs_info['available_year_min'], obs_info['available_year_max']))
print("  selected time records: {}".format(obs_info['selected_count']))
print("  in-bounds model cells: {} / {}".format(np.count_nonzero(inBounds), nCells))
print("  observational density used for conversion: {} kg/m^3".format(obs_density))

def calcMelt(deltaT):
    """Compute cellwise and regional melt diagnostics for a given deltaT value."""
    TFdraft = np.zeros((nCells,))
    meanTFcell = np.zeros((nCells,))

    for iCell in range(nCells):
        if floatMask[iCell] == 1:
           # Linear interpolation of the thermal forcing on the ice draft depth:
           TFdraft[iCell] = np.interp(lowerSurface[iCell], np.flip(zOcean), np.flip(TFocean[iCell,:])) # flip b/c z is ordered from sfc to deep but interp needs to be increasing
    meanTF = np.zeros((nRegions,))
    for reg in range(nRegions):
        ind = np.nonzero(floatMask * regionCellMasks[:,reg])[0]
        meanTF[reg] = (TFdraft[ind]*areaCell[ind]).sum() / areaCell[ind].sum()
        meanTFcell[ind] = meanTF[reg]

    melt = gamma0 * coef * (TFdraft + deltaT) * np.absolute(meanTFcell + deltaT) # m/yr

    totalMelt = np.zeros((nRegions,))
    for reg in range(nRegions):
        ind = np.nonzero(floatMask * regionCellMasks[:,reg])[0]
        totalMelt[reg] = (melt[ind]*areaCell[ind]).sum() * rhoi / 1.0e12 #convert to Gt/yr

    return melt, totalMelt, TFdraft, meanTF

print("Considering dTs", dTs)
allMelts = np.zeros((len(dTs), nRegions))
for i, dT in enumerate(dTs):
    print(f"Calculating with dT={dT}")
    melt, allMelts[i,:], TFdraft, meanTF = calcMelt(dT)
    if False: # Generally don't want to produce a melt map for every dT, but this can be useful for debugging / special cases
       fig, axs = plt.subplots(1,1, figsize=(9,9), num=i)
       plt.scatter(xCell, yCell, s=1, c=melt, cmap='jet')
       plt.title(f'dT={dT}')
       plt.colorbar()

# Find dT that results in observed regional melt rate for each region
bestdTCells = np.zeros((nCells,))
regionCells = np.zeros((nCells,), 'i')
bestdT = np.zeros((nRegions,))
for reg in range(nRegions):
    bestdT[reg] = np.interp(obsTargetMelt[reg], allMelts[:,reg], dTs)
    print(f"{rNames[reg]}: target={obsTargetMelt[reg]:.3f} Gt/yr, best dT={bestdT[reg]:.3f}")
    bestdTCells[regionCellMasks[:,reg]==1] = bestdT[reg]
    # Also write out a region mask.
    # Note that regionCellMasks has a separate 0/1 mask for each region, whereas the ISMIP6 region mask is a single integer field where each cell is marked with the region to which it belongs.
    regionCells[regionCellMasks[:,reg]==1] = reg+1
#np.save(f'standard_bestdTs_{gamma0}.npy', bestdT)

nrow=4
ncol=4
fig4, axs4 = plt.subplots(nrow, ncol, figsize=(13, 11), num=4)
fig4.suptitle(f'melt sensitivity, gamma={gamma0}')
for reg in range(nRegions):
    plt.sca(axs4.flatten()[reg])
    plt.xlabel('delta T')
    plt.ylabel('total basin melt (Gt)')
    #plt.xticks(np.arange(22)*xtickSpacing)
    plt.grid()
    axs4.flatten()[reg].set_title(f'{rNames[reg]}: {bestdT[reg]:.5f}')
    if reg == 0:
        axX = axs4.flatten()[reg]
    else:
        axs4.flatten()[reg].sharex(axX)
    meltHere = obsTargetMelt[reg]
    plt.plot(dTs, meltHere * np.ones(dTs.shape), 'k-')
    plt.plot(dTs, allMelts[:,reg], 'b-')
fig4.tight_layout()

fig5, axs5 = plt.subplots(nrow, ncol, figsize=(13, 11), num=5)
fig5.suptitle(f'melt sensitivity, gamma={gamma0}')
for reg in range(nRegions):
    plt.sca(axs5.flatten()[reg])
    plt.xlabel('mean TF')
    plt.ylabel('total basin melt (Gt)')
    #plt.xticks(np.arange(22)*xtickSpacing)
    plt.grid()
    axs5.flatten()[reg].set_title(f'{rNames[reg]}: best TF={bestdT[reg]+meanTF[reg]:.3f}')
    if reg == 0:
        axX = axs5.flatten()[reg]
    else:
        axs5.flatten()[reg].sharex(axX)
    meltHere = obsTargetMelt[reg]
    plt.plot(dTs+meanTF[reg], meltHere * np.ones(dTs.shape), 'k-')
    plt.plot(dTs+meanTF[reg] - bestdT[reg], allMelts[:,reg], 'b-')
    plt.plot(meanTF[reg]*np.ones((2,)), [0,meltHere*2.0], 'r--')

    TFs = dTs+meanTF[reg] #+ bestdTstd[reg]
    ind = np.nonzero(floatMask * regionCellMasks[:,reg])[0]
    plt.plot(TFs, coef * gamma0 * (TFs**2 + 2.0*bestdT[reg]*TFs + bestdT[reg]**2) * areaCell[ind].sum() * rhoi / 1.0e12, 'm:')

    #np.save(f'standard_allMelts_{gamma0}_{reg}.npy', allMelts)
fig5.tight_layout()
#np.save(f'meanTF.npy', meanTF)


# write new file
foutName=f'basin_and_coeff_DeltaT_quadratic_non_local_gamma{int(gamma0)}.nc'
fout = Dataset(foutName, 'w')
fout.createDimension('nCells', nCells)
dTOut = fout.createVariable('ismip6shelfMelt_deltaT', 'd', ('nCells',))
basinOut = fout.createVariable('ismip6shelfMelt_basin', 'i', ('nCells',))
gammaOut = fout.createVariable('ismip6shelfMelt_gamma0', 'd')
dTOut[:] = bestdTCells
basinOut[:] = regionCells
gammaOut[:] = gamma0
fout.close()

plt.show()
