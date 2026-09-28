#!/usr/bin/env python
'''
Script to plot common time-series from one or more landice regionalStats files.
Currently only useful for whole-AIS simulations.

Observational Dataset Options:
  This script supports multiple observational datasets for validation. Use CLI
  options to select which datasets to plot:

  --obs-melt: Shelf melt rates (single choice)
  --obs-mass-balance: Mass balance data (SMB, outflow, net from single source)

  Run with --list-datasets to see all available datasets with citations.

  Default Behavior:
    - Shelf melt: Adusumilli 2020
    - Mass balance: Rignot 2019 (provides outflow only)

  Examples:
    # Use defaults
    python plot_regionalStats.py -1 globalStats.nc

    # Use Rignot 2013 for melt rates
    python plot_regionalStats.py --obs-melt rignot2013 -1 globalStats.nc

    # Use all Rignot 2008 data
    python plot_regionalStats.py --obs-mass-balance rignot2008 --obs-melt rignot2013 -1 globalStats.nc


    # List all available datasets with citations
    python plot_regionalStats.py --list-datasets

Matt Hoffman, 8/23/2022
Updated with CLI-configurable datasets: Sept 2026
'''

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import numpy as np
from netCDF4 import Dataset
from argparse import ArgumentParser
import matplotlib.pyplot as plt
plt.rcParams["text.hinting"] = "no_hinting"

rhoi = 910.0

# ============================================================================
# OBSERVATIONAL DATASETS REGISTRY
# ============================================================================
# This registry contains all observational datasets available for comparison.
# Each dataset includes citation information and data for all 16 ISMIP6 basins.

OBSERVATIONAL_DATASETS = {
    'mass_balance': {
        'rignot2008': {
            'citation': 'Rignot, E., Bamber, J., van den Broeke, M. et al. (2008). Nature Geosci 1, 106-110',
            'year': 2000,
            'variables': ['input', 'outflow', 'net'],
            'description': 'SMB, ice discharge, and net mass balance for year 2000',
            'data': {
                'ISMIP6BasinAAp': {'input': [60,9], 'outflow': [60,7], 'net': [0, 11]},
                'ISMIP6BasinApB': {'input': [39,5], 'outflow': [40,2], 'net': [-1,5]},
                'ISMIP6BasinBC': {'input': [73, 10], 'outflow': [77,4], 'net': [-4, 11]},
                'ISMIP6BasinCCp': {'input': [81, 13], 'outflow': [87,7], 'net':[-7,15]},
                'ISMIP6BasinCpD': {'input': [198,37], 'outflow': [207,13], 'net': [-8,39]},
                'ISMIP6BasinDDp': {'input': [93,14], 'outflow': [94,6], 'net': [-2,16]},
                'ISMIP6BasinDpE': {'input': [20,1], 'outflow': [22,3], 'net': [-2,4]},
                'ISMIP6BasinEF': {'input': [61+110,(10**2+7**2)**0.5], 'outflow': [49+80,(4**2+2**2)**0.5], 'net': [11+31,(11**2+7**2)**0.5]},
                'ISMIP6BasinFG': {'input': [108,28], 'outflow': [128,18], 'net': [-19,33]},
                'ISMIP6BasinGH': {'input': [177,25], 'outflow': [237,4], 'net': [-61,26]},
                'ISMIP6BasinHHp': {'input': [51,16], 'outflow': [86,10], 'net': [-35,19]},
                'ISMIP6BasinHpI': {'input': [71,21], 'outflow': [78,7], 'net': [-7,23]},
                'ISMIP6BasinIIpp': {'input': [15,5], 'outflow': [20,3], 'net': [-5,6]},
                'ISMIP6BasinIppJ': {'input': [8,4], 'outflow': [9,2], 'net': [-1,4]},
                'ISMIP6BasinJK': {'input': [93+142, (8**2+11**2)**0.5], 'outflow': [75+145,(4**2+7**2)**0.5], 'net': [18-4,(9**2+13**2)**0.5]},
                'ISMIP6BasinKA': {'input': [42+26,(8**2+7**2)**0.5], 'outflow': [45+28,(4**2+2**2)**0.5], 'net':[-3-1,(9**2+8**2)**0.5]}
            }
        },
        'rignot2019': {
            'citation': 'Rignot, E., Mouginot, J., Scheuchl, B., et al. (2019). PNAS 116(4), 1095-1103',
            'year': '2009-2017',
            'variables': ['input', 'outflow', 'net'],
            'description': 'Ice discharge (D09-17) for 2009-2017 period; uses Rignot 2008 SMB',
            'data': {
                'ISMIP6BasinAAp': {'input': [99.4, 1.0], 'outflow': [101.0, 1.3], 'net': [-1.6, 1.6]},
                'ISMIP6BasinApB': {'input': [116.7, 0.9], 'outflow': [112.9, 1.1], 'net': [3.8, 1.4]},
                'ISMIP6BasinBC': {'input': [75.3, 4.4], 'outflow': [76.8, 3.6], 'net': [-1.5, 5.7]},
                'ISMIP6BasinCCp': {'input': [136.2, 1.5], 'outflow': [144.6, 1.8], 'net': [-8.3, 2.3]},
                'ISMIP6BasinCpD': {'input': [247.1, 1.9], 'outflow': [268.0, 1.3], 'net': [-20.8, 2.3]},
                'ISMIP6BasinDDp': {'input': [125.0, 1.4], 'outflow': [132.2, 1.3], 'net': [-7.2, 1.9]},
                'ISMIP6BasinDpE': {'input': [47.3, 0.2], 'outflow': [49.6, 0.3], 'net': [-2.3, 0.4]},
                'ISMIP6BasinEF': {'input': [194.1, 4.2], 'outflow': [160.9, 5.7], 'net': [33.2, 7.1]},
                'ISMIP6BasinFG': {'input': [110.4, 1.8], 'outflow': [128.8, 2.0], 'net': [-18.4, 2.7]},
                'ISMIP6BasinGH': {'input': [199.3, 2.3], 'outflow': [299.9, 2.5], 'net': [-100.6, 3.4]},
                'ISMIP6BasinHHp': {'input': [63.1, 0.7], 'outflow': [78.6, 1.1], 'net': [-15.5, 1.3]},
                'ISMIP6BasinHpI': {'input': [136.9, 1.2], 'outflow': [143.3, 1.8], 'net': [-6.4, 2.2]},
                'ISMIP6BasinIIpp': {'input': [131.3, 3.9], 'outflow': [159.1, 2.0], 'net': [-27.8, 4.4]},
                'ISMIP6BasinIppJ': {'input': [24.8, 1.2], 'outflow': [24.8, 0.9], 'net': [0.0, 1.5]},
                'ISMIP6BasinJK': {'input': [259.6, 10.7], 'outflow': [257.2, 9.8], 'net': [2.4, 14.5]},
                'ISMIP6BasinKA': {'input': [53.9, 0.7], 'outflow': [50.8, 1.1], 'net': [3.0, 1.3]},
            }
        }
    },
    'shelf_melt': {
        'rignot2013': {
            'citation': 'Rignot, E., Jacobs, S., Mouginot, J., and Scheuchl, B. (2013). Science 341(6143), 266-270',
            'year': 2013,
            'description': 'Ice-shelf melting around Antarctica',
            'data': {
                'ISMIP6BasinAAp': [57.5, None],
                'ISMIP6BasinApB': [24.6, None],
                'ISMIP6BasinBC': [35.5, None],
                'ISMIP6BasinCCp': [107.9, None],
                'ISMIP6BasinCpD': [102.3, None],
                'ISMIP6BasinDDp': [22.8, None],
                'ISMIP6BasinDpE': [22.9, None],
                'ISMIP6BasinEF': [70.3, None],
                'ISMIP6BasinFG': [152.9, None],
                'ISMIP6BasinGH': [290.9, None],
                'ISMIP6BasinHHp': [76.3, None],
                'ISMIP6BasinHpI': [152.3, None],
                'ISMIP6BasinIIpp': [32.9, None],
                'ISMIP6BasinIppJ': [4.3, None],
                'ISMIP6BasinJK': [155.4, None],
                'ISMIP6BasinKA': [10.4, None]
            }
        },
        'adusumilli2020': {
            'citation': 'Adusumilli, S., Fricker, H. A., Medley, B., et al. (2020). Nature Geoscience 13(9), 616-620',
            'year': '2010-2018',
            'description': 'Interannual variations in meltwater flux from ice shelves',
            'data': {
                'ISMIP6BasinAAp': [108.6, None],
                'ISMIP6BasinApB': [7.4, 7.1],
                'ISMIP6BasinBC': [48.9, 39.9],
                'ISMIP6BasinCCp': [56.4, 52.3],
                'ISMIP6BasinCpD': [77.1, 12.1],
                'ISMIP6BasinDDp': [18.4, 8.9],
                'ISMIP6BasinDpE': [14.1, 5.7],
                'ISMIP6BasinEF': [104.3, 84.2],
                'ISMIP6BasinFG': [141.5, 41.4],
                'ISMIP6BasinGH': [327.9, 43.8],
                'ISMIP6BasinHHp': [191.4, 69.7],
                'ISMIP6BasinHpI': [143.3, 57.5],
                'ISMIP6BasinIIpp': [68.5, 99.6],
                'ISMIP6BasinIppJ': [35.3, 31.8],
                'ISMIP6BasinJK': [54.8, 123.5],
                'ISMIP6BasinKA': [42.3, 40.5]
            }
        }
    },
    'basin_names': {
        'ISMIP6BasinAAp': 'Dronning Maud Land',
        'ISMIP6BasinApB': 'Enderby Land',
        'ISMIP6BasinBC': 'Amery-Lambert',
        'ISMIP6BasinCCp': 'Phillipi, Denman',
        'ISMIP6BasinCpD': 'Totten',
        'ISMIP6BasinDDp': 'Mertz',
        'ISMIP6BasinDpE': 'Victoria Land',
        'ISMIP6BasinEF': 'Ross',
        'ISMIP6BasinFG': 'Getz',
        'ISMIP6BasinGH': 'Thwaites/PIG',
        'ISMIP6BasinHHp': 'Bellingshausen',
        'ISMIP6BasinHpI': 'George VI',
        'ISMIP6BasinIIpp': 'Larsen A-C',
        'ISMIP6BasinIppJ': 'Larsen E',
        'ISMIP6BasinJK': 'FRIS',
        'ISMIP6BasinKA': 'Brunt-Stancomb'
    }
}


# ============================================================================
# HELPER FUNCTIONS FOR DATASET MANAGEMENT
# ============================================================================

def list_available_datasets():
    """
    Print formatted list of all available observational datasets.
    Called when --list-datasets flag is used.
    """
    print("\n" + "="*70)
    print("AVAILABLE OBSERVATIONAL DATASETS")
    print("="*70)

    print("\nMass Balance Datasets (--obs-mass-balance):")
    print("-" * 70)
    for key, dataset in sorted(OBSERVATIONAL_DATASETS['mass_balance'].items()):
        print(f"\n  {key}:")
        print(f"    Variables: {', '.join(dataset['variables'])}")
        print(f"    Year/Period: {dataset['year']}")
        print(f"    Description: {dataset['description']}")
        print(f"    Citation: {dataset['citation']}")

    print("\n" + "-" * 70)
    print("Shelf Melt Datasets (--obs-melt):")
    print("-" * 70)
    for key, dataset in sorted(OBSERVATIONAL_DATASETS['shelf_melt'].items()):
        print(f"\n  {key}:")
        print(f"    Year/Period: {dataset['year']}")
        print(f"    Description: {dataset['description']}")
        print(f"    Citation: {dataset['citation']}")

    print("\n" + "="*70)
    print("Usage Notes:")
    print("  --obs-mass-balance provides SMB, outflow, and net from single source")
    print("  Use 'none' for any option to skip plotting that variable")
    print("="*70 + "\n")


def get_dataset_label(dataset_type, dataset_key):
    """
    Generate short label for dataset (for plot legends).

    Parameters
    ----------
    dataset_type : str
        'mass_balance' or 'shelf_melt'
    dataset_key : str
        Dataset key (e.g., 'rignot2013')

    Returns
    -------
    str
        Short label (e.g., "Rignot 2013")
    """
    if dataset_key is None or dataset_key == 'none':
        return ''

    if dataset_type in OBSERVATIONAL_DATASETS:
        if dataset_key in OBSERVATIONAL_DATASETS[dataset_type]:
            citation = OBSERVATIONAL_DATASETS[dataset_type][dataset_key]['citation']
            # Extract first author and year from citation
            # Format: "Author, X., et al. (YEAR). ..."
            try:
                author = citation.split(',')[0]
                year_part = citation.split('(')[1].split(')')[0]
                # Handle year ranges like "2009-2017" - use end year
                if '-' in year_part:
                    year = year_part.split('-')[1]
                else:
                    year = year_part
                return f"{author} {year}"
            except:
                return dataset_key

    return dataset_key


def validate_dataset_selections(mass_balance, melt_dataset):
    """
    Validate that selected datasets exist and are compatible.

    Exits with error message if validation fails.
    """
    errors = []

    # Check mass balance dataset exists
    if mass_balance not in OBSERVATIONAL_DATASETS['mass_balance']:
        errors.append(f"Unknown mass balance dataset: '{mass_balance}'")
        errors.append(f"  Available: {', '.join(OBSERVATIONAL_DATASETS['mass_balance'].keys())}")

    # Check melt dataset exists
    if melt_dataset not in OBSERVATIONAL_DATASETS['shelf_melt']:
        errors.append(f"Unknown shelf melt dataset: '{melt_dataset}'")
        errors.append(f"  Available: {', '.join(OBSERVATIONAL_DATASETS['shelf_melt'].keys())}")

    # Exit if any errors
    if errors:
        print("\n" + "="*60)
        print("ERROR: Dataset validation failed")
        print("="*60)
        for error in errors:
            print(f"  {error}")
        print("\nUse --list-datasets to see all available options.")
        print("="*60 + "\n")
        sys.exit(1)


def build_basin_info(mass_balance_dataset, melt_dataset):
    """
    Build ISMIP6basinInfo dictionary from selected observational datasets.

    This function creates a dictionary compatible with the existing plotting code
    by merging data from multiple dataset sources based on CLI selections.

    Parameters
    ----------
    mass_balance_dataset : str
        Dataset key for mass balance (provides SMB, outflow, net) - required
    melt_dataset : str
        Dataset key for shelf melt rates - required

    Returns
    -------
    dict
        Dictionary with structure:
        {basin_key: {'name': str, 'input': [mean, unc], 'outflow': [mean, unc],
                     'net': [mean, unc], 'shelfMelt': [[mean, unc]]}}

    Notes
    -----
    - shelfMelt is a list containing a single [mean, uncertainty] pair
    """
    basin_info = {}

    # Get all basin keys
    all_basins = OBSERVATIONAL_DATASETS['basin_names'].keys()

    for basin_key in all_basins:
        basin_info[basin_key] = {
            'name': OBSERVATIONAL_DATASETS['basin_names'][basin_key]
        }

        # Add mass balance data (SMB, outflow, net) from single dataset
        dataset = OBSERVATIONAL_DATASETS['mass_balance'][mass_balance_dataset]
        if basin_key in dataset['data']:
            # Copy all available variables from this dataset
            if 'input' in dataset['data'][basin_key]:
                basin_info[basin_key]['input'] = dataset['data'][basin_key]['input']
            if 'outflow' in dataset['data'][basin_key]:
                basin_info[basin_key]['outflow'] = dataset['data'][basin_key]['outflow']
            if 'net' in dataset['data'][basin_key]:
                basin_info[basin_key]['net'] = dataset['data'][basin_key]['net']

        # Add shelf melt dataset (single dataset in a list for compatibility)
        if basin_key in OBSERVATIONAL_DATASETS['shelf_melt'][melt_dataset]['data']:
            basin_info[basin_key]['shelfMelt'] = [
                OBSERVATIONAL_DATASETS['shelf_melt'][melt_dataset]['data'][basin_key]
            ]

    return basin_info


# ============================================================================
# CLI ARGUMENT PARSING
# ============================================================================

print("** Gathering information.  (Invoke with --help for more details. All arguments are optional)")
parser = ArgumentParser(description=__doc__)

# Existing arguments (preserved)
parser.add_argument("-1", dest="file1inName", help="input filename",
                    default="globalStats.nc", metavar="FILENAME")
parser.add_argument("-2", dest="file2inName", help="input filename",
                    metavar="FILENAME")
parser.add_argument("-3", dest="file3inName", help="input filename",
                    metavar="FILENAME")
parser.add_argument("-4", dest="file4inName", help="input filename",
                    metavar="FILENAME")
parser.add_argument("-u", dest="units",
                    help="units for mass/volume: m3, kg, Gt", default="Gt",
                    metavar="UNITS")
parser.add_argument("-n", dest="fileRegionNames",
                    help="region name filename. If not specified, will attempt to read region names from file 1.",
                    metavar="FILENAME")

# New dataset options
parser.add_argument("--obs-melt", dest="obsMeltDataset",
                    help="Shelf melt dataset (choose one). Options: rignot2013, adusumilli2020",
                    default="adusumilli2020")
parser.add_argument("--obs-mass-balance", dest="obsMassBalanceDataset",
                    help="Mass balance dataset (provides SMB, outflow, and net together). Options: rignot2008, rignot2019",
                    default="rignot2019")
parser.add_argument("--list-datasets", action="store_true",
                    help="List available datasets and exit")

options = parser.parse_args()

# Handle --list-datasets
if options.list_datasets:
    list_available_datasets()
    sys.exit(0)

# Parse dataset selections
selected_melt_dataset = options.obsMeltDataset
selected_mass_balance = options.obsMassBalanceDataset

# Validate selections
validate_dataset_selections(selected_mass_balance, selected_melt_dataset)

# Build basin info from selections
ISMIP6basinInfo = build_basin_info(selected_mass_balance, selected_melt_dataset)

# ============================================================================
# MAIN SCRIPT (REMAINDER UNCHANGED FROM ORIGINAL)
# ============================================================================

print("Using ice density of {} kg/m3 if required for unit conversions".format(rhoi))

# Build string for titles about the runs in use
runinfo=f'solid={options.file1inName}'
if options.file2inName:
    runinfo = f'{runinfo}\ndotted={options.file2inName}'
if options.file3inName:
    runinfo = f'{runinfo}\ndashed={options.file3inName}'
if options.file4inName:
    runinfo = f'{runinfo}\ndashdot={options.file4inName}'

if options.units == "m3":
   massUnit = "m$^3$"
elif options.units == "kg":
   massUnit = "kg"
elif options.units == "Gt":
   massUnit = "Gt"
else:
   sys.exit("Unknown mass/volume units")
print("Using volume/mass units of: ", massUnit)

# Get nRegions and yr from first file
f = Dataset(options.file1inName, 'r')
nRegions = len(f.dimensions['nRegions'])
yr = f.variables['daysSinceStart'][:]/365.0

# Get region names from file
if options.fileRegionNames:
   fn = Dataset(options.fileRegionNames, 'r')
   rNamesIn = fn.variables['regionNames'][:]
else:
   rNamesIn = f.variables['regionNames'][:]
# Process region names
rNamesOrig = list()
for r in range(nRegions):
    thisString = rNamesIn[r, :].tobytes().decode('utf-8').strip()  # convert from char array to string
    rNamesOrig.append(''.join(filter(str.isalnum, thisString)))  # this bit removes non-alphanumeric chars

# Parse region names to more usable names, if available
rNames = [None]*nRegions
for r in range(nRegions):
    if rNamesOrig[r] in ISMIP6basinInfo:
        rNames[r] = ISMIP6basinInfo[rNamesOrig[r]]['name']
    else:
        rNames[r] = rNamesOrig[r]

if nRegions <= 4:
    ncol = 2
elif nRegions <= 9:
    ncol = 3
elif nRegions <= 16:
    ncol = 4
elif nRegions <= 25:
    ncol = 5
else:
    sys.exit("ERROR: More than 25 regions found.  Attempting to plot this many regions is likely a bad idea.")
nrow = np.ceil(nRegions / ncol).astype('int') # Set nrow to have enough rows to plot number of regions based on ncol calculated above

# Set up Figure 1: volume stats overview
fig1, axs1 = plt.subplots(nrow, ncol, figsize=(13, 11), num=1)
fig1.suptitle(f'Mass change summary\n{runinfo}', fontsize=9)
for reg in range(nRegions):
   plt.sca(axs1.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('volume change ({})'.format(massUnit))
   plt.grid()
   axs1.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs1.flatten()[reg]
   else:
      axs1.flatten()[reg].sharex(axX)
   # plot obs if applicable
   if rNamesOrig[reg] in ISMIP6basinInfo and 'net' in ISMIP6basinInfo[rNamesOrig[reg]]:
       [mn, sig] = ISMIP6basinInfo[rNamesOrig[reg]]['net']
       label = f'grd obs ({get_dataset_label("mass_balance", selected_mass_balance)})'
       axs1.flatten()[reg].fill_between(yr, yr*(mn-sig), yr*(mn+sig),
                                       color='b', alpha=0.2, label=label)

# Set up Figure 2: grounded MB
fig2, axs2 = plt.subplots(nrow, ncol, figsize=(13, 11), num=2)
fig2.suptitle(f'Grounded mass change\n{runinfo}', fontsize=9)
for reg in range(nRegions):
   plt.sca(axs2.flatten()[reg])
   if reg // nrow == nrow-1:
      plt.xlabel('Year')
   if reg % ncol == 0:
      plt.ylabel('volume change ({})'.format(massUnit))
   plt.grid()
   axs2.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs2.flatten()[reg]
   else:
      axs2.flatten()[reg].sharex(axX)
   # plot obs if applicable
   if rNamesOrig[reg] in ISMIP6basinInfo:
       [mn, sig] = ISMIP6basinInfo[rNamesOrig[reg]]['input']
       label = f'SMB obs ({get_dataset_label("mass_balance", selected_mass_balance)})'
       axs2.flatten()[reg].fill_between(yr, yr*(mn-sig), yr*(mn+sig),
                                       color='b', alpha=0.2, label=label)

       [mn, sig] = ISMIP6basinInfo[rNamesOrig[reg]]['outflow']
       label = f'outflow obs ({get_dataset_label("mass_balance", selected_mass_balance)})'
       axs2.flatten()[reg].fill_between(yr, -yr*(mn-sig), -yr*(mn+sig),
                                       color='g', alpha=0.2, label=label)

       [mn, sig] = ISMIP6basinInfo[rNamesOrig[reg]]['net']
       label = f'net obs ({get_dataset_label("mass_balance", selected_mass_balance)})'
       axs2.flatten()[reg].fill_between(yr, yr*(mn-sig), yr*(mn+sig),
                                       color='k', alpha=0.2, label=label)


# Set up Figure 3: floating MB
fig3, axs3 = plt.subplots(nrow, ncol, figsize=(13, 11), num=3)
fig3.suptitle(f'Floating mass change\n{runinfo}', fontsize=9)
for reg in range(nRegions):
   plt.sca(axs3.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('volume change ({})'.format(massUnit))
   plt.grid()
   axs3.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs3.flatten()[reg]
   else:
      axs3.flatten()[reg].sharex(axX)

# Set up Figure 4: area change
fig4, axs4 = plt.subplots(nrow, ncol, figsize=(13, 11), num=4)
fig4.suptitle(f'Area change\n{runinfo}', fontsize=9)
for reg in range(nRegions):
   plt.sca(axs4.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('Area change (km^2)')
   plt.grid()
   axs4.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs4.flatten()[reg]
   else:
      axs4.flatten()[reg].sharex(axX)


# Set up Figure 5
fig5, axs5 = plt.subplots(2,1, figsize=(13, 11), num=5)
fig5.suptitle(f'regional contributions\n{runinfo}', fontsize=9)
mnTot=0.0
sigTot = 0.0
for reg in range(nRegions):
    if rNamesOrig[reg] in ISMIP6basinInfo:
        [mn, sig] = ISMIP6basinInfo[rNamesOrig[reg]]['net']
        mnTot += mn
        sigTot += sig**2

sigTot = sigTot**0.5
label = f'net obs ({get_dataset_label("mass_balance", selected_mass_balance)})'
axs5.flatten()[0].fill_between(yr, yr*(mnTot-sigTot), yr*(mnTot+sigTot),
                               color='k', alpha=0.2, label=label)
axs5.flatten()[1].fill_between(yr, yr*(mnTot-sigTot), yr*(mnTot+sigTot),
                               color='k', alpha=0.2, label=label)

plt.sca(axs5.flatten()[0])
plt.xlabel('Year')
plt.ylabel('Mass change (Gt)')
plt.grid()
plt.sca(axs5.flatten()[1])
plt.xlabel('Year')
plt.ylabel('VAF mass change (Gt)')
plt.grid()


# Set up Figure 6: melt rate vs obs
fig6, axs6 = plt.subplots(nrow, ncol, figsize=(13, 11), num=6)
fig6.suptitle(f'Ice-shelf melt rate\n{runinfo}', fontsize=9)
for reg in range(nRegions):
   plt.sca(axs6.flatten()[reg])
   plt.xlabel('Year')
   plt.ylabel('Ice-shelf melt rate (Gt/yr)')
   plt.grid()
   axs6.flatten()[reg].set_title(rNames[reg])
   if reg == 0:
      axX = axs6.flatten()[reg]
   else:
      axs6.flatten()[reg].sharex(axX)
   if rNamesOrig[reg] in ISMIP6basinInfo:
       # Get the melt data
       melt_data = ISMIP6basinInfo[rNamesOrig[reg]]['shelfMelt'][0]
       melt_mean, melt_unc = melt_data[0], melt_data[1]

       # Get label from selected dataset
       label = f'melt obs ({get_dataset_label("shelf_melt", selected_melt_dataset)})'

       if melt_unc is not None:
           # Plot uncertainty band + center line
           axs6.flatten()[reg].fill_between(
               yr,
               np.ones(yr.shape) * (melt_mean - melt_unc),
               np.ones(yr.shape) * (melt_mean + melt_unc),
               color='k', alpha=0.2, label=label
           )
           axs6.flatten()[reg].plot(
               yr, np.ones(yr.shape) * melt_mean,
               color='k', linestyle='-', linewidth=1.5
           )
       else:
           # Plot line only (no uncertainty)
           axs6.flatten()[reg].plot(
               yr, np.ones(yr.shape) * melt_mean,
               color='k', linestyle='-', linewidth=1.5, label=label
           )

# Set up unit conversion factors to be used when reading variables
if options.units == "m3":
    volUnitFactor = 1.0
    massUnitFactor = 1.0 / rhoi
elif options.units == "kg":
    volUnitFactor = rhoi
    massUnitFactor = 1.0
elif options.units == "Gt":
    volUnitFactor = rhoi / 1.0e12
    massUnitFactor = 1.0 / 1.0e12
else:
    sys.exit("ERROR: Unknown unit specified")


def plotStat(fname, sty, addToLegend=False):
    print("Reading and plotting file: {}".format(fname))

    name = fname

    f = Dataset(fname,'r')
    yr = f.variables['daysSinceStart'][:]/365.0
    dt = f.variables['deltat'][:]/(3600.0*24.0*365.0) # in yr
    dtnR = np.tile(dt.reshape(len(dt),1), (1,nRegions))  # repeated per region with dim of nt,nRegions
    nRegionsLocal = len(f.dimensions['nRegions'])
    if nRegionsLocal != nRegions:
        sys.exit(f"ERROR: Number of regions in file {fname} does not match number of regions in first input file!")

    # Fig 1: summary plot
    vol = f.variables['regionalIceVolume'][:] * volUnitFactor
    lbl ='total' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs1.flatten()[r].plot(yr, vol[:,r] - vol[0,r], label=lbl, linestyle=sty, color='k')

    VAF = f.variables['regionalVolumeAboveFloatation'][:] * volUnitFactor
    VAF = VAF[:,:] - VAF[0,:]
    lbl ='VAF' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs1.flatten()[r].plot(yr, VAF[:,r] - VAF[0,r], label=lbl, linestyle=sty, color='m')

    volGround = f.variables['regionalGroundedIceVolume'][:] * volUnitFactor
    volGround = volGround[:,:] - volGround[0,:]
    lbl ='grd' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs1.flatten()[r].plot(yr, volGround[:,r] - volGround[0,r], label=lbl, linestyle=sty, color='b')

    volFloat = f.variables['regionalFloatingIceVolume'][:] * volUnitFactor
    volFloat = volFloat[:,:] - volFloat[0,:]
    lbl ='flt' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs1.flatten()[r].plot(yr, volFloat[:,r] - volFloat[0,r], label=lbl, linestyle=sty, color='g')


    # Fig 2: Grd MB ------------
    lbl ='vol chg' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, volGround[:,r] - volGround[0,r], label=lbl, linestyle=sty, color='k', linewidth=2)

    grdSMB = f.variables['regionalSumGroundedSfcMassBal'][:] * massUnitFactor
    cumGrdSMB = np.cumsum(grdSMB*dtnR, axis=0)
    lbl ='SMB' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, cumGrdSMB[:,r], label=lbl, linestyle=sty, color='b')

    grdBMB = f.variables['regionalSumGroundedBasalMassBal'][:] * massUnitFactor
    cumGrdBMB = np.cumsum(grdBMB*dtnR, axis=0)
    lbl ='BMB' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, cumGrdBMB[:,r], label=lbl, linestyle=sty, color='c', lw=0.5)

    if 'regionalSumCalvingFluxGrounded' in f.variables:
       grdCalv = f.variables['regionalSumCalvingFluxGrounded'][:] * massUnitFactor
    else:
       grdCalv = grdBMB * 0.0  # set to zero if the stats file is missing this field
    cumGrdCalv = np.cumsum(grdCalv*dtnR, axis=0)
    lbl ='Grounded calving' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, -1.0*cumGrdCalv[:,r], label=lbl, linestyle=sty, color='lime')

    FM = f.variables['regionalSumFaceMeltingFlux'][:] * massUnitFactor
    cumFM = np.cumsum(FM*dtnR, axis=0)
    lbl ='FM' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, -1.0*cumFM[:,r], label=lbl, linestyle=sty, color='m')

    GLflux = f.variables['regionalSumGroundingLineFlux'][:] * massUnitFactor
    cumGLflux = np.cumsum(GLflux*dtnR, axis=0)
    lbl ='GL flux' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, -1.0*cumGLflux[:,r], label=lbl, linestyle=sty, color='g')

    GLMigflux = f.variables['regionalSumGroundingLineMigrationFlux'][:] * massUnitFactor
    cumGLMigflux = np.cumsum(GLMigflux*dtnR, axis=0)
    lbl ='GL mig flux' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, -1.0*cumGLMigflux[:,r], label=lbl, linestyle=sty, color='y')

    # sum of components
    grdSum = grdSMB + grdBMB - grdCalv - FM - GLflux - GLMigflux # note negative sign on two GL terms - they are both positive grounded to floating
    cumGrdSum = np.cumsum(grdSum*dtnR, axis=0)
    lbl ='sum' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs2.flatten()[r].plot(yr, cumGrdSum[:,r], label=lbl, linestyle=sty, color='hotpink', linewidth=0.75)

    grdSum2 = grdSum + GLMigflux  # version with migration flux removed - note the sign convention
    cumGrdSum2 = np.cumsum(grdSum2*dtnR, axis=0)
    lbl ='sum, no GLmig' if addToLegend else '_nolegend_'
    for r in range(nRegions):
        axs2.flatten()[r].plot(yr, cumGrdSum2[:,r], label=lbl, linestyle=':', color='hotpink', linewidth=0.75)


    # Fig 3: Flt MB ---------------
    lbl ='vol chg' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, volFloat[:,r] - volFloat[0,r], label=lbl, linestyle=sty, color='k', linewidth=2)

    fltSMB = f.variables['regionalSumFloatingSfcMassBal'][:] * massUnitFactor
    cumFltSMB = np.cumsum(fltSMB*dtnR, axis=0)
    lbl = 'SMB' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumFltSMB[:,r], label=lbl, linestyle=sty, color='b')

    lbl = 'GL flux' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumGLflux[:,r], label=lbl, linestyle=sty, color='g')

    lbl = 'GL mig flux' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumGLMigflux[:,r], label=lbl, linestyle=sty, color='y')

    clv = f.variables['regionalSumCalvingFlux'][:] * massUnitFactor
    cumClv = np.cumsum(clv*dtnR, axis=0)
    lbl = 'calving' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, -1.0*cumClv[:,r], label=lbl, linestyle=sty, color='m', linewidth=1)

    BMB = f.variables['regionalSumFloatingBasalMassBal'][:] * massUnitFactor
    cumBMB = np.cumsum(BMB*dtnR, axis=0)
    lbl = 'BMB' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumBMB[:,r], label=lbl, linestyle=sty, color='r', linewidth=1)

    # sum of components
    fltSum = fltSMB + GLflux + GLMigflux - clv + BMB
    cumFltSum = np.cumsum(fltSum*dtnR, axis=0)
    lbl = 'sum' if addToLegend else '_nolegend_'
    for r in range(nRegions):
       axs3.flatten()[r].plot(yr, cumFltSum[:,r], label=lbl, linestyle=sty, color='hotpink', linewidth=0.75)
    fltSum2 = fltSMB + GLflux - clv + BMB
    cumFltSum2 = np.cumsum(fltSum2*dtnR, axis=0)
    lbl = 'sum, no GLmig' if addToLegend else '_nolegend_'
    for r in range(nRegions):
        axs3.flatten()[r].plot(yr, cumFltSum2[:,r], label=lbl, linestyle=':', color='hotpink', linewidth=0.75)


    # Fig 4: area change  ---------------

    areaTot = f.variables['regionalIceArea'][:]/1000.0**2
    areaGrd = f.variables['regionalGroundedIceArea'][:]/1000.0**2
    areaFlt = f.variables['regionalFloatingIceArea'][:]/1000.0**2
    for r in range(nRegions):
        axs4.flatten()[r].plot(yr, areaTot[:,r] - areaTot[0,r], label=("total area" if addToLegend else '_nolegend_'), linestyle=sty, color='k')
        axs4.flatten()[r].plot(yr, areaGrd[:,r] - areaGrd[0,r], label=("grd area" if addToLegend else '_nolegend_'), linestyle=sty, color='b')
        axs4.flatten()[r].plot(yr, areaFlt[:,r] - areaFlt[0,r], label=("flt area" if addToLegend else '_nolegend_'), linestyle=sty, color='g')

    # Fig. 5:  select global stats ---------
    for r in range(nRegions):
        if rNamesOrig[r] == 'ISMIP6BasinGH':
           indTG = r
           break
    axs5.flatten()[0].plot(yr, volGround.sum(axis=1), label='total', color='b', linestyle=sty)
    volGroundnoTG = np.delete(volGround, indTG, 1)
    axs5.flatten()[0].plot(yr, volGroundnoTG.sum(axis=1), label='no TG/PIG', color='c', linestyle=sty)
    axs5.flatten()[1].plot(yr, VAF.sum(axis=1), label='total', color='b', linestyle=sty)
    VAFnoTG = np.delete(VAF, indTG, 1)
    axs5.flatten()[1].plot(yr, VAFnoTG.sum(axis=1), label='no TG/PIG', color='c', linestyle=sty)

    # Fig. 6:  melt rates ---------
    for r in range(nRegions):
        axs6.flatten()[r].plot(yr, -BMB[:,r], label=("BMB" if addToLegend else '_nolegend_'), linestyle=sty, color='b')

    f.close()


plotStat(options.file1inName, sty='-', addToLegend=True)


if(options.file2inName):
    plotStat(options.file2inName, sty=':')

if(options.file3inName):
    plotStat(options.file3inName, sty='--')

if(options.file4inName):
    plotStat(options.file4inName, sty='-.')


axs1.flatten()[-1].legend(loc='best', prop={'size': 5})
axs2.flatten()[-1].legend(loc='best', prop={'size': 5})
axs3.flatten()[-1].legend(loc='best', prop={'size': 5})
axs4.flatten()[-1].legend(loc='best', prop={'size': 6})
axs5.flatten()[0].legend(loc='best', prop={'size': 6})
axs6.flatten()[-1].legend(loc='best', prop={'size': 6})

print("Generating plot.")
fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()
fig5.tight_layout()
fig6.tight_layout()
plt.show()
