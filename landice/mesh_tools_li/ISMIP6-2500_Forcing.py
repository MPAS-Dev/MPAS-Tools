import os 
import sys
import argparse
import numpy as np
import pandas as pd
import xarray as xr
from cftime import datetime
from mpas_tools.io import write_netcdf

def __xtime2cftime(xtime_str, format="%Y-%m-%d_%H:%M:%S"): 
    """Convert a single xtime value to a cftime datetime object
    """

    stripped_str = xtime_str.tobytes().decode("utf-8").strip().strip('\x00')
    
    return datetime.strptime(stripped_str, format)
    
def parse_xtime(ds): 
    """Parse `xtime` vaules to cftime datetimes and create new coordinate
       from parsed values, named `Time`
    """

    times = list(map(__xtime2cftime, ds.xtime.values))
    
    time_da = xr.DataArray(times, [('Time', times)])
    
    return ds.assign_coords({"Time":time_da})

def generate_samples(ds, sampling_start, sampling_end, rng, repeat_period=200):
    """Generate array of sample indices from `Time` dimension of `ds`
    """

    sampling_period = sampling_end - sampling_start

    # 0 index assumes `Time` coordinate is sorted
    forcing_start = int(ds.Time.dt.year[0])
    
    # get the offset for the indexes
    idx_offset = sampling_start - forcing_start
    
    samples = idx_offset + rng.integers(0, sampling_period, repeat_period)
        
    return samples

def extend_forcing(src_ds, sample_idxs): 
    """Create new "extended" dataset using the `sample_idxs`
    """
    
    # create new time index based on hardcoded (for now) input start and 
    # end years with an calendar occurance on the first of every year
    time_da = xr.cftime_range("2300-01-01", "2500-01-01",
                              freq="YS", inclusive="left", calendar="noleap")
    # conver the CFtimeindex to a left justified string
    xtime_da = time_da.strftime("%Y-%m-%d_%H:%M:%S").str.ljust(64)

    # create a "new" extended datatset with the xtime variable
    new_ds = xr.Dataset({"xtime": ("Time", xtime_da)})
    # convert xtime from string to character array
    new_ds["xtime"] = new_ds.xtime.astype("S")
    
    for var in src_ds:
        # do not copy xtime, as it was created above
        if var == "xtime": continue 
    
        # only sample variables w/ Time dimension
        if "Time" in src_ds[var].dims:
            new_ds[var] = src_ds[var].isel(Time=sample_idxs).drop_vars("Time")
        else: 
            new_ds[var] = src_ds[var].copy()
   
    return new_ds


def cli_parser(argv): 
    """Command line interfacr parser
    """
    parser = argparse.ArgumentParser(prog='ISMIP6 2500 Extensions')

    parser.add_argument('-i', '--input', type=str,
                        help="input forcing to sample from")
    parser.add_argument('-o', '--output_filename', type=str,
                        help="output filename of extended forcing")
    parser.add_argument('-s', '--seed', type=int, default=4727,
                        help="seed for random number generator")
    parser.add_argument('--sampling_start', type=int, default=2270,
                        help="start of sampling window in reference dataset")
    parser.add_argument('--sampling_end', type=int, default=2300,
                        help="end of sampling window in reference dataset")
    parser.add_argument('--repeat_period', type=int, default=200,
                        help="length of extended forcing window")
    
    args, _ = parser.parse_known_args(argv)
    
    assert os.path.exists(args.input)

    path = os.path.dirname(args.input)

    output_filename = os.path.join(path, args.output_filename)
    
    return (args.input, output_filename, args.seed,
            args.sampling_start, args.sampling_end, args.repeat_period)


if __name__ == "__main__": 
    
    # parse the command line arguments
    input_fp, output_fp, seed, start, end, period = cli_parser(sys.argv[1:]) 
    
    # open the reference dataset and parse the xtime variable
    ref_forcing = xr.open_dataset(input_fp)
    ref_forcing = parse_xtime(ref_forcing)
    
    # initialize the random number generator and create sample index array
    rng = np.random.default_rng(seed)
    sample_idxs = generate_samples(ref_forcing, start, end, rng, period)
    
    # find the sample years from the sample indices
    sample_years = ref_forcing.isel(Time=sample_idxs).Time.dt.year.values

    # generate the extended forcing file and write it to disk
    new_forcing = extend_forcing(ref_forcing, sample_idxs)
    write_netcdf(new_forcing, output_fp)

    print("\n" + "*"*75)
    print(f"ISMIP6 2500 forcing created by randomly sampling 2300 forcing "
          f"files from {start}-{end} \n"
          f"\nSampled file:  {input_fp}"
          f"\nExtened file: {output_fp}")
    print("*"*75)
