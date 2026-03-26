#!/usr/bin/env python3
import argparse
import xarray as xr


TIME_DIM = "time_counter"
DEPTH_DIM = "deptht"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute monthly climatology over all years for variables "
            "without the depth dimension"
        )
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Input NetCDF file"
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Output NetCDF file (monthly climatology)"
    )
    parser.add_argument(
        "nb_year",
        type=int,
        nargs="?",
        default=-1,
        help="Number of years from the end over which to compute the mean (default: entire period)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    ds = xr.open_dataset(args.input_file)

    # Select last N years if requested
    if args.nb_year > 0:
        ntime = ds.dims["time_counter"]

        if args.nb_year > ntime:
            raise ValueError(
                f"Requested {args.nb_year} years but dataset only contains {ntime} time steps."
            )

        ds = ds.isel(time_counter=slice(-args.nb_year*12, None))
        print(f"Averaging over last {args.nb_year} years")
    else:
        print("Averaging over full time period")

    vars_2d = [
        v for v in ds.data_vars
        if TIME_DIM in ds[v].dims
        and DEPTH_DIM not in ds[v].dims
    ]

    ds = ds.assign_coords(
        month=ds[TIME_DIM].dt.month
    )

    ds_clim = (
        ds[vars_2d]
        .groupby("month")
        .mean(dim=TIME_DIM, skipna=True)
    )

    for coord in ("nav_lat", "nav_lon"):
        if coord in ds:
            ds_clim[coord] = ds[coord]

    ds_clim.attrs = ds.attrs.copy()
    ds_clim.to_netcdf(args.output_file)

    print(f"File created : {args.output_file}")


if __name__ == "__main__":
    main()

