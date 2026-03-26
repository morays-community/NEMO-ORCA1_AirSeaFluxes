#!/usr/bin/env python3
import argparse
import xarray as xr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute temporal mean of all time-dependent fields in a NetCDF file"
    )

    parser.add_argument(
        "input_file",
        type=str,
        help="Input NetCDF file"
    )

    parser.add_argument(
        "output_file",
        type=str,
        help="Output NetCDF file (time mean)"
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

    # Open dataset (lazy loading, safe for large files)
    ds = xr.open_dataset(args.input_file)

    # Select last N years if requested
    if args.nb_year > 0:
        ntime = ds.dims["time_counter"]

        if args.nb_year > ntime:
            raise ValueError(
                f"Requested {args.nb_year} years but dataset only contains {ntime} time steps."
            )

        if args.nb_year > 0:
            ds = ds.isel(time_counter=slice(-args.nb_year, None))  # dernières années
        elif args.nb_year < 0:
            ds = ds.isel(time_counter=slice(0, -args.nb_year))    # premières années
        print(f"Averaging over last {args.nb_year} years")
    else:
        print("Averaging over full time period")

    # Select variables that depend on time
    vars_with_time = [
        v for v in ds.data_vars
        if "time_counter" in ds[v].dims
    ]

    # Compute temporal mean
    ds_mean = ds[vars_with_time].mean(
        dim="time_counter",
        skipna=True
    )

    # Preserve horizontal coordinates if present
    for coord in ("nav_lat", "nav_lon"):
        if coord in ds:
            ds_mean[coord] = ds[coord]

    # Preserve vertical coordinate if present
    if "deptht" in ds:
        ds_mean["deptht"] = ds["deptht"]

    # Copy global attributes
    ds_mean.attrs = ds.attrs.copy()
    ds_mean.attrs["description"] = "Temporal mean"
    ds_mean.attrs["averaging"] = (
        f"Mean over last {args.nb_year} years"
        if args.nb_year > 0
        else "Mean over full time period"
    )

    # Write output file
    ds_mean.to_netcdf(args.output_file)

    print(f"Created file: {args.output_file}")


if __name__ == "__main__":
    main()

