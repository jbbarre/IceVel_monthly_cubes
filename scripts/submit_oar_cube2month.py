import os
import sys
import time
import glob
import gc
from pathlib import Path

import numpy as np
import xarray as xr

from utils.logger import get_logger

from cube_helpers import CubeCalibrationHelper

from cube_class import Cube
from cube2monthly_map_using_linear_regression import (
    cube2monthly_map_using_linear_regression,
    check_vy_for_s2,
)



def load_all_cubes(
    file_pattern: str, rootdir: str = "/summer/ice_speed/surface_flow_velocity/"
):

    logger("Start loading cubes")

    search_patterns = [
        ("SENTINEL2/2018/10d", ""),  # Base S2 cube
        ("RADARSAT-2_", range(2016, 2022)),
        ("LANDSAT", range(2013, 2022), range(16, 416, 16)),
        ("SENTINEL1", range(2014, 2022), range(6, 30, 6)),
        ("SENTINEL2", range(2016, 2022), range(5, 405, 5)),
    ]

    cube_list = []

    # Load base S2 cube first
    base_path = Path(rootdir) / "ANTARCTICA/SENTINEL2/2018/10d/MOSAIC/cubes"
    base_files = list(base_path.glob(file_pattern))
    if not base_files:
        raise FileNotFoundError(f"No base cube found for {file_pattern} at {base_path}")
    logger(f"Loading base cube: {base_files[0]}")
    cube_list.append(Cube.from_netcdf(base_files[0]))

    # Load all other cubes
    for source in search_patterns[1:]:
        satellite, *ranges = source
        logger(f"Loading {satellite}...")

        if satellite.startswith("RADARSAT-2_"):
            for yr in ranges[0]:
                path = Path(rootdir) / f"ANTARCTICA/RADARSAT-2_{yr}/MOSAIC/cubes"
                files = list(path.glob(file_pattern))
                for f in files:
                    logger(f" -> {f}")
                    cube_list.append(Cube.from_netcdf(f))

        elif satellite == "LANDSAT":
            for yr in ranges[0]:
                for cy in ranges[1]:
                    path = Path(rootdir) / f"ANTARCTICA/LANDSAT/{yr}/{cy}d/MOSAIC/cubes"
                    files = list(path.glob(file_pattern))
                    for f in files:
                        logger(f" -> {f}")
                        cube_list.append(Cube.from_netcdf(f))

        elif satellite == "SENTINEL1":
            for yr in ranges[0]:
                for cy in ranges[1]:
                    path = (
                        Path(rootdir) / f"ANTARCTICA/SENTINEL1/{yr}/{cy}d/MOSAIC/cubes"
                    )
                    files = list(path.glob(file_pattern))
                    for f in files:
                        logger(f" -> {f}")
                        cube_list.append(Cube.from_netcdf(f))

        elif satellite == "SENTINEL2":
            for yr in ranges[0]:
                for cy in ranges[1]:
                    path = (
                        Path(rootdir) / f"ANTARCTICA/SENTINEL2/{yr}/{cy}d/MOSAIC/cubes"
                    )
                    files = list(path.glob(file_pattern))
                    for f in files:
                        logger(f" -> {f}")
                        cube_list.append(Cube.from_netcdf(f))

    if not cube_list:
        raise RuntimeError("No cubes were loaded. Check your file pattern.")

    # Concatenate all cubes along time axis
    logger(f"Concatenating {len(cube_list)} cubes...")

    # TODO. add help to Cube : describe() and sieve_empty_layers()
    # c.describe()
    # c.sieve_empty_layers(sieve_n=1000)
    # c.describe()

    vx_all = np.concatenate([c.vx for c in cube_list], axis=-1)
    vy_all = np.concatenate([c.vy for c in cube_list], axis=-1)
    errx_all = np.concatenate([c.errx for c in cube_list])
    erry_all = np.concatenate([c.erry for c in cube_list])
    t_all = np.concatenate([c.t for c in cube_list])

    combined_cube = Cube(vx_all, vy_all, errx_all, erry_all, t_all)
    combined_cube.mask_invalid()

    return combined_cube


def submit_cube2monthly(i: str, j: str):

    logger("Starting the process")
    start_time = time.time()
    cube_name_pattern = "x{0:05d}_y{1:05d}".format(int(i), int(j))

    output_files = [
        f"cube_monthly_{cube_name_pattern}_yr{year}.nc" for year in range(2015, 2022)
    ]
    if all(Path(f).exists() for f in output_files):
        logger("All monthly files already exist. Skipping.")
        return

    file_pattern = "c_" + cube_name_pattern + "_post0450_yr*nc"

    cube = load_all_cubes(file_pattern)

    logger(" . Check vy for S2")
    cube = check_vy_for_s2(cube)

    logger(" . Start calibration")
    calibrated_cube = CubeCalibrationHelper.calibrate_using_low_speed(cube)
    
    # Free up memory
    cube = None
    gc.collect()
    # TODO : continue from here 
    # TODO: 04/18/2025
    # TODO:vx and vy reordering (x,y,z) may cause trouble 
    
    for year in range(2015, 2022):
        out_file = f"cube_monthly_{cube_name_pattern}_yr{year}.nc"
        if not Path(out_file).exists():
            logger(f" . Generating {out_file}...")
            ds = cube2monthly_map_using_linear_regression(
                calibrated_cube, year, year, deltamax=36
            )

            comp = dict(zlib=True, complevel=5)
            encoding = {var: comp for var in ds.data_vars}

            ds.to_netcdf(out_file, encoding=encoding)
            ds.close()
            gc.collect()

    logger(f"Finished processing cube in {time.time() - start_time:.2f} seconds.")


if __name__ == "__main__":

    logger = get_logger("IceVel_monthlycubes")

    logger.info("Starting processing...")
    
    if len(sys.argv) == 3:
        submit_cube2monthly(sys.argv[1], sys.argv[2])
    else:

        submit_cube2monthly("02695", "07105")
