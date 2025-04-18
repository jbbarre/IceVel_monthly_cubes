"""
cube_helpers.py

Helper classes and functions for Cube processing:
- Slicing helpers (vx_, vy_ access)
- Calibration helpers
- Masking utilities
- Future extensions
"""

import numpy as np
from typing import Dict
from osgeo import osr
import logging

from cube_class import Cube

logger = logging.getLogger("IceVel_monthlycubes")


class CubeSlicingHelper:
    """Methods to slice VX and VY from  the entire cube

    - Full time series at a pixel (i,j): vx[i, j, :]
    - Full velocity map at time z: vx[:, :, z]
    - Time series along vertical x=i: vx[i, :, :]
    - Time series along horizontal y=j: vx[:, j, :]
    - Reorder cube to (time, x, y): np.moveaxis(vx, -1, 0)

    Note:
    - Cube data shape: (x, y, time)
    - vx and vy arrays are organized similarly.
    
    """

    @staticmethod
    def vx_(cube:Cube, i=None, j=None):
        vx = np.asarray(cube.vx)

        if i is None and j is None:
            return np.moveaxis(vx, -1, 0)  # Reorganize (x, y, time) --> (time, x, y)
        elif i is not None and j is None:
            return np.moveaxis(vx[i, :, :], -1, 0)[:, np.newaxis, :]
        elif i is None and j is not None:
            return np.moveaxis(vx[:, j, :], -1, 0)[:, :, np.newaxis]
        elif i is not None and j is not None:
            return np.moveaxis(vx[i, j, :], -1, 0)

    @staticmethod
    def vy_(cube:Cube, i=None, j=None):
        vy = np.asarray(cube.vy)

        if i is None and j is None:
            return np.moveaxis(vy, -1, 0) # Reorganize (x, y, time) --> (time, x, y)
        elif i is not None and j is None:
            return np.moveaxis(vy[i, :, :], -1, 0)[:, np.newaxis, :]
        elif i is None and j is not None:
            return np.moveaxis(vy[:, j, :], -1, 0)[:, :, np.newaxis]
        elif i is not None and j is not None:
            return np.moveaxis(vy[i, j, :], -1, 0)



class CubeCalibrationHelper:
    """Methods for calibrating cubes using low-speed fields."""

    @staticmethod
    def calibrate_using_low_speed(cube: Cube, inplace=False, low_speed_thresh=100, intermediate_speed_thresh=200):
        """
        Adjust Cube velocity fields based on low-speed median reference.
        """
        if inplace:
            c2 = cube
        else:
            c2 = cube.copy()


        vx_data = np.asarray(cube.vx)
        vy_data = np.asarray(cube.vy)

        vx_med = np.nanmedian(vx_data, axis=-1)
        vy_med = np.nanmedian(vy_data, axis=-1)
        v_med = np.sqrt(vx_med**2 + vy_med**2)

        for idx in range(vx_data.shape[-1]):
            vx_slice = np.asarray(c2.vx[..., idx])
            vy_slice = np.asarray(c2.vy[..., idx])

            mask_low = (v_med < low_speed_thresh)
            mask_intermediate = (v_med < intermediate_speed_thresh)

            if np.count_nonzero(mask_low) > 100:
                difx = np.nanmedian(vx_slice[mask_low] - vx_med[mask_low])
                if np.abs(difx) > 10:
                    vx_slice -= difx

                dify = np.nanmedian(vy_slice[mask_low] - vy_med[mask_low])
                if np.abs(dify) > 10:
                    vy_slice -= dify

            elif np.count_nonzero(mask_intermediate) > 100:
                difx = np.nanmedian(vx_slice[mask_intermediate] - vx_med[mask_intermediate])
                if np.abs(difx) > 30:
                    vx_slice -= difx

                dify = np.nanmedian(vy_slice[mask_low] - vy_med[mask_low])  # typo kept for logical compatibility
                if np.abs(dify) > 30:
                    vy_slice -= dify

            # Save corrections
            c2.vx[..., idx] = vx_slice
            c2.vy[..., idx] = vy_slice

        return c2
    


class CubeProjectionHelper:
    """
    Helper class to handle spatial reference system (SRS) operations for a Cube.
    Builds an OSR SpatialReference object from Cube metadata.
    """

    def __init__(self, meta: Dict[str, str]):
        """
        Initialize the CubeProjectionHelper.

        Args:
            meta: A dictionary containing projection metadata from a Cube.
        """
        self.meta = meta
        self.srs = self._build_spatial_ref()

    def _build_spatial_ref(self) -> osr.SpatialReference:
        """
        Build a SpatialReference object from metadata.

        Returns:
            osr.SpatialReference: The built spatial reference.

        Raises:
            ValueError: If no valid projection information is found.
        """
        srs = osr.SpatialReference()

        if "proj4" in self.meta and self.meta["proj4"]:
            srs.ImportFromProj4(self.meta["proj4"])
        elif "projection" in self.meta and "EPSG:" in self.meta["projection"]:
            epsg_code = int(self.meta["projection"].replace("EPSG:", ""))
            srs.ImportFromEPSG(epsg_code)
        else:
            raise ValueError("No valid projection information found in metadata (need 'proj4' or 'projection').")

        return srs

    def export_wkt(self) -> str:
        """
        Export the spatial reference as a WKT string.

        Returns:
            str: WKT string representation.
        """
        return self.srs.ExportToWkt()

    def export_proj4(self) -> str:
        """
        Export the spatial reference as a proj4 string.

        Returns:
            str: proj4 string representation.
        """
        return self.srs.ExportToProj4()

    def get_epsg_code(self) -> int:
        """
        Attempt to retrieve the EPSG code from the SRS.

        Returns:
            int: EPSG code if available.

        Raises:
            ValueError: If EPSG code is not found.
        """
        authority_name = self.srs.GetAuthorityName(None)
        authority_code = self.srs.GetAuthorityCode(None)

        if authority_name == 'EPSG' and authority_code is not None:
            return int(authority_code)
        else:
            raise ValueError("EPSG code not found in spatial reference.")

