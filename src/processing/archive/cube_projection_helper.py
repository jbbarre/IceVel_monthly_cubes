from typing import Dict
from osgeo import osr

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
