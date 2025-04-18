from typing import Optional, Union
import numpy as np
import xarray as xr
from pathlib import Path

try:
    import dask.array as da
except ImportError:
    da = None

ArrayLike = Union[np.ndarray, "dask.array.Array"]

class Cube:
    __slots__ = ("vx", "vy", "errx", "erry", "t", "meta")

    def __init__(
        self,
        vx: ArrayLike,
        vy: ArrayLike,
        errx: ArrayLike,
        erry: ArrayLike,
        t: ArrayLike,
        meta: Optional[dict] = None,
    ):
        self.vx = vx
        self.vy = vy
        self.errx = errx
        self.erry = erry
        self.t = t
        self.meta = meta or {}
        self._validate()

    @classmethod
    def from_netcdf(cls, filepath: Union[str, Path], chunks: Optional[dict] = None) -> "Cube":
        ds = xr.open_dataset(filepath, chunks=chunks)
        return cls(
            vx=ds["vx"].data,
            vy=ds["vy"].data,
            errx=ds["errx"].data,
            erry=ds["erry"].data,
            t=ds["t"].data,
            meta=dict(ds.attrs),
        )

    def to_netcdf(self, filepath: Union[str, Path]):
        ds = xr.Dataset(
            {
                "vx": ("x", "y", "time", self.vx),
                "vy": ("x", "y", "time", self.vy),
                "errx": ("time", self.errx),
                "erry": ("time", self.erry),
                "t": ("time", self.t),
            },
            attrs=self.meta,
        )
        ds.to_netcdf(filepath)

    def copy(self) -> "Cube":
        return Cube(
            vx=self.vx.copy(),
            vy=self.vy.copy(),
            errx=self.errx.copy(),
            erry=self.erry.copy(),
            t=self.t.copy(),
            meta=self.meta.copy(),
        )

    def mask_invalid(self, threshold: float = 1e5):
        invalid_mask = (
            np.isnan(self.vx) | np.isnan(self.vy) |
            np.isinf(self.vx) | np.isinf(self.vy) |
            (np.abs(self.vx) > threshold) | (np.abs(self.vy) > threshold) |
            ((self.vx == 0) & (self.vy == 0))
        )
        if hasattr(self.vx, 'where'):
            self.vx = self.vx.where(~invalid_mask)
            self.vy = self.vy.where(~invalid_mask)
        else:
            self.vx[invalid_mask] = np.nan
            self.vy[invalid_mask] = np.nan

    def chunk_cube(self, chunks: dict = {"x": -1, "y": -1, "time": 128}):
        """
        Convert vx and vy to dask arrays, chunked along time.
        Call AFTER concatenation.
        """
        self.vx = xr.DataArray(self.vx, dims=("x", "y", "time")).chunk(chunks).data
        self.vy = xr.DataArray(self.vy, dims=("x", "y", "time")).chunk(chunks).data

    def compute_weighted_mean(self) -> xr.Dataset:
        wx = 1.0 / (self.errx ** 2)
        wy = 1.0 / (self.erry ** 2)

        wx_expand = xr.DataArray(wx, dims=("time",))
        wy_expand = xr.DataArray(wy, dims=("time",))

        vx_weighted_sum = (self.vx * wx_expand).sum(dim="time", skipna=True)
        vx_weighted_total = wx_expand.sum(dim="time", skipna=True)
        vy_weighted_sum = (self.vy * wy_expand).sum(dim="time", skipna=True)
        vy_weighted_total = wy_expand.sum(dim="time", skipna=True)

        weighted_vx = vx_weighted_sum / vx_weighted_total
        weighted_vy = vy_weighted_sum / vy_weighted_total

        std_vx = np.sqrt(((self.vx - weighted_vx.expand_dims(time=self.t.shape[0])) ** 2 * wx_expand).sum(dim="time") / (vx_weighted_total))
        std_vy = np.sqrt(((self.vy - weighted_vy.expand_dims(time=self.t.shape[0])) ** 2 * wy_expand).sum(dim="time") / (vy_weighted_total))

        return xr.Dataset({
            "weighted_vx": weighted_vx,
            "weighted_vy": weighted_vy,
            "std_vx": std_vx,
            "std_vy": std_vy,
        })

    def compute_median_velocity(self) -> xr.Dataset:
        median_vx = self.vx.median(dim="time", skipna=True)
        median_vy = self.vy.median(dim="time", skipna=True)

        return xr.Dataset({
            "median_vx": median_vx,
            "median_vy": median_vy,
        })

    def apply_custom_mask(self, mask: ArrayLike):
        if hasattr(self.vx, 'where'):
            self.vx = self.vx.where(mask)
            self.vy = self.vy.where(mask)
        else:
            self.vx[~mask] = np.nan
            self.vy[~mask] = np.nan

    def _validate(self):
        assert self.vx.shape == self.vy.shape, "vx and vy must have same shape"
        assert self.vx.ndim == 3, "vx must be 3D (x, y, time)"
        assert self.errx.ndim == 1, "errx must be 1D"
        assert self.erry.ndim == 1, "erry must be 1D"
        assert self.t.ndim == 1, "t must be 1D"
        assert self.vx.shape[-1] == self.errx.shape[0] == self.erry.shape[0] == self.t.shape[0], "Mismatch along time dimension"

    @property
    def is_chunked(self) -> bool:
        return hasattr(self.vx, "chunks")
