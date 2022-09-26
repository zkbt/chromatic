"""
    3. Edit the `to_xarray_raw_light_curves` function so that it will
    write out a Rainbow to whatever format you want.
"""

# import the general list of packages
from ...imports import *

import xarray as xr
import json
from astropy.utils.misc import JsonCustomEncoder

# define list of the only things that will show up in imports
__all__ = ["to_xarray_raw_light_curves"]

required_attrs = ["author", "contact", "code", "data_origin"]
required_coords = ["wavelength", "time"]
required_data_vars = ["raw_flux", "raw_flux_error"]

chromatic_to_ers = dict(
    flux="raw_flux",
    uncertainty="raw_flux_error",
    ok="quality_flag",
    wavelength_lower="start_wavelength",
    wavelength_upper="end_wavelength",
    # wavelength="central_wavelength",
    # time="time_flux",
)


def to_xarray_raw_light_curves(self, filepath, overwrite=True):
    """
    Write a Rainbow to a file in the xarray_raw_light_curves format.

    Parameters
    ----------
    self : Rainbow
        The object to be saved.

    filepath : str
        The path to the file to write.
    """

    # warn about missing metadata
    for k in required_attrs:
        if k not in self.metadata:
            cheerfully_suggest(
                f"""
            The required metadata keyword `{k}` was not found.
            Before saving, please set it with `rainbow.{k} = ?`
            """
            )

    self._make_sure_wavelength_edges_are_defined()

    # populate the attrs
    attrs_dict = {}
    for k, v in self.metadata.items():
        # if isinstance(v, dict):
        attrs_dict[k] = json.dumps(v, cls=JsonCustomEncoder)
        # else:
        #    attrs_dict[k] = v

    # populate the data_vars and coords
    data_vars_dict = {}
    coords_dict = {}

    def make_DataArray(key, **kw):
        quantity = self.get(key)

        try:
            data = quantity.value
            units = quantity.unit.to_string()
        except AttributeError:
            data = quantity
            units = ""
        attrs = self.metadata.get(f"metadata-for-{k}", {})
        attrs.update(units=units)

        return xr.DataArray(
            name=chromatic_to_ers.get(key, key), data=data, attrs=attrs, **kw
        )

    for k in ["time", "wavelength"]:
        da = make_DataArray(key=k, coords={k: self.get(k)}, dims=[k])
        coords_dict[da.name] = da

    # warn about missing data_vars
    for k in required_coords:
        if k not in coords_dict:
            cheerfully_suggest(
                f"""
            The required coord keyword `{k}` was not found.
            Before saving, please set it with `rainbow.{k} = ?`
            """
            )

    for k in self.fluxlike:
        da = make_DataArray(key=k, coords=coords_dict, dims=["wavelength", "time"])
        data_vars_dict[da.name] = da

    for k in self.timelike:
        if k in coords_dict:
            continue
        da = make_DataArray(
            key=k, coords={c: coords_dict[c] for c in ["time"]}, dims=["time"]
        )
        data_vars_dict[da.name] = da

    for k in self.wavelike:
        if k in coords_dict:
            continue
        da = make_DataArray(
            key=k,
            coords={c: coords_dict[c] for c in ["wavelength"]},
            dims=["wavelength"],
        )
        data_vars_dict[da.name] = da

    # populate the dataset
    ds = xr.Dataset(data_vars=data_vars_dict, coords=coords_dict, attrs=attrs_dict)

    # warn about missing data_vars
    for k in required_data_vars:
        if k not in ds:
            cheerfully_suggest(
                f"""
            The required data_var keyword `{k}` was not found.
            Before saving, please set it with `rainbow.{k} = ?`
            """
            )

    # save out to file
    if overwrite:
        try:
            os.remove(filepath)
        except FileNotFoundError:
            pass
    ds.to_netcdf(filepath, mode="w")  # , engine="h5netcdf", invalid_netcdf=True)
