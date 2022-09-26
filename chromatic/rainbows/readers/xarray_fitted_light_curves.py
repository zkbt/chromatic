"""

    3. Edit the `from_xarray_fitted_light_curves` function so that it will
    load a chromatic light curve file in your format and,
    for some Rainbow object `rainbow`, populate at least:

        + rainbow.timelike['time']
        + rainbow.wavelike['wavelength']
        + rainbow.fluxlike['flux']

    You'll need to replace the cartoon functions on each
    line with the actual code needed to load your file.

    (This template assumes that only one file needs to be
    loaded. If you need to load multiple segments, or each
    time point is stored in its own file or something, then
    check out `stsci.py` for an example of loading and
    stitching together multiple input files. You'll probably
    want to change `filepath` to accept a glob-friendly string
    like `my-neato-formatted-files-*.npy` or some such.)

    4. Edit the `readers/__init__.py` file to import your
    `from_xarray_fitted_light_curves` function to be accessible when people
    are trying to create new Rainbows. Add an `elif` statement
    to the `guess_reader` function that will help guess which
    reader to use from some aspect(s) of the filename.

    (This `guess_reader` function also accepts a `format=`
    keyword that allows the user to explicitly specify that
    the xarray_fitted_light_curves reader should be used.)

    5. Submit a pull request to the github repository for
    this package, so that other folks can use your handy
    new reader too!
"""

# import the general list of packages
from ...imports import *
from ..writers.xarray_fitted_light_curves import xr, json, chromatic_to_ers

ers_to_chromatic = {v: k for k, v in chromatic_to_ers.items()}

# define list of the only things that will show up in imports
__all__ = ["from_xarray_fitted_light_curves"]


def from_xarray_fitted_light_curves(self, filepath):
    """
    Populate a Rainbow from a file in the xarray_fitted_light_curves format.

    Parameters
    ----------

    self : self
        The object to be populated.

    filepath : str
        The path to the file to load.
    """

    import xarray as xr

    ds = xr.open_dataset(filepath)

    def make_Quantity(da):
        """
        Convert a data array into a chromatic quantity (with astropy units).

        """
        unit_string = da.attrs.get("units", "")
        if unit_string != "":
            unit = u.Unit(unit_string)
        else:
            unit = 1

        return da.data * unit

    self.wavelike["wavelength"] = make_Quantity(
        ds[ers_to_chromatic.get("wavelength", "wavelength")]
    )
    self.timelike["time"] = make_Quantity(ds[ers_to_chromatic.get("time", "time")])

    for key, da in ds.items():
        chromatic_key = ers_to_chromatic.get(key, key)
        self._put_array_in_right_dictionary(chromatic_key, make_Quantity(da))
        for k, v in da.attrs.items():
            if k != "units":
                metadata_key = f"metadata-for-{chromatic_key}"
                try:
                    self.metadata[metadata_key]
                except KeyError:
                    self.metadata[metadata_key] = dict()
                self.metadata[metadata_key][k] = v

    # kludge to convert to corrected flux
    k = "systematics_model"
    if k not in self.fluxlike:
        self.fluxlike[k] = self.flux / self.get("corrected_flux", 1)

    # kludge to convert to corrected flux
    # k = "corrected_flux_error"
    # self.fluxlike[k] = self.get(k, self.uncertainty / self.get("systematics_model", 1))

    for k, v in ds.attrs.items():

        self.metadata[k] = json.loads(v)
