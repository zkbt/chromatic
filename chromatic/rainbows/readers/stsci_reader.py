from ...imports import *
import h5py as h5
from astropy.io import ascii
import os

os.chdir('/Users/willwaalkes/Desktop/PhD_Thesis/chromatic/chromatic/tests/example-extracted-datasets/eureka-extraction')


def read_stsci(filename):
    """
    Read STScI's output table (time-series spectra).
    """

    # load the event data
    data = ascii.read(filename)


    # pull out some variables
    t = data['bjdtdb']
    f = data['optspec']
    e = data['opterr']
    w = data['wave_1d']
    # this assumes the wavelength axis is the same for all exposures

    timelike = {}
    timelike["time"] = t * u.day  # TODO: check time units

    wavelike = {}
    wavelike["wavelength"] = w * u.micron  # TODO: check wavelength units

    fluxlike = {}
    fluxlike["flux"] = f.transpose()
    fluxlike["error"] = e.transpose()

    return wavelike, timelike, fluxlike
    # TO-DO: extract data properly, organize wavelengths to make a cube


def from_stsci(rainbow, filename, **kwargs):
    """
    Populate a Rainbow from a stsci pipeline output.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated. (This is intended
        to be used as a class method of Rainbow or
        of a class derived from Rainbow, as a way of
        initializing an object from files.)
    filename : str
        The path to the file.
    """

    # load the Eureka event
    wavelike, timelike, fluxlike = read_stsci(filename)

    # populate the rainbow
    rainbow._initialize_from_dictionaries(
        wavelike=wavelike, timelike=timelike, fluxlike=fluxlike
    )
