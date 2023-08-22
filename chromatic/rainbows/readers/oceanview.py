"""
This module serves as a template for creating a new Rainbow
reader. If you want to add the ability to read chromatic
light curves from a new kind of file format, a good process
would be to do something like the following

    1. Copy this `template.py` file into a new file in the
    `readers/` directory, ideally with a name that's easy
    to recognize, such as `readers/oceanview.py` (assuming
    `oceanview` is the name of your format)

    2. Start by finding and replacing `oceanview` in this
    template with the name of your format.

    3. Edit the `from_oceanview` function so that it will
    load a chromatic light curve file in your format and,
    for some Rainbow object `rainbow`, populate at least:

        + self.timelike['time']
        + self.wavelike['wavelength']
        + self.fluxlike['flux']

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
    `from_oceanview` function to be accessible when people
    are trying to create new Rainbows. Add an `elif` statement
    to the `guess_reader` function that will help guess which
    reader to use from some aspect(s) of the filename.

    (This `guess_reader` function also accepts a `format=`
    keyword that allows the user to explicitly specify that
    the oceanview reader should be used.)

    5. Submit a pull request to the github repository for
    this package, so that other folks can use your handy
    new reader too!
"""

# import the general list of packages
from ...imports import *

# define list of the only things that will show up in imports
__all__ = ["from_oceanview"]


# set the filename we want to open
def read_oceanview_spectrum(filename):
    """
    Read an ascii spectrum file from OceanView into
    a `Spectrum` objects from `rainbow-connection`.

    Parameters
    ----------
    filename : string
        The filename to load.

    """

    # create empty dictionary
    metadata = {}

    # print(f'Trying to read {filename}.')

    # open the file as a plain text file
    with open(filename) as f:

        # set up a search for the line indicating the start of the arrays
        i_first_line_of_data = 0
        this_line = ""
        string_to_look_for = ">>>>>Begin Spectral Data<<<<<"

        # print('The metadata lines found in the file header include:')

        # loop over all lines in the header until that line is found
        while string_to_look_for not in this_line:

            # read one line of the text file into a string
            this_line = f.readline()

            # advance a counter indicating we want to skip that line
            i_first_line_of_data = i_first_line_of_data + 1

            # include a break to avoid an endless loop
            if i_first_line_of_data > 10000:
                print(
                    f"""
                In loading the file...

                '{filename}'

                ...we tried to find a line that says

                '{string_to_look_for}'

                but didn't find it anywhere. Are you
                sure this is an ascii spectrograph
                file saved by OceanView?
                """
                )
                break

            # extract metadata from the header lines
            delimiter = ": "
            if delimiter in this_line:
                # split the data into keys (= labels) and values (= data)
                key, value = this_line.split(delimiter)

                # store the key and value in a dictionary
                metadata[key] = value.strip()

                # print out that line of metadata
                # print(f'{key:>40} = {metadata[key]}')

    # print(f'All rows after line {i_first_line_of_data} will be assumed to be wavelengths/fluxes.')

    # read that file into a table
    table = ascii.read(
        filename,  # which text file to load?
        data_start=i_first_line_of_data,  # where to start reading data
        names=["wavelength", "flux"],
    )  # what to call the columns

    # use the column names to extract arrays from the table
    w = table["wavelength"]
    f = table["flux"]

    # initialize with the wavelength and flux arrays
    wavelength = w * u.nm
    flux = f * u.photon

    return wavelength, flux, metadata


def from_oceanview(self, filepath):
    """
    Populate a Rainbow from a file in the oceanview format.

    Parameters
    ----------

    rainbow : Rainbow
        The object to be populated. This function is meant
        to be called during the initialization of a Rainbow
        object. You can/should assume that the `rainbow` object
        being passed will already contain the four core
        dictionaries `timelike`, `wavelike`, `fluxlike`, and
        `metadata`. This function should, at a minimum, add
        the following items
            + `self.timelike['time']`
            + `self.wavelike['wavelength']`
            + `self.fluxlike['flux']`
        and optionally, additional entries like
            + `self.metadata['some-useful-parameter']`
            + `self.fluxlike['uncertainty']`
            + `self.fluxlike['ok']`

    filepath : str
        The path to the file to load.
    """

    # read in your file, however you like
    files = expand_filenames(filepath)
    wfm_spectra = [read_oceanview_spectrum(f) for f in tqdm(files)]

    # populate a 1D array of wavelengths (with astropy units of length)
    self.wavelike["wavelength"] = wfm_spectra[0][0]

    # populate a 2D (row = wavelength, col = time) array of fluxes
    flux = u.Quantity([wfm[1] for wfm in wfm_spectra]).T
    self.fluxlike["flux"] = flux
    self.fluxlike["uncertainty"] = np.sqrt(flux.to_value("ph")) * u.ph

    # populate a 1D array of times (with astropy units of time)
    self.timelike["time"] = (
        np.arange(flux.shape[1])
        * float(wfm_spectra[0][2]["Integration Time (sec)"])
        * u.s
    )
