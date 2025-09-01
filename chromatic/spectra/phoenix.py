from .library import *


class Library_PHOENIX_NextGen(Library):
    # a string naming this library
    _library_name = "phoenix+nextgen"

    def __init__(self, *args, **kwargs):
        """
        Initialize a PHOENIX model library to provide easy
        access to model stellar spectra at resolutions up
        from the PHOENIX NextGen grid (Husser et al. 2013).
        """

        super().__init__(*args, **kwargs)
        self._cache_label = f"{self._cache_label}-{self._library_name}"

    def _download_raw_data(self, metallicity=0.0, cache=True):
        """
        Make sure the raw data from the online PHOENIX database
        are downloaded to your local computer. (Most users
        shouldn't have to interact with this, unless you
        want to do something particularly fancy.)

        You must be connected to the internet for this to work.
        It may take quite a long time!

        Parameters
        ----------
        metallicity : float
            The stellar metallicity (= log10[metals/solar])
        cache : bool
            Once it's downloaded, should we keep it for next time?
        """

        # create a dictionary to store the local
        self._raw_local_paths = {}

        # where are the raw data located at online?
        self._raw_base_url = "ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS"

        # get the one shared wavelength array
        self._raw_wavelengths_filename = "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
        self._raw_wavelengths_url = "/".join(
            [self._raw_base_url, self._raw_wavelengths_filename]
        )
        cheerfully_suggest(
            f"""
        Downloading (or finding locally) the shared wavelength grid from
        {self._raw_wavelengths_url}
        """
        )
        self._raw_local_paths["wavelengths"] = download_file_with_warning(
            self._raw_wavelengths_url, pkgname=self._cache_label, cache=cache
        )

        # get the index of files in this metallicity's directory
        metallicity_string = self._stringify_metallicity(metallicity)
        self._raw_directory = f"PHOENIX-ACES-AGSS-COND-2011/Z{metallicity_string}"
        self._raw_directory_url = "/".join([self._raw_base_url, self._raw_directory])
        self._raw_local_paths[f"index-Z={metallicity_string}"] = (
            download_file_with_warning(
                self._raw_directory_url, pkgname=self._cache_label, cache=cache
            )
        )
        cheerfully_suggest(
            f"""
        Downloading (or finding locally) the index of files from
        {self._raw_directory_url}
        """
        )

        # get the complete list of spectrum files
        self._raw_spectrum_filenames = list(
            [
                x
                for x in ascii.read(
                    self._raw_local_paths[f"index-Z={metallicity_string}"]
                )
                .columns[-1]
                .data
                if ".fits" in x
            ]
        )
        self._raw_spectrum_urls = [
            "/".join([self._raw_base_url, self._raw_directory, f])
            for f in self._raw_spectrum_filenames
        ]
        N = len(self._raw_spectrum_urls)
        cheerfully_suggest(
            f"""
        Downloading (or finding locally) {N} very large files from
        {self._raw_directory_url}

        If the files aren't already downloaded,
        this might take an annoyingly long time!
        If it crashes due to a timeout,
        try restarting.
        """
        )
        self._current_raw_metallicity = metallicity
        self._raw_downloaded = {}
        for url, file in tqdm(
            zip(self._raw_spectrum_urls, self._raw_spectrum_filenames), leave=False
        ):
            self._raw_downloaded[file] = download_file_with_warning(
                url, pkgname=self._cache_label, cache=cache
            )

    def _load_raw_wavelength(self):
        """
        Load in the raw wavelength array.

        Returns
        -------
        wavelength : Quantity
            The wavelengths associated with this grid,
            with astropy units of microns.
        """

        # default to the preloaded raw wavelength
        try:
            return self._raw_wavelength

        # load the raw wavelengths and save them for next time
        except AttributeError:
            wavelength_filename = self._raw_local_paths["wavelengths"]
            hdu = fits.open(wavelength_filename)
            wavelength_without_unit = hdu[0].data
            wavelength_unit = u.Angstrom
            wavelength = wavelength_without_unit * wavelength_unit
            self._raw_wavelength = wavelength.to("micron")
            return self._raw_wavelength

    def _load_raw_spectrum(self, filename):
        """
        Load in a raw spectrum array.

        Parameters
        ----------
        filename : string
            The filename of the raw PHOENIX spectrum
        Returns
        -------
        spectrum : Quantity
            The spectrum, with astropy units of W/(m**2 nm)
        """
        hdus = fits.open(filename)
        flux_without_unit = hdus[0].data
        flux_unit = u.Unit("erg/(s * cm**2 * cm)")
        flux = flux_without_unit * flux_unit
        return flux.to("W/(m**2 nm)")

    def _stringify_metallicity(self, metallicity):
        """
        Convert a metallicity into a PHOENIX-style string.

        Parameters
        ----------
        metallicity : float
            [Fe/H]-style metallicity (= 0.0 for solar)

        Returns
        -------
        s : string
            The metallicity, as a PHOENIX-style string.
        """
        if metallicity <= 0:
            return f"-{np.abs(metallicity):03.1f}"
        else:
            return f"+{metallicity:03.1f}"

    def _get_Tgz_from_filename(self, filename):
        """
        A helper to get the temperature, logg, and metallicity
        from a PHOENIX spectrum model filename.

        Parameters
        ----------
        filename : string
            The filename of PHOENIX model.

        Returns
        -------
        temperature : float
            Temperature (K)
        logg : float
            log10([surface gravity]/[cm/s**2])
        metallicity : float
            log10([metallicity]/[solar metallicity])
        """

        f = os.path.basename(filename)
        temperature = float(f[3:8])
        logg = float(f[9:13])
        metallicity = float(f[13:17].replace("-0.0", "0.0"))
        return temperature, logg, metallicity

    def _get_filename_from_Tgz(self, temperature, logg, metallicity):
        """
        A helper to get a PHOENIX spectrum model filename
        from temperature, logg, and metallicity.

        Parameters
        ----------
        T : float
            Temperature (K)
        logg : float
            log10([surface gravity]/[cm/s**2])
        metallicity : float
            log10([metallicity]/[solar metallicity])

        Returns
        -------
        filename : string
            The filename of PHOENIX model (excluding directory).
        """

        f"lte{temperature:05.0f}-{logg:04.2f}{self._stringify_metallicity(metallicity)}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
        return f

    def _create_grid(self, R, metallicity=0.0, remake=False):
        """
        Create a pre-processed grid for a single resolution.

        Parameters
        ----------
        R : float
            The resolution of the grid.
        metallicity : float
            The stellar metallicity.
        remake : bool
            Should we remake the library even if a file exists?
        """

        # make sure that directory exists
        try:
            os.mkdir(self._directory_for_new_grids)
        except FileExistsError:
            pass

        try:
            assert self._current_raw_metallicity == metallicity
            assert len(self._raw_spectrum_filenames) == len(self._raw_downloaded)
        except (AttributeError, AssertionError):
            self._download_raw_data(metallicity=metallicity)

        # skip this resolution if already made
        filename = os.path.join(
            self._directory_for_new_grids,
            self._get_grid_filename(R, metallicity=metallicity),
        )
        if os.path.exists(filename) and (not remake):
            print(
                textwrap.dedent(
                    f"""
                a grid for R={R}, metallicity={metallicity} exists at
                {filename}
                so we're not remaking it
                """
                )
            )
            return

        print(
            f"Creating a new grid for R={R}, metallicity={metallicity}. Its details are..."
        )
        shared = {}
        shared["grid"] = "PHOENIX-ACES-AGSS-COND-2011"
        shared["url"] = "https://phoenix.astro.physik.uni-goettingen.de/?page_id=15"
        shared["citation"] = "2013A&A…553A…6H"
        shared["photons"] = self._are_the_units_photons
        shared["R"] = R
        shared["metallicity"] = metallicity
        shared["filename"] = os.path.basename(filename)
        unbinned_w = self._load_raw_wavelength()
        for k, v in shared.items():
            print(f"{k:>20} = {v}")

        d = {}
        for k, v in tqdm(list(self._raw_downloaded.items()), leave=False):

            # load the unbinned spectrum
            unbinned_f = self._load_raw_spectrum(v)

            # convert from W to photon/s
            if self._are_the_units_photons:
                photon_energy = (con.h * con.c / unbinned_w) / u.photon
                unbinned_f = (unbinned_f / photon_energy).to("ph/(s m**2 nm)")
            else:
                unbinned_f = unbinned_f.to("W/(m**2 nm)")
            # figure out its stellar inputs
            T, logg, Z = self._get_Tgz_from_filename(k)
            if R == "original":
                w = unbinned_w
                f = unbinned_f
            else:
                binned = bintoR(unbinned_w, unbinned_f, R=R, drop_nans=False)
                w, f = binned["x"], binned["y"]

            if "wavelength" not in shared:
                shared["wavelength"] = w.value

            assert len(f) == len(shared["wavelength"])

            key = (T, logg, Z)
            d[key] = f.value

        # pull out the unique values of the keys
        for i, k in enumerate(self._keys_for_indexing):
            shared[k] = np.unique([x[i] for x in d])

        shared["wavelength_unit"] = w.unit.to_string()
        shared["spectrum_unit"] = f.unit.to_string()

        # save everything to an easy-to-load file
        np.save(filename, [shared, d], allow_pickle=True)
        print(f"That grid has been saved to {filename}.\n")


phoenix_nextgen_library = Library_PHOENIX_NextGen(photons=True)


def get_phoenix_photons(
    temperature=5780,
    logg=4.43,
    metallicity=0.0,
    R=100,
    wavelength=None,
    wavelength_edges=None,
    visualize=False,
):
    """
    Get a PHOENIX model spectrum for an arbitrary temperature, logg, metallicity.

    Calculate the surface flux from a thermally emitted surface,
    according to PHOENIX model spectra, in units of photons/(s * m**2 * nm).

    Parameters
    ----------
    temperature : float, optional
        Temperature, in K (with no astropy units attached).
    logg : float, optional
        Surface gravity log10[g/(cm/s**2)] (with no astropy units attached).
    metallicity : float, optional
        Metallicity log10[metals/solar] (with no astropy units attached).
    R : float, optional
        Spectroscopic resolution (lambda/dlambda). Currently, this must
        be in one of [3,10,30,100,300,1000,3000,10000,30000,100000], but
        check back soon for custom wavelength grids. There is extra
        overhead associated with switching resolutions, so if you're
        going to retrieve many spectra, try to group by resolution.
        (If you're using the `wavelength` or `wavelength_edges` option
        below, please be ensure your requested R exceeds that needed
        to support your wavelengths.)
    wavelength : Quantity, optional
        A grid of wavelengths on which you would like your spectrum.
        If this is None, the complete wavelength array will be returned
        at your desired resolution. Otherwise, the spectrum will be
        returned exactly at those wavelengths. Grid points will be
        cached for this new wavelength grid to speed up applications
        that need to retreive lots of similar spectra for the same
        wavelength (like many optimization or sampling problems).
    wavelength_edges : Quantity, optional
        Same as `wavelength` (see above!) but defining the wavelength
        grid by its edges instead of its centers. The returned spectrum
        will have 1 fewer element than `wavelength_edges`.

    Returns
    -------
    wavelength : Quantity
        The wavelengths, at the specified resolution.
    photons : Quantity
        The surface flux in photon units
    """
    return phoenix_nextgen_library.get_spectrum(
        temperature=temperature,
        logg=logg,
        metallicity=metallicity,
        R=R,
        wavelength=wavelength,
        wavelength_edges=wavelength_edges,
        visualize=visualize,
    )
