from ..imports import *
from ..resampling import bintoR, bintogrid
from ..version import __version__
import astropy.config.paths
from astropy.utils.data import is_url_in_cache, cache_total_size


class PHOENIXLibrary:
    # downloaded files will be stored in "~/.{_cache_label}"
    _cache_label = "chromatic"

    # these coordinates are needed to index the grid (in this order)
    _keys_for_indexing = ["temperature", "logg", "metallicity"]

    # the resolutions available as precalculated grids
    _available_resolutions = [
        3,
        10,
        30,
        100,
        300,
        1000,
        3000,
        10000,
        30000,
        100000,
    ]

    def get_cache_dir(self):
        return astropy.config.paths.get_cache_dir(self._cache_label)

    def get_cache_size(self):
        return (cache_total_size(self._cache_label) * u.byte).to(u.gigabyte)

    def __init__(self, directory=".", photons=True):
        """
        Initialize a PHOENIX model library to provide easy
        access to model stellar spectra at resolutions up
        to R=100000.

        Parameters
        ----------
        directory : str
            The path to a directory where *new* library grids
            will be created. This will only be used if you're
            calling the super secret `._create_grids` functions
            to download raw PHOENIX spectrum files and compile
            them into a new library. Most users can just
            not worry about this!
        photons : bool
            Should the units be in photons, rather than power?
            If True, spectrum units will be photons/(s * m**2 * nm)
            If False, spectrum units will be J/(s * m**2 * nm)
        """

        # should we use photon units (e.g. photons/s/m**2/nm)?
        # (if not, we'll use flux per wavelength units like W/m**2/nm)
        self._are_the_units_photons = photons

        # figure out the directory where the models will be stored
        self._directory_for_new_grids = os.path.join(
            os.path.join(directory),
            "chromatic-model-stellar-spectra",
        )

        self._local_paths = {}

    def _download_grid(self, R, metallicity=0.0, cache=True):
        """
        Download the preprocessed grid for a particular
        resolution, or load it from a cached local file.

        Parameters
        ----------
        R : float
            The resolution of the grid to retrieve
        metallicity : float
            The stellar metallicity (= log10[metals/solar])
        cache : bool
            Once it's downloaded, should we keep it for next time?

        Returns
        -------
        filename : str
            The path to the downloaded local file for the grid.
        """

        # make sure the resolution is reasonable
        if R not in self._available_resolutions:
            raise ValueError(
                f"""
            Your requested resolution of R={R} is not in
            {self._available_resolutions}
            Please pick one of those choices to download.
            """
            )

        # make sure it's a possible R
        assert R in self._available_resolutions

        # download (or find the cached local file)
        basename = self._get_grid_filename(R, metallicity=metallicity)
        url = f"https://casa.colorado.edu/~bertathompson/chromatic/{basename}"

        if is_url_in_cache(url, pkgname=self._cache_label):
            self._local_paths[R] = download_file(
                url, pkgname=self._cache_label, cache=True
            )
        else:

            threshold = 1000
            if R <= threshold:
                expected_time = f"Because the resolution is R<={threshold}, this should be pretty quick."
            else:
                expected_time = f"Because the resolution is R>{threshold}, this might be annoyingly slow."
            print(
                f"""
            Downloading pre-processed grid for R={R}, metallicity={metallicity} from
            {url}
            {expected_time}
            """
            )

            self._local_paths[R] = download_file_with_warning(
                url, pkgname=self._cache_label, cache=cache, show_progress=True
            )
        return self._local_paths[R]

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
        self._raw_local_paths[
            f"index-Z={metallicity_string}"
        ] = download_file_with_warning(
            self._raw_directory_url, pkgname=self._cache_label, cache=cache
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

    def _get_grid_filename(self, R, metallicity=0.0):
        """
        Get the filename of a pre-processed grid.

        Parameters
        ----------
        R : float
            The resolution of the grid.
        metallicity : float
            The stellar metallicity.

        Returns
        -------
        filename : str
            The filepath to the pre-processed grid.
        """
        if self._are_the_units_photons:
            unit_string = "photons"
        else:
            unit_string = "flux"

        filename = f"phoenix_{unit_string}_metallicity={metallicity:3.1f}_R={R:.0f}.npy"
        return filename

    def _create_grids(self, remake=False):
        """
        Create pre-processed grids for all resolutions.

        Parameters
        ----------
        metallicity : float
            The stellar metallicity.
        remake : bool
            Should we remake the library even if a file exists?
        """
        for metallicity in self._available_metallicities:
            for R in self._available_resolutions:
                self._create_grid(
                    R,
                    metallicity=metallicity,
                    remake=remake,
                )

    def _create_upload_command_for_grids(self, destination):
        """
        A tiny wrapper to generate a UNIX command that
        can upload locally-generated library files to
        an online source.

        Returns
        -------
        command : string
            A copy-paste `rsync` command to upload new files.
        """
        source = os.path.join(self._directory_for_new_grids, "phoenix_*.npy")
        command = f"rsync -v --progress {source} {destination}"
        return command

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
        shared["chromatic-version"] = __version__
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

    def _find_smallest_R(self, R):
        """
        Make sure the smallest grid necessary to support
        a particular resolution is loaded into the library.

        Parameters
        ----------
        R : float
            The resolution of the grid.
        """

        ok = np.array(self._available_resolutions) >= R
        if np.any(ok):
            i = np.flatnonzero(ok)[0]
        else:
            cheerfully_suggest(
                f"""
            Your request resolution of R={R} is higher than the maximum.
            It will be rounded down to R={np.max(self._available_resolutions)}.
            """
            )
            i = -1

        if i is not None:
            smallest_sufficient_R = self._available_resolutions[i]
        else:
            cheerfully_suggest(
                f"""
            Your requested resolution of R={R}
            is higher than the largest possible value
            with uniform binning (R={np.max(self._available_resolutions)}).
            """
            )
            smallest_sufficient_R = np.max(self._available_resolutions)
        return smallest_sufficient_R

    def _get_local_grid(self, R, metallicity=0.0, directory="."):
        bespoke = os.path.join(
            self._directory_for_new_grids,
            self._get_grid_filename(R, metallicity=metallicity),
        )
        if os.path.exists(bespoke):
            return bespoke
        else:
            return self._download_grid(R, metallicity=metallicity)

    def _load_grid(self, R, metallicity=0.0):
        """
        Load a pre-processed grid for a particular resolution.

        Parameters
        ----------
        R : float
            The resolution of the grid.
        metallicity : float
            The stellar metallicity.
        """
        try:
            assert self.metadata["R"] == R
            assert metallicity in self.metadata["metallicity"]
            return
        except (AttributeError, AssertionError):
            pass

        metadata, models = np.load(
            self._get_local_grid(R, metallicity=metallicity), allow_pickle=True
        )[()]

        try:
            for k in ["R", "photons", "wavelength"]:
                assert np.all(self.metadata[k] == metadata[k])
            for k in ["metallicity", "filename", "chromatic-version"]:
                self.metadata[k] = np.hstack([self.metadata[k], metadata[k]])
            self.models.update(**models)

        except (AttributeError, AssertionError):
            self.metadata = metadata
            self.models = models
            self.wavelength_cached_models = {}

        self.units = {
            k: u.Unit(self.metadata[f"{k}_unit"]) for k in ["wavelength", "spectrum"]
        }
        self.wavelength = self.metadata["wavelength"] * self.units["wavelength"]

    def _find_bounds(self, value, key, inputs={}):
        """
        Find a mask where `possible` is immediately
        below or above `value`.

        Parameters
        ----------
        value : float
            The exact values to find.
        possible : array
            The possible grid point values to compare to.

        Returns
        -------
         :
        """
        possible = self.metadata[key]

        # return one index where it's exact
        if value in possible:
            return [value]

        # figure out above and below grid points otherwise
        try:
            above = np.min(possible[possible > value])
        except ValueError:
            above = None
        try:
            below = np.max(possible[possible <= value])
        except ValueError:
            below = None

        if (above is None) or (below is None):
            raise ValueError(
                f"""
            Your requested coordinate of
            {repr(inputs)}
            is outside the limits {np.min(possible)}<={key}<{np.max(possible)}.
            Please rewrite your code to avoid this happening
            or proceed very, very, very, very cautiously.
            """
            )
        else:
            return below, above

    def _get_interpolation_weights(self, value, bounds):
        """
        Get the interpolation weights for `value`
        relative to the `bounds` below and above.
        """

        # if there's only one value, the weight should be 1
        if np.size(bounds) == 1:
            return [1]
        elif np.size(bounds) == 2:
            span = bounds[1] - bounds[0]
            weight_below = (bounds[1] - value) / span
            weight_above = (value - bounds[0]) / span
            return weight_below, weight_above

    _available_metallicities = [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0]

    def _wavelengths_to_hashable(self, w):
        """
        Convert an arbitrary array of wavelengths into
        something that can be used as a hash in a dictionary.

        Parameters
        ----------
        w : Quantity
            The wavelength array.

        Returns
        -------
        h : tuple
            A static tuple that summarizes the wavelength array
            with its minimum, maximum, and number of elements.
            (Technically this won't be perfectly unique for
            non-uniform grids, but it's probably good enough
            for most uses.)
        """
        return (np.min(w), np.max(w), len(w))

    def _wavelengths_to_R(self, w):
        """
        Calculate the minimum spectral resolution R that
        would be sufficient to support a set of wavelengths.

        Parameters
        ----------
        w : Quantity
            The wavelength array.

        Returns
        -------
        R : float
            The maxmium value of w/dw (= 1/dlnw) in the array.
        """
        return 1 / np.min(np.diff(np.log(w.value)))

        return (np.min(w), np.max(w), len(w))

    def _get_spectrum_from_grid(self, key, wavelength=None, wavelength_edges=None):
        if (wavelength is None) and (wavelength_edges is None):
            return self.models[key]
        else:
            # make sure ask only for one type of wavelength
            assert (wavelength is None) or (wavelength_edges is None)

            # figure out the key associated with these wavelengths
            if wavelength is None:
                wavelength_key = self._wavelengths_to_hashable(wavelength_edges)
            elif wavelength_edges is None:
                wavelength_key = self._wavelengths_to_hashable(wavelength)

            # make sure the dictionary exists
            try:
                self.wavelength_cached_models[wavelength_key]
            except KeyError:
                self.wavelength_cached_models[wavelength_key] = {}

            # get or populate the necessary entry in the grid
            try:
                return self.wavelength_cached_models[wavelength_key][key]
            except KeyError:
                self.wavelength_cached_models[wavelength_key][key] = bintogrid(
                    self.wavelength,
                    self.models[key],
                    newx=wavelength,
                    newx_edges=wavelength_edges,
                )["y"]
                return self.wavelength_cached_models[wavelength_key][key]

    def get_spectrum(
        self,
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

        # this kludgy business is to try to avoid loading multiple metallicities unless absolutely necessary
        try:
            # if no wavelength grids are provided, just go with R
            if (wavelength is None) and (wavelength_edges is None):
                assert self.metadata["R"] == R
            # if wavelength grids are provided, make sure the grid has high enough R
            else:
                # make sure ask only for one type of wavelength
                assert (wavelength is None) or (wavelength_edges is None)

                # figure out the minimum R that'd be enough
                if wavelength is None:
                    necessary_R = self._wavelengths_to_R(wavelength_edges)
                elif wavelength_edges is None:
                    necessary_R = self._wavelengths_to_R(wavelength)

                try:
                    assert self.metadata.get("R", 0) >= necessary_R
                except (AttributeError, AssertionError):
                    R = self._find_smallest_R(necessary_R)
                    assert False

            assert (metallicity >= np.min(self.metadata["metallicity"])) and (
                metallicity <= np.max(self.metadata["metallicity"])
            )
        except (AttributeError, AssertionError):
            if metallicity in self._available_metallicities:
                self._load_grid(
                    self._find_smallest_R(R=R),
                    metallicity=metallicity,
                )
            else:
                for m in self._available_metallicities:
                    self._load_grid(
                        self._find_smallest_R(R=R),
                        metallicity=m,
                    )

        # strip units
        if isinstance(temperature, u.Quantity):
            temperature = temperature.value

        # store the inputs as a convenient dictionary to pass around
        inputs = dict(temperature=temperature, logg=logg, metallicity=metallicity)

        # figure out the
        bounding_temperature = self._find_bounds(
            temperature, key="temperature", inputs=inputs
        )
        bounding_logg = self._find_bounds(logg, key="logg", inputs=inputs)
        bounding_metallicity = self._find_bounds(
            metallicity, key="metallicity", inputs=inputs
        )

        if visualize:
            self.plot_available(**inputs)

        N = (
            np.size(bounding_temperature)
            * np.size(bounding_logg)
            * np.size(bounding_metallicity)
        )
        spectra = []
        if N == 1:
            weights = [1]
            key = (bounding_temperature[0], bounding_logg[0], bounding_metallicity[0])
            spectrum = self._get_spectrum_from_grid(
                key, wavelength=wavelength, wavelength_edges=wavelength_edges
            ).flatten()
        else:
            logT, bounding_logT = np.log(temperature), np.log(bounding_temperature)
            weight_temperature = self._get_interpolation_weights(logT, bounding_logT)
            weight_logg = self._get_interpolation_weights(logg, bounding_logg)
            weight_metallicity = self._get_interpolation_weights(
                metallicity, bounding_metallicity
            )

            weights = np.zeros(N)
            i = 0
            for wt, t in zip(weight_temperature, bounding_temperature):
                for wg, g in zip(weight_logg, bounding_logg):
                    for wz, z in zip(weight_metallicity, bounding_metallicity):
                        weights[i] = wt * wg * wz
                        key = (t, g, z)
                        try:
                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore")
                                this_log_spectrum = np.log(
                                    self._get_spectrum_from_grid(
                                        key,
                                        wavelength=wavelength,
                                        wavelength_edges=wavelength_edges,
                                    )
                                )
                        except KeyError:
                            raise ValueError(
                                f"""
                            The grid point coordinate {key} is needed to interpolate to your
                            requested coordinate {inputs},
                            but it's not available. This problem is probably caused by
                            requesting something near the jagged edge of the
                            grid's availability.

                            Please run `.plot_available()` to see what models are possible.
                            """
                            )
                        spectra.append(this_log_spectrum)
                        # print(f"{weights[i]} * [{key}]")
                        i += 1
            weight_sum = np.sum(weights)
            assert np.isclose(weight_sum, 1)
            spectrum = np.exp(np.sum(weights[:, np.newaxis] * spectra, axis=0))

        if wavelength is None:
            if wavelength_edges is None:
                wavelength = self.wavelength
            else:
                wavelength = 0.5 * (wavelength_edges[:-1] + wavelength_edges[1:])

        if visualize:
            fi = plt.figure(figsize=(8, 3))
            for w, s in zip(weights, spectra):
                plt.plot(wavelength, np.exp(s), alpha=w)
            plt.plot(wavelength, spectrum, color="black")
            plt.xlabel(
                f"Wavelength ({self.units['wavelength'].to_string('latex_inline')})"
            )
            plt.ylabel(
                f"Model Flux\n({self.units['spectrum'].to_string('latex_inline')})"
            )
            plt.suptitle(", ".join([f"{k}={v}" for k, v in inputs.items()]))

        return wavelength, spectrum * self.units["spectrum"]

    def plot_available(self, temperature=None, logg=None, metallicity=None):
        """
        Make plots indicating the available grid points
        in the loaded library of stellar spectra.

        Parameters
        ----------
        temperature : float
            Temperature, in K (with no astropy units attached).
        logg : float
            Surface gravity log10[g/(cm/s**2)] (with no astropy units attached).
        metallicity : float
            Metallicity log10[metals/solar] (with no astropy units attached).
        """
        N = len(self._keys_for_indexing)
        fi, ax = plt.subplots(
            N, N, figsize=(8, 8), sharex="col", sharey="row", constrained_layout=True
        )
        labels = dict(
            temperature="Temperature (K)",
            logg="$log_{10}[g/(cm/s^2)]$",
            metallicity="[Fe/H] (metallicity)",
        )
        gridkw = dict(marker=".", alpha=0.3)
        for i, ky in enumerate(self._keys_for_indexing):
            for j, kx in enumerate(self._keys_for_indexing):
                y = [k[i] for k in self.models]
                x = [k[j] for k in self.models]
                plt.sca(ax[i, j])
                plt.scatter(x, y, **gridkw)
                plt.scatter(locals()[kx], locals()[ky])
                # if j == 0:
                plt.ylabel(labels[ky])
                # if i == (N - 1):
                plt.xlabel(labels[kx])
        plt.suptitle(f"R={self.metadata['R']}")

    def plot_time_required(self, iterations=5):
        """
        Run benchmarking tests to determine how long it takes
        to perform different actions with the spectral library.

        Parameters
        ----------
        iterations : int
            How many times should we run each action?
            (definitely do something more than 1, because
            there are some overheads that behave differently
            between the first time and subsequent times you
            you generate a model spectrum)
        """
        timings = {}
        for k in [
            "R",
            "load library",
            "get spectrum at grid point",
            "get interpolated spectrum",
            "get interpolated spectrum at specific wavelengths",
        ]:
            timings[k] = []

        for R in tqdm(self._available_resolutions, leave=False):
            for i in range(iterations):
                timings["R"].append(R)

                # how long does it take to load?
                start = get_current_seconds()
                self._load_grid(R)
                dt = get_current_seconds() - start
                timings["load library"].append(dt)

                # how long does it take to retrieve a spectrum that exists?
                start = get_current_seconds()
                self.get_spectrum(temperature=3000, logg=4.5, metallicity=0.0, R=R)
                dt = get_current_seconds() - start
                timings["get spectrum at grid point"].append(dt)

                # how long does it take to retrieve a spectrum that needs to be interpolated?
                start = get_current_seconds()
                self.get_spectrum(temperature=3456, logg=5.67, metallicity=0.0, R=R)
                dt = get_current_seconds() - start
                timings["get interpolated spectrum"].append(dt)

                # how long does it take to retrieve a spectrum on a particular wavelength grid
                wavelength = np.exp(np.arange(-1, 1, 1.01 / R)) * u.micron
                start = get_current_seconds()
                self.get_spectrum(
                    temperature=3456, logg=5.67, metallicity=0.0, wavelength=wavelength
                )
                dt = get_current_seconds() - start
                timings["get interpolated spectrum at specific wavelengths"].append(dt)

        t = Table(timings)
        plt.figure(figsize=(8, 4), dpi=300)
        for k in timings.keys():
            if k != "R":
                plt.loglog(timings["R"], timings[k], marker="o", label=k, alpha=0.3)
        plt.legend(bbox_to_anchor=(1, 1), frameon=False)
        plt.xlabel(r"R = $\lambda/\Delta\lambda$")
        plt.ylabel("Time Required")
        return t


phoenix_library = PHOENIXLibrary(photons=True)


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
    return phoenix_library.get_spectrum(
        temperature=temperature,
        logg=logg,
        metallicity=metallicity,
        R=R,
        wavelength=wavelength,
        wavelength_edges=wavelength_edges,
        visualize=visualize,
    )
