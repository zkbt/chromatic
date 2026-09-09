from ...imports import *
from .extract import *
from .visualize import *

# Details of PandExo files are [here](https://natashabatalha.github.io/PandExo/jwstdict.html?highlight=rawdata).


from exoatlas import *
from exoatlas.models import Mamajek
from pathlib import Path

def print_dict(d, prefix=' '):
    for k in d:
        if isinstance(d[k], dict):
            print(f'{prefix}{k} = ' + '{')
            print_dict(d[k], prefix=f'{prefix} ')
            print(f'{prefix}' + '}')
        else:
            try:
                assert len(d[k]) > 50
                f = np.array(d[k]).flatten()
                print(f'{prefix}{k} = [{f[0]} ... {f[-1]}]')
            except AttributeError, AssertionError, TypeError:
                print(f'{prefix}{k} = {d[k]}')

class Pandexo:

    # make an interpolator to get logg from Teff, defined at class level for faster runtime
    _mama = Mamajek()
    _logM_from_logT = _mama.tofrom('logM')('logT')
    _logR_from_logT = _mama.tofrom('logR')('logT')
    def _estimate_logg_from_Teff(self, Teff):
        logT = np.log10(Teff)
        M = 10**self._logM_from_logT(logT)*u.Msun
        R = 10**self._logR_from_logT(logT)*u.Rsun
        logg = np.log10((con.G*M/R**2).to_value('cm/s**2'))
        return logg


    def __init__(self, instrument_mode='NIRSpec Prism', **kw):

        try:
            import pandexo.engine.justdoit as jdi
            self.jdi = jdi
            import pandeia.engine
            pandeia.engine.pandeia_version()
            print(f'Using pandeia reference data from', os.environ['pandeia_refdata'])
        except ModuleNotFoundError, ImportError:
            raise ImportError(f'''
            Alas, you're trying to use Pandexo and Pandeia to calculate the predicted 
            noise in a JWST observation, but we're having trouble importing them. 
            We didn't include them as automatic dependencies in this package 
            because they are a smidge complicated to install, so we encourage 
            you to please install them via the installation instructions at 

            https://natashabatalha.github.io/PandExo/
            
            and then try again. Thanks for your patience! Good luck!
            ''')

        # load empty system dictionary
        self.exo_dict = self.jdi.load_exo_dict()

        # load default instrument dictionary
        available_instrument_modes = self.jdi.print_instruments(verbose=False)
        if instrument_mode in available_instrument_modes:
            self.inst_dict = self.jdi.load_mode_dict(instrument_mode)
            self.instrument_mode = instrument_mode
        else:
            raise ValueError(f'''
            Instrument mode "{instrument_mode}" is not available.
            Please choose from {available_instrument_modes}.
            ''')

        self.set_observation(**kw)
        self.set_star(**kw)
        self.set_planet(**kw)
        self.set_instrument(**kw)

    def set_observation(self, **kw):
        self.exo_dict['observation']['sat_level'] = 80 # (call >80% full well saturated)
        self.exo_dict['observation']['sat_unit'] = '%'
        self.exo_dict['observation']['noccultations'] = 1 # (number of transits)
        self.exo_dict['observation']['R'] = None # (no binning)
        self.exo_dict['observation']['baseline_unit'] = 'frac' # ('frac' : fraction of time in transit versus out = in/out)
        self.exo_dict['observation']['baseline'] = 1 # (equal time in/out of transit)
        self.exo_dict['observation']['noise_floor'] = 0  # (any excess noise floor to consider)

    def set_star(self, mag=10.0, Teff=5500, logg='auto', **kw):
        self.exo_dict['star']['type'] = 'phoenix' # (use PHOENIX stellar model)
        self.exo_dict['star']['mag'] = mag # (magnitude of the system)
        self.exo_dict['star']['ref_wave'] = 1.25 # (for J mag = 1.25, H = 1.6, K = 2.22, all in micron)
        self.exo_dict['star']['temp'] = Teff # (stellar Teff, in K)
        self.exo_dict['star']['metal'] = 0.0 # (stellar metallicity as log Fe/H)
        if logg == 'auto':
            self.exo_dict['star']['logg'] = self._estimate_logg_from_Teff(Teff) # (stellar log[g/(cm/s**2)])

    def directory(self):
        instrument = self.instrument_mode.replace(' ', '-')
        ngroups = self.inst_dict["configuration"]["detector"]["ngroup"]
        subarray = self.inst_dict['configuration']['detector']['subarray']
        d = Path(f'{instrument}-{ngroups}-{subarray}')
        d.mkdir(parents=True, exist_ok=True)
        return d

    def filename(self):
        magnitude = self.exo_dict['star']['mag']
        teff = self.exo_dict['star']['temp']
        instrument = self.instrument_mode.replace(' ', '-')
        ngroups = self.inst_dict["configuration"]["detector"]["ngroup"]
        subarray = self.inst_dict['configuration']['detector']['subarray']
        return f'{instrument}-{ngroups}-{subarray}-J={magnitude:.2f}-Teff={teff:.0f}K.p'

    def set_planet(self, duration=1*u.hour, **kw):
        self.exo_dict['planet']['type'] = 'constant' # (tells pandexo you want a fixed transit depth)
        self.exo_dict['planet']['transit_duration'] = duration.to_value('s') # (transit duration in seconds)
        self.exo_dict['planet']['td_unit'] = 's'
        self.exo_dict['planet']['f_unit'] = 'fp/f*' # (make eclipse with zero depth)
        self.exo_dict['planet']['temp'] = 0
        self.exo_dict['planet']['radius'] = 1
        self.exo_dict['planet']['r_unit'] = 'R_jup'
        self.exo_dict['star']['radius'] = 1
        self.exo_dict['star']['r_unit'] = 'R_sun'


    def set_instrument(self, ngroup='optimize', **kw):
        instrument_name = self.inst_dict['configuration']['instrument']['instrument']


        self.inst_dict["configuration"]["detector"]["ngroup"] = ngroup

        self.inst_dict['background'] = 'ecliptic'
        self.inst_dict['background_level'] = 'high'

        if 'subarray' in kw:
            available_subarrays = self.jdi.subarrays(instrument_name)
            if kw['subarray'] in available_subarrays:
                self.inst_dict['configuration']['detector']['subarray'] = kw['subarray']
            else:
                raise ValueError(f"{kw['subarray']} is not in {available_subarrays}")

        if 'filter' in kw:
            available_filters = self.jdi.filters(instrument_name)
            if kw['filter'] in available_filters:
                self.inst_dict['configuration']['instrument']['filter'] = kw['filter']
            else:
                raise ValueError(f"{kw['filter']} is not in {available_filters}")

        if 'disperser' in kw:
            available_dispersers = self.jdi.dispersers(instrument_name)
            if kw['disperser'] in available_dispersers:
                self.inst_dict['configuration']['instrument']['disperser'] = kw['disperser']
            else:
                raise ValueError(f"{kw['disperser']} is not in {available_dispersers}")

    def run(self, **kw):
        if len(kw) > 0:
            self.set_star(**kw)
            self.set_planet(**kw)
            self.set_observation(**kw)
            self.set_instrument(**kw)
        if os.path.exists(os.path.join(self.directory(), self.filename())):
            print(f'{self.filename()} already exists; skipping!')
            return

        print(f'Running pandexo for {self.filename()}')
        return self.jdi.run_pandexo(self.exo_dict, self.inst_dict,
                                    save_file=True,
                                    output_path=self.directory(),
                                    output_file=self.filename())


def match_one_array_to_another(incoming_x, outgoing_x):
    '''
    Find the indices to an incoming array that
    align (perfectly) with values in an
    outgoing array.

    Parameters
    ----------
    incoming_x : array_like
        An array of x-values (typically an independent variable
        like wavelength) from an array that we are trying to
        merge into an outgoing array.
    outgoing_x : array_like
        The array of x-values (nearly identical to the incoming
        ones, but possibly ordered differently) onto which
        we want to map the incoming values.

    Returns
    -------
    i_incoming : array_like
        The indices of the incoming array that successfully
        match to an element of the outgoing array.
    i_outgoing : array_like
        The indices of the outgoing array to which they match.
    '''
    i_outgoing = np.arange(len(outgoing_x))
    outgoing_has_match = np.zeros(len(outgoing_x)).astype(bool)
    i_incoming = np.zeros(len(outgoing_x)).astype(int)
    for i in i_outgoing:
        try:
            i_incoming[i] = np.flatnonzero(incoming_x == outgoing_x[i])[0]
            # (could change to isclose, if needed)
            outgoing_has_match[i] = True
            # print(i, i_incoming[i])
        except IndexError:
            outgoing_has_match[i] = False
            # print(f'{i} has no match!')

    i_incoming = i_incoming[outgoing_has_match]
    i_outgoing = i_outgoing[outgoing_has_match]

    return i_incoming, i_outgoing

def read_pandexo(filename : str | dict, extract : bool=False):
    """
    Read a PandExo output file, including S/N estimates.

    When calling PandExo, please use the following settings:
    - set "Baseline" to "Fraction of time: in/out" = 1
    - set "Number of Transits" = 1
    - set "Constant Minimum Noise" = 0 (unless you have a good reason to believe otherwise)

    Parameters
    ----------
    filename : str, dict
        The filepath to a `.p` pickle file generated by Pandexo,
        or simply the dictionary output of Pandexo's self.jdi.run_pandexo
    extract : bool
        Should we try to extract a S/N from a 2D image?
        (This currently likely works only for MIRI/LRS.
        It definitely won't work for NIRISS/SOSS.)

    Returns
    -------
    results : dict
        A dictionary containing noise estimates and intermediate ingredients.
        The typical keys are:
            `1D` = tabular results along the wavelength axis
            `2D` = image results along the wavelength axis and one spatial axis
            `3D` = cube results along the wavelength axis and two spatial axes
    """

    # read the pickle file
    if isinstance(filename, dict):
        d = filename
    else:
        f = open(filename, "rb")
        d = pickle.load(f)

    # get the overall pandexo depth and uncertainty results
    s = d["FinalSpectrum"]
    spectra = {}

    # all results are per-pixel
    spectra["pixel_number"] = np.arange(len(s["wave"]))
    spectra["wavelength"] = s["wave"]
    spectra["planet_model"] = s["spectrum"]
    spectra["planet_realization"] = s["spectrum_w_rand"]
    spectra["depth_uncertainty"] = s["error_w_floor"]
    spectra["snr_per_transit"] = 1 / spectra["depth_uncertainty"]
    spectra['full_saturation_mask'] = s['full_saturation_mask']

    # get some metadata about the timing.
    t = d["timing"]
    metadata = {}

    metadata["pandexo_input"] = d["input"]
    metadata["transit_duration"] = t["Transit Duration"]
    metadata["observation_duration"] = t["Transit+Baseline, no overhead (hrs)"]
    metadata["time_per_integration"] = t["Time/Integration incl reset (sec)"]
    metadata["observing_efficiency"] = 0.01 * t["Observing Efficiency (%)"]
    metadata["number_of_groups_per_integration"] = t["APT: Num Groups per Integration"]
    metadata["number_of_integrations_per_transit"] = t["Num Integrations In Transit"]
    metadata["number_of_groups_per_integration"] = t["APT: Num Groups per Integration"]
    metadata["time_per_group"] = t["Seconds per Frame"]

    # make sure that the number of transits is 1
    assert t["Number of Transits"] == 1

    # make sure the in/out transit ratio is about 1
    assert metadata["transit_duration"] / metadata["observation_duration"] > 0.4
    assert metadata["transit_duration"] / metadata["observation_duration"] < 0.6

    # collect some warnings as metadata
    w = d["PandeiaOutTrans"]["warnings"]
    for k in w:
        if "saturated" in k:
            metadata[k] = w[k]

    # collect the raw pandeia results, all per-pixel
    # Details [here](https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-outputs-overview/jwst-etc-downloads)
    pandeia_results = d["PandeiaOutTrans"]["1d"]
    for k in pandeia_results:
        # skip inputs for calculations that are not necessarily on the pixel grid
        if k in ["wave_calc", "target", "fp", "bg", "bg_rate", "total_flux"]:
            continue
        # most inputs have a wavelength axis embedded with them
        if len(pandeia_results[k]) == 2:
            this_wave, this_y = pandeia_results[k]
            i_this, i_spectra = match_one_array_to_another(this_wave, spectra['wavelength'])
            spectra[k] = np.ones(len(spectra['wavelength'])) * np.nan
            spectra[k][i_spectra] = this_y[i_this]

    # estimate from pure photon noise (no read noise or 1/f!)
    t_integration = metadata["time_per_integration"] * metadata["observing_efficiency"]
    n_integrations = metadata["number_of_integrations_per_transit"]
    convert_integration_to_transit = np.sqrt(n_integrations) / np.sqrt(2)

    N_photons = spectra["extracted_flux"] * t_integration
    sigma_N_photons = np.sqrt(spectra["extracted_flux_plus_bg"] * t_integration)
    spectra["snr_per_transit_from_photons_only"] = (
            N_photons / sigma_N_photons * convert_integration_to_transit
    )

    # get the per-integration noise from the direct ETC result
    spectra["snr_per_integration_from_etc"] = (
            spectra["extracted_flux"]
            / spectra["extracted_noise"]
            / convert_integration_to_transit
    )
    spectra["snr_per_integration_from_photons_only"] = N_photons / sigma_N_photons

    # make sure the number of integrations per transit lines up
    estimated_integrations_per_transit = (
            metadata["transit_duration"] * 60 * 60 / metadata["time_per_integration"]
    )
    assert np.isclose(
        estimated_integrations_per_transit,
        metadata["number_of_integrations_per_transit"],
        atol=2,
    )

    # make sure the number of groups makes sense
    N_groups = metadata["time_per_integration"] / metadata["time_per_group"] - 1
    assert np.isclose(N_groups, metadata["number_of_groups_per_integration"])
    # assert np.isclose(metadata["observing_efficiency"], (N_groups - 1) / (N_groups + 1))

    t = Table(spectra, meta=metadata)

    # get some images
    images = {}
    disperser = d["input"]["Disperser"]
    for k in ["detector", "snr", "saturation"]:
        images[k] = trim_image(d["PandeiaOutTrans"]["2d"][k], disperser=disperser)
    images["snr"] /= np.sqrt(metadata["number_of_integrations_per_transit"])

    if False:
        try:
            t["snr_extracted_from_image"] = extract_sn_from_image(images)
        except:
            print('for some reason, could not extract S/N from image')
    t["snr_per_integration_from_pandexo_depth"] = 1 / (
            t["depth_uncertainty"]
            * np.sqrt(metadata["number_of_integrations_per_transit"])
            / np.sqrt(2)
    )

    return {"1D": t, "2D": images}