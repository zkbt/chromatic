from ...imports import *
import json
from .extract import *
from .visualize import *


def read_etc(directory, t_integration=None, extract=False):
    """
    Read the outputs of the JWST ETC into a table and some images.

    Parameters
    ----------
    directory : str
        The filepath to the directory downloaded from the ETC.
    t_integration : float
        The integration time in seconds (used to provide a photon-only S/N).
    extract : bool
        Should we try to extract a S/N from a 2D image?

    Returns
    -------
    table : array
        1D extracted results
    images : array
        2D image results

    """

    # find all the FITS files in the directory
    fits_filenames = glob.glob(os.path.join(directory, "*/*.fits"))
    input_filename = os.path.join(directory, "input.json")

    # load the metadata inputs
    with open(input_filename) as f:
        metadata = json.load(f)
    metadata.pop("background")

    spectra = {}
    images = {}
    image_slice = set_image_slice(metadata["configuration"]["instrument"]["disperser"])

    for f in fits_filenames:

        if "lineplot" in f:
            k = os.path.basename(f).split(".fits")[0].replace("lineplot_", "")
            if k in ["wave_calc", "target", "fp", "bg", "bg_rate", "total_flux"]:
                continue
            print(k)
            if k == "wave_pix" in k:
                spectra["wavelength"] = fits.open(f)[1].data["wavelength"]
            else:
                spectra[k] = fits.open(f)[1].data[k]
        if "image" in f:
            k = os.path.basename(f).split(".fits")[0].replace("image_", "")
            images[k] = fits.open(f)[0].data[:, image_slice]

    n_integrations_per_exposure = metadata["configuration"]["detector"]["nint"]
    spectra["snr_per_exposure"] = spectra["extracted_flux"] / spectra["extracted_noise"]
    spectra["snr_per_integration"] = spectra["snr_per_exposure"] / np.sqrt(
        n_integrations_per_exposure
    )

    if t_integration is not None:
        signal = spectra["extracted_flux"] * t_int
        noise = np.sqrt(
            (spectra["extracted_flux"] + spectra["extracted_bg_total"]) * t_int
        )
        spectra["snr_per_exposure_from_photons_only"] = signal / noise
        spectra["snr_per_integration_from_photons_only"] = (
            signal / noise / np.sqrt(n_integrations_per_exposure)
        )
        metadata["time_per_integration"] = t_integration

    t = Table(spectra, meta=metadata)
    if extract:
        t["snr_extracted_from_image"] = extract_sn_from_image(images)

    return t, images
