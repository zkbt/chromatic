from .setup_tests import *
from ..imports import *
from ..spectra import *


def test_spectral_library_R(cmap=one2another("indigo", "tomato"), N=5):
    fi, ax = plt.subplots(
        3, 1, figsize=(8, 6), dpi=300, sharex=True, sharey=True, constrained_layout=True
    )
    for R, a in zip([10, 100, 1000], ax):
        plt.sca(a)
        for i, T in enumerate(
            np.round(np.logspace(np.log10(2300), np.log10(12000), N))
        ):
            color = cmap(i / (N - 1))
            w, f = get_phoenix_photons(temperature=T, R=R)
            plt.loglog(w, f, color=color)
            w, f = get_planck_photons(temperature=T, R=R)
            plt.plot(w, f, color=color, linestyle=":")
            plt.ylim(1e21, 1e25)
            plt.xlim(0.05, 5)
            plt.text(0.98, 0.92, f"R={R}", transform=a.transAxes, ha="right", va="top")
    fi.supxlabel(f"Wavelength ({w.unit.to_string('latex_inline')})")
    fi.supylabel(f"Surface Flux ({f.unit.to_string('latex_inline')})")
    plt.savefig(
        os.path.join(
            test_directory, "demonstration-of-spectral-library-with-constant-R.pdf"
        )
    )


def test_spectral_library_wavelengths(cmap=one2another("indigo", "tomato"), N=5):
    fi, ax = plt.subplots(
        3, 1, figsize=(8, 6), dpi=300, sharex=True, sharey=True, constrained_layout=True
    )
    for how_many, a in zip([10, 100, 1000], ax):
        plt.sca(a)
        wavelength = np.logspace(-1, 0, how_many) * u.micron
        for i, T in enumerate(
            np.round(np.logspace(np.log10(2300), np.log10(12000), N))
        ):
            color = cmap(i / (N - 1))
            w, f = get_phoenix_photons(temperature=T, wavelength=wavelength)
            plt.loglog(w, f, color=color)
            w, f = get_planck_photons(temperature=T, wavelength=wavelength)
            plt.plot(w, f, color=color, linestyle=":")
            plt.ylim(1e21, 1e25)
            plt.xlim(0.05, 2)
            plt.text(
                0.98, 0.92, f"N={how_many}", transform=a.transAxes, ha="right", va="top"
            )
    fi.supxlabel(f"Wavelength ({w.unit.to_string('latex_inline')})")
    fi.supylabel(f"Surface Flux ({f.unit.to_string('latex_inline')})")
    plt.savefig(
        os.path.join(
            test_directory,
            "demonstration-of-spectral-library-with-custom-wavelengths.pdf",
        )
    )


def test_spectral_library_loads_correct_R():
    for i in range(2):
        w = np.linspace(0.7, 0.71, 1000) * u.micron
        wamb, s_amb = get_phoenix_photons(
            temperature=int(3900.0), wavelength=w, logg=4.4, metallicity=0.0
        )
        wspot, s_spot = get_phoenix_photons(
            temperature=int(3100.0), wavelength=w, logg=4.4, metallicity=0.0
        )
        plt.figure()
        plt.plot(wspot, s_spot)
        plt.plot(wamb, s_amb)
    plt.savefig(
        os.path.join(
            test_directory,
            "demonstration-of-spectral-library-loading-minimum-necessary-resolution.pdf",
        )
    )
