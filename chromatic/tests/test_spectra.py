from ..spectra import *


def test_planck_flux():
    w = np.logspace(-2, 4, 100000) * u.micron
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for T in [5, 50, 500, 5000, 50000] * u.K:
            bolometric_analytic = con.sigma_sb * T**4
            f = calculate_planck_flux(w, T)
            bolometric_numerical = np.trapz(f, w).decompose()
            print(T, bolometric_analytic, bolometric_numerical)
            assert np.isclose(bolometric_numerical, bolometric_analytic, rtol=1e-2)
