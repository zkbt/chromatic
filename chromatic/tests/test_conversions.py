from ..rainbows import *

closekw = dict(rtol=1e-10)


def test_to_df():

    r = SimulatedRainbow(dt=1 * u.minute, R=50).inject_noise(signal_to_noise=10)

    r_df = r.to_df()

    # ensure the length of the df is the same as the original rainbow
    assert len(r_df) == r.nflux

    # ensure we have the right column names
    columnnames = r_df.columns
    for colname in ["Time (d)", "Wavelength (micron)", "Flux", "Flux Uncertainty"]:
        assert colname in columnnames

    # check the values in the df match the rainbow
    assert np.isclose(
        r_df["Time (d)"].values[0], r.time.to_value("h")[0] / 24.0, **closekw
    )  # default=days
    assert r_df["Wavelength (micron)"].values[0] == r.wavelength.to_value("micron")[0]
    assert r_df["Flux"].values[0] == r.flux[0, 0]
    assert r_df["Flux Uncertainty"].values[0] == r.uncertainty[0, 0]

    # test the timeformat parameter
    for t_unit in ["h", "hour", "day", "minute", "second", "s"]:
        r_df = r.to_df(t_unit=t_unit)
        assert f"Time ({t_unit})" in r_df.columns


def test_to_nparray():
    r = SimulatedRainbow(dt=1 * u.minute, R=50).inject_noise(signal_to_noise=100)

    rflux, rfluxu, rtime, rwavel = r.to_nparray()

    assert len(rtime) == r.ntime
    assert len(rwavel) == r.nwave
    assert len(rtime) * len(rwavel) == r.nflux
    assert len(rflux.flatten()) == r.nflux
    assert len(rfluxu.flatten()) == r.nflux
    assert np.shape(rflux) == r.shape
    assert np.shape(rfluxu) == r.shape

    assert np.all(rflux == r.flux)
    assert np.all(rfluxu == r.uncertainty)
    assert np.all(rwavel == r.wavelength.to_value("micron"))

    # issues with rounding errors:
    assert np.all(
        np.isclose(rtime, r.time.to_value("h") / 24.0, **closekw)
    )  # the default is days

    # test if the hours format works
    rflux, rfluxu, rtime, rwavel = r.to_nparray(t_unit="h")
    assert np.all(np.isclose(rtime, r.time.to_value("h"), **closekw))

    # test if the minutes format works
    rflux, rfluxu, rtime, rwavel = r.to_nparray(t_unit="min")
    assert np.all(np.isclose(rtime, r.time.to_value("h") * 60, **closekw))

    # test if the minutes format works
    rflux, rfluxu, rtime, rwavel = r.to_nparray(t_unit="s")
    assert np.all(np.isclose(rtime, r.time.to_value("h") * 3600, **closekw))

def test_to_spectra():
    times = np.arange(1, 30) * u.day
    r = SimulatedRainbow(dt=1 * u.minute, R=50, time=times).inject_noise(signal_to_noise=10)

    r_specs = r.to_spectra()

    # ensure we have a spectrum for each discrete time
    assert len(r_specs) == len(r.time)

    # ensure the length of each quanity in each spectrum is the same as the original rainbow WLs
    for index, obj in enumerate(r_specs):
        assert len(obj.spectral_axis) == len(r.wavelength)
        assert len(obj.flux) == len(r.flux[:, index])
        assert len(obj.uncertainty) == np.shape( r.uncertainty[:, index])[0]


    # check the values in the spectrums match the rainbow
    for index, obj in enumerate(r_specs):
        assert np.concat(obj.flux.value) == r.flux[:, index]
        assert np.concat(obj.uncertainty.value) == r.uncertainty[:, index]

    # test the wlformat parameter
    for w_unit, phys in zip(["micron", "nm", "Angstrom"], [u.micron, u.nm, u.Angstrom]):
        r_spec = r.to_spectra(w_unit=w_unit)
        assert r_spec[0].spectral_axis.unit == phys

    # test the time format parameter
    for t_unit, phys in zip(["h", "hour", "day", "minute", "second", "s"], [u.hour, u.hour, u.day, u.minute, u.second, u.second]):
        r_spec = r.to_spectra(t_unit=t_unit)
        assert r_spec[0].time.unit == phys

    # test the flux format parameter
    for f_unit, phys in zip(["W/m2/um", "erg/s/cm2/AA"], [u.Unit("W/m2/um"), u.Unit("erg/s/cm2/AA")]):
        r_spec = r.to_spectra(f_unit=f_unit)
        assert r_spec[0].flux.unit == phys