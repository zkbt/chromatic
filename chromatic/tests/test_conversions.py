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
