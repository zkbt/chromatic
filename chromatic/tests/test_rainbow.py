from ..rainbows import *


def test_basic_rainbow():

    nw, nt = 20, 40
    r = Rainbow(
        wavelength=np.linspace(0.5, 5, nw) * u.micron,
        time=np.linspace(-1, 1, nt) * u.hour,
        flux=np.ones((nw, nt)),
    )

    assert r.shape == (nw, nt)
    assert r.nwave == nw
    assert r.ntime == nt
    assert r.nflux == nw * nt

def test_bin_in_time():

    s = SimulatedRainbow(dt=5*u.minute, signal_to_noise=100)
    b = s.bin_in_time(dt=1*u.hour)
    assert(b.ntime < s.ntime)

    assert(s.bin_in_time(dt=None, time=None) == s)

    fi, ax = plt.subplots(2, 1, sharex=True)
    imshowkw = dict( vmin=0.98, vmax=1.02)
    s.imshow(ax=ax[0], **imshowkw)
    b.imshow(ax=ax[1], **imshowkw)
    plt.show()

def test_bin_in_wavelength():

    s = SimulatedRainbow(dw=50*u.nm)
    b = s.bin_in_wavelength(dw=500*u.nm)
    assert(b.nwave < s.nwave)

    assert(s.bin_in_wavelength(R=None, dw=None, wavelength=None) == s)

    fi, ax = plt.subplots(2, 1, sharex=True)
    imshowkw = dict( vmin=0.98, vmax=1.02)
    s.imshow(ax=ax[0], **imshowkw)
    b.imshow(ax=ax[1], **imshowkw)
    plt.show()

def test_bin():

    s = SimulatedRainbow(dt=10*u.minute, R=25)
    b = s.bin(dt=0.5*u.hour, dw=0.5*u.micron)
    assert(b.nwave < s.nwave)
    assert(b.ntime < s.ntime)

    assert(s.bin(dw=None, dt=None, R=None, time=None, wavelength=None) == s)

    fi, ax = plt.subplots(2, 1, sharex=True)
    imshowkw = dict( vmin=0.98, vmax=1.02)
    s.imshow(ax=ax[0], **imshowkw)
    plt.title('Unbinned')
    b.imshow(ax=ax[1], **imshowkw)
    plt.title('Binned')
    plt.tight_layout()
    plt.show()
