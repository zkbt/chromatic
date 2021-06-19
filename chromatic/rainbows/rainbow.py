from ..imports import *
from ..talker import Talker


class Rainbow:
    """
    Rainbow objects represent the flux of an object
    as a function of both wavelength and time.
    """

    def __init__(self, wavelength=None, time=None, flux=None, uncertainty=None, **kw):
        """
        Initialize a generic Rainbow object.

        Parameters
        ----------
        wave : astropy.unit.Quantity
            A 1D array of wavelengths, in any unit.
        time : astropy.unit.Quantity or astropy.unit.Time
            A 1D array of times, in any unit.
        flux : np.array
            A 2D array of flux values.
        uncertainty : np.array
            A 2D array of uncertainties, associated with the flux.

        """

        # wavelike quanities are 1D arrays with nwave elements
        self.wavelike = {}
        self.wavelike["wavelength"] = wavelength

        # timelike quantities are 1D arrays with ntime elements
        self.timelike = {}
        self.timelike["time"] = time

        # fluxlike quantities are 2D arrays with nwave x time elements
        self.fluxlike = {}
        self.fluxlike["flux"] = flux
        self.fluxlike["uncertainty"] = uncertainty

    @property
    def wavelength(self):
        return self.wavelike["wavelength"]

    # @wavelength.setter
    # def wavelength(self, value):
    #    self.wavelike['wavelength'] = value

    @property
    def time(self):
        return self.timelike["time"]

    @property
    def flux(self):
        return self.fluxlike["flux"]

    @property
    def uncertainty(self):
        return self.fluxlike["uncertainty"]

    @property
    def shape(self):
        return (self.nwave, self.ntime)

    @property
    def nwave(self):
        if self.wavelength is None:
            return 0
        else:
            return len(self.wavelength)

    @property
    def ntime(self):
        if self.time is None:
            return 0
        else:
            return len(self.time)

    @property
    def nflux(self):
        return np.prod(self.shape)

    def __repr__(self):
        n = self.__class__.__name__
        return f"<{n} ({self.nwave}w, {self.ntime}t)>"

    def imshow(
        self,
        ax=None,
        w_unit="micron",
        t_unit="hour",
        aspect="auto",
        colorbar=True,
        origin="upper",
        **kw,
    ):

        if ax is None:
            ax = plt.gca()

        w_unit, t_unit = u.Unit(w_unit), u.Unit(t_unit)

        if self.wscale == "linear":
            extent = [
                (min(self.time) / t_unit).decompose(),
                (max(self.time) / t_unit).decompose(),
                (min(self.wavelength) / w_unit).decompose(),
                (max(self.wavelength) / w_unit).decompose(),
            ]
        elif self.wscale == "log":
            extent = [
                (min(self.time) / t_unit).decompose(),
                (max(self.time) / t_unit).decompose(),
                np.log10(min(self.wavelength) / w_unit),
                np.log10(max(self.wavelength) / w_unit),
            ]
        else:
            raise RuntimeError("Can't imshow without knowing wscale.")

        with quantity_support():
            plt.sca(ax)
            plt.imshow(self.flux, extent=extent, aspect=aspect, origin=origin)
            if self.wscale == "linear":
                plt.ylabel(f"Wavelength ({w_unit.to_string('latex_inline')})")
            elif self.wscale == "log":
                plt.ylabel(
                    r"log$_{10}$" + f"[Wavelength/({w_unit.to_string('latex_inline')})]"
                )
            plt.xlabel(f"Time ({t_unit.to_string('latex_inline')})")
            if colorbar:
                plt.colorbar(ax=ax)
        return ax
