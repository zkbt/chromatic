from ...imports import *

__all__ = ["imshow"]


def imshow(
    self,
    ax=None,
    quantity="flux",
    w_unit="micron",
    t_unit="day",
    colorbar=True,
    aspect="auto",
    **kw,
):
    """
    imshow flux as a function of time (x = time, y = wavelength, color = flux).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes into which to make this plot.
    quantity : str
        The fluxlike quantity to imshow.
        (Must be a key of `rainbow.fluxlike`).
    w_unit : str, astropy.unit.Unit
        The unit for plotting wavelengths.
    t_unit : str, astropy.unit.Unit
        The unit for plotting times.
    colorbar : bool
        Should we include a colorbar?
    aspect : str
        What aspect ratio should be used for the imshow?
    """

    # self.speak(f'imshowing')
    if ax is None:
        ax = plt.subplot()

    w_unit, t_unit = u.Unit(w_unit), u.Unit(t_unit)

    if self.wscale == "linear":
        extent = [
            (min(self.time) / t_unit).decompose(),
            (max(self.time) / t_unit).decompose(),
            (max(self.wavelength) / w_unit).decompose(),
            (min(self.wavelength) / w_unit).decompose(),
        ]
        ylabel = f"Wavelength ({w_unit.to_string('latex_inline')})"

    elif self.wscale == "log":
        extent = [
            (min(self.time) / t_unit).decompose(),
            (max(self.time) / t_unit).decompose(),
            np.log10(max(self.wavelength) / w_unit),
            np.log10(min(self.wavelength) / w_unit),
        ]
        ylabel = r"log$_{10}$" + f"[Wavelength/({w_unit.to_string('latex_inline')})]"

    else:
        message = f"""
        The wavelength scale for this rainbow is '{self.wscale}'.
        It's hard to imshow something with a wavelength axis
        that isn't uniform in linear or logarithmic space, so
        we're giving up and just using the wavelength index
        as the wavelength axis.

        If you want a real wavelength axis, one solution would
        be to bin your wavelengths to a more uniform grid with
        `rainbow.bin(R=...)` (for logarithmic wavelengths) or
        `rainbow.bin(dw=...)` (for linear wavelengths)
        """
        warnings.warn(message)
        extent = [
            (min(self.time) / t_unit).decompose(),
            (max(self.time) / t_unit).decompose(),
            self.nwave,
            0,
        ]
        ylabel = "Wavelength Index"

    # define some default keywords
    imshow_kw = dict(interpolation="nearest")
    imshow_kw.update(**kw)
    with quantity_support():
        plt.sca(ax)
        plt.imshow(
            self.fluxlike[quantity],
            extent=extent,
            aspect=aspect,
            origin="upper",
            **imshow_kw,
        )
        plt.ylabel(ylabel)
        plt.xlabel(f"Time ({t_unit.to_string('latex_inline')})")
        if colorbar:
            plt.colorbar(
                ax=ax,
                label=u.Quantity(self.fluxlike[quantity]).unit.to_string(
                    "latex_inline"
                ),
            )
    return ax
