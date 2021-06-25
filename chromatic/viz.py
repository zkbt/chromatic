# Code to create plots of rainbows

import matplotlib.pyplot as plt
import numpy as np


def wavelength_plot(
    flux, time_q, wavelength_q, step_size=None, ax=None, plotkw={}, fontkw={}
):
    """
    Plot a time series for each wavelength bin

    Parameters
    ----------
    flux : astropy quantity
        time bins as an array of quantities
    time_q : quantity
        DESCRIPTION.
    wavelength_q : quantity
        wavelength bins as as array of quantities
    step_size : float, default=None
        Spacing between each time series.
        None uses half the standard dev of entire flux data.

    Returns
    -------
    None.

    """

    time = time_q.value
    time_unit = time_q.unit.to_string("latex_inline")
    wavelength = wavelength_q.value
    waveunit = wavelength_q.unit.to_string("latex_inline")

    fcadences = 0.05
    start_pcad = int(np.ceil(len(time) * fcadences))
    end_pcad = int(np.floor(len(time) * fcadences))

    cont_time = np.zeros(len(time))
    cont_time[0:start_pcad] = 1
    cont_time[end_pcad:-1] = 1

    min_time = np.min(time)

    if step_size is None:
        step_size = 3 * np.std(flux)

    nsteps = len(wavelength)

    if ax is None:
        fi, ax = plt.subplots(1, 1)
    # plt.figure(figsize = (8,12) )

    for i, w in enumerate(wavelength):

        lc = flux[i, :]
        # normalize this light curve to one.
        cont_level = np.nanmedian(lc[cont_time == 1])
        plot_flux = -i * step_size + lc / cont_level

        color = (0, 0.3 - 0.3 * (i / nsteps), i / nsteps)

        plt.plot(time, plot_flux, ".--", color=color, **plotkw)
        plt.annotate(
            f"{wavelength[i]:4.2f} {waveunit}",
            (min_time, 1 - i * step_size + 0.4 * step_size),
            color=color,
            **fontkw,
        )

    plt.xlabel(f"Time ({time_unit})")
    plt.ylabel("Relative Flux (+ offsets)")
