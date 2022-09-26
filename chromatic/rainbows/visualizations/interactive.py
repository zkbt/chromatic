from ...imports import *

# don't make Altair a necessary part of chromatic
try:
    import altair as alt

    alt.data_transformers.disable_max_rows()
except Exception as e:
    print(e)
    cheerfully_suggest(
        "Issue importing Altair, cannot make interactive plot :(! \n \
                  You can install Altair using: pip install altair"
    )

__all__ = ["imshow_interact"]


# Convert this grid to columnar data expected by Altair
def imshow_interact(
    self,
    quantity="Flux",
    t_unit="d",
    w_unit="micron",
    cmap="viridis",
    ylim=[],
    ylog=None,
    filename=None,
):
    """
    Display interactive spectrum plot for chromatic Rainbow with a
    wavelength-averaged 2D quantity defined by the user. The user
    can interact with the 3D spectrum to choose the wavelength range
    over which the average is calculated.

    Parameters
    ----------
    self : Rainbow object
        chromatic Rainbow object to plot
    quantity : str
        (optional, default='flux')
        The quantity to imshow, currently either `flux` or `uncertainty`
    ylog : boolean
        (optional, default=None)
        Boolean for whether to take log10 of the y-axis data.
        If None, will be guessed from the data.
    t_unit : str
        (optional, default='d')
        The time unit to use (seconds, minutes, hours, days etc.)
    w_unit : str
        (optional, default='micron')
        The wavelength unit to use
    cmap : str
        (optional, default='viridis')
        The color scheme to use from Vega documentation
    ylim : list
        (optional, default=[])
        If the user wants to define their own ylimits on the lightcurve plot
    """

    # preset the x and y axes as Time (in units defined by the user) and Wavelength
    xlabel = f"Time ({t_unit})"
    ylabel = f"Wavelength ({w_unit})"

    # allow the user to plot flux or uncertainty
    if quantity.lower() == "flux":
        z = "Flux"
    elif quantity.lower() == "uncertainty":
        z = "Flux Uncertainty"
    elif quantity.lower() == "error":
        z = "Flux Uncertainty"
    elif quantity.lower() == "flux_error":
        z = "Flux Uncertainty"
    elif quantity.lower() == "flux_uncertainty":
        z = "Flux Uncertainty"
    else:
        # if the quantity is not one of the predefined values:
        cheerfully_suggest("Unrecognised Quantity!")
        return

    # convert rainbow object to pandas dataframe
    source = self.to_df(t_unit=t_unit, w_unit=w_unit)[[xlabel, ylabel, z]]

    # if there are >10,000 data points Altair will be very laggy/slow. This is probably unbinned, therefore
    # encourage the user to bin the Rainbow before calling this function in future/
    N_warning = 100000
    if len(source) > N_warning:
        cheerfully_suggest(
            f"""
        The dataset {self} has {N_warning} data points/
        The interactive plot may lag. Try binning first!
        """
        )

    if (self._is_probably_normalized() == False) and "model" not in self.fluxlike:
        cheerfully_suggest(
            """
        It looks like you might be trying to use `imshow_interact` with an
        unnormalized Rainbow object. You might consider normalizing first
        with `rainbow.normalize().imshow_interact()`.
        """
        )

    # The unbinned Rainbow is sometimes in log scale, therefore plotting will be ugly with uniform axis spacing
    # ylog tells the function to take the log10 of the y-axis data
    try:
        ylog = ylog or (self.wscale == "log")
    except AttributeError:
        ylog = ylog or False

    if ylog:
        source[ylabel] = np.log10(source[ylabel])
        source = source.rename(columns={ylabel: f"log10({ylabel})"})
        ylabel = f"log10({ylabel})"

    if len(ylim) > 0:
        domain = ylim
    else:
        domain = [source[z].min() - 0.005, source[z].max() + 0.005]

    # Add interactive part
    brush = alt.selection(type="interval", encodings=["y"])

    # Define the 3D spectrum plot
    spectrum = (
        alt.Chart(source, width=280, height=230)
        .mark_rect(
            clip=False,
            width=280 / len(self.timelike["time"]),
            height=230 / len(self.wavelike["wavelength"]),
        )
        .encode(
            x=alt.X(
                f"{xlabel}:Q",
                scale=alt.Scale(
                    zero=False,
                    nice=False,
                    domain=[np.min(source[xlabel]), np.max(source[xlabel])],
                ),
            ),
            y=alt.Y(
                f"{ylabel}:Q",
                scale=alt.Scale(
                    zero=False,
                    nice=False,
                    domain=[np.max(source[ylabel]), np.min(source[ylabel])],
                ),
            ),
            fill=alt.Color(
                f"{z}:Q",
                scale=alt.Scale(
                    scheme=cmap,
                    zero=False,
                    domain=domain,
                ),
            ),
            tooltip=[f"{xlabel}", f"{ylabel}", f"{z}"],
        )
    )

    # gray out the background with selection
    background = spectrum.encode(color=alt.value("#ddd")).add_selection(brush)

    # highlights on the transformed data
    highlight = spectrum.transform_filter(brush)

    # Layer the various plotting parts
    spectrum_int = alt.layer(background, highlight, data=source)

    # Add the 2D averaged lightcurve (or uncertainty)
    lightcurve = (
        alt.Chart(source, width=280, height=230, title=f"Mean {z} for Wavelength Range")
        .mark_point(filled=True, size=20, color="black")
        .encode(
            x=alt.X(
                f"{xlabel}:Q",
                scale=alt.Scale(
                    zero=False,
                    nice=False,
                    domain=[
                        np.min(source[xlabel])
                        - (0.02 * np.abs(np.min(source[xlabel]))),
                        1.02 * np.max(source[xlabel]),
                    ],
                ),
            ),
            y=alt.Y(
                f"mean({z}):Q",
                scale=alt.Scale(zero=False, domain=domain),
                title="Mean " + z,
            ),
        )
        .transform_filter(brush)
    )

    # display the interactive Altair plot
    (spectrum_int | lightcurve).display()
    if filename is not None:
        (spectrum_int | lightcurve).save(self._label_plot_file(filename))
