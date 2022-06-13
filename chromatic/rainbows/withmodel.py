from .rainbow import *
from .visualizations import _add_panel_labels


class RainbowWithModel(Rainbow):
    """
    `RainbowWithModel` objects have a fluxlike `model`
    attached to them, meaning that they can

    This class definition inherits from `Rainbow`.
    """

    @property
    def residuals(self):
        """
        Calculate the residuals on the fly,
        to make sure they're always up to date.
        """
        return self.flux - self.model

    def _validate_core_dictionaries(self):
        super()._validate_core_dictionaries()
        try:
            model = self.get("model")
            assert np.shape(model) == np.shape(self.flux)
        except (AttributeError, AssertionError):
            message = """
            No fluxlike 'model' was found attached to this
            `RainbowWithModel` object. The poor thing,
            its name is a lie! Please connect a model.
            The simplest way to do so might look like...
            `rainbow.model = np.ones(rainbow.shape)`
            ...or similarly with a more interesting array.
            """
            warnings.warn(message)

    def imshow_data_with_models(
        self,
        models=["model"],
        vlimits_data=[0.98, 1.02],
        vspan_residuals=0.02,
        figsize=(8, 3),
        label="inside",
        labelkw={},
        **kw,
    ):
        """
        Produce a multi-panel imshow figure comparing
        the data to model components. The left panel
        will be the data ('flux'), the right panel will
        be the residuals ('flux' - 'model'), and the middle
        panel(s) will be the model(s) connecting them.

        Parameters
        ----------
        models : list
            The list of model quantities to display,
            each in its own panel. A popular option might be:
            `models = ['systematics_model', 'planet_model']
            Whatever you include here does not affect the
            values of the Residuals panel, which are always
            exactly 'flux' - 'model'.
        vlimits_data : list
            The value limits for the data and model colormap,
            specified as `vlimits_data = [vmin, vmax]`.
            If `None` or `[None, None]`, the limits will be
            replaced with the [1,99] percentiles of the flux.
        vspan_residuals : float
            The distance away from zero where the colormap
            for the residuals will cut off. If `None`, the
            span will be set to 3 standard deviations.
        figsize : tuple
            The figure size for the multipanel plot.
            It should probably be wider than it is tall.
        label : bool
            Should we add (a), (b), (c) labels to the panels,
            and where? False or None give no labels, 'inside'
            sets them in the upper left inside corner, 'outside'
            sets them above the upper left corner.
        labelkw : dict
            Keywords to pass along to the `plt.text` command
            that makes the labels.
        kw : dict
            All additional keywords will be passed along
            to `Rainbow.imshow` for all imshow panels.
        """
        # make sure color limits are set for data and model
        percentiles = [1, 99]
        if vlimits_data is None:
            vlimits_data = np.percentile(self.flux, percentiles)
        elif vlimits_data[0] is None:
            vlimits_data[0] = np.percentile(self.flux, percentiles[0])
        elif vlimits_data[1] is None:
            vlimits_data[1] = np.percentile(self.flux, percentiles[1])

        # make sure color limits are set for the residuals
        if vspan_residuals is None:
            vspan_residuals = np.nanstd(self.residuals) * 3

        # set up the panels
        row_data, row_colorbar = 1, 0
        N = len(models) + 2
        fi, ax = plt.subplots(
            2,
            N,
            figsize=figsize,
            dpi=600,
            gridspec_kw=dict(height_ratios=[1, 20]),
            sharex="row",
            sharey="row",
            constrained_layout=True,
        )

        # a helper function to make better titles
        def create_title(s):
            if s == "flux":
                return "Data"
            else:
                return s.replace("_", " ").title()

        # plot the first three panels
        for i, k in enumerate(["flux"] + models):
            self.imshow(
                ax=ax[row_data, i],
                quantity=k,
                vmin=vlimits_data[0],
                vmax=vlimits_data[1],
                colorbar=False,
                **kw,
            )
            ax[row_data, i].set_title(create_title(k))

        # set up the colorbar to span the data and models
        gs = ax[row_colorbar, 0].get_gridspec()
        for a in ax[row_colorbar, :-1]:
            a.remove()
        ax_colorbar = fi.add_subplot(gs[row_colorbar, :-1])
        plt.colorbar(cax=ax_colorbar, orientation="horizontal")

        # plot the residuals
        k = "residuals"
        self.imshow(
            ax[row_data, -1],
            k,
            vmin=-vspan_residuals,
            vmax=vspan_residuals,
            colorbar=False,
            **kw,
        )
        ax[row_data, -1].set_title("Residuals")

        # set up the colorbar for the residuals
        ax_colorbar_residuals = plt.subplot(gs[row_colorbar, -1])
        cbar = plt.colorbar(cax=ax_colorbar_residuals, orientation="horizontal")
        for a in [ax_colorbar, ax_colorbar_residuals]:
            a.tick_params(labelsize=8)
        x = vspan_residuals / 2
        cbar.set_ticks([-x, 0, x], labels=[f"{-x:.2}", "0", f"{x:.2}"])

        # clear all the axis labels and replace with shared
        for a in ax[row_data, :]:
            xlabel, ylabel = a.get_xlabel(), a.get_ylabel()
            a.set_ylabel("")
            a.set_xlabel("")
        fi.supxlabel(xlabel)
        fi.supylabel(ylabel)

        if label:
            _add_panel_labels(ax[row_data, :], preset=label, **labelkw)
