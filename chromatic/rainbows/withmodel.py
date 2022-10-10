from .rainbow import *


class RainbowWithModel(Rainbow):
    """
    `RainbowWithModel` objects have a fluxlike `model`
    attached to them, meaning that they can

    This class definition inherits from `Rainbow`.
    """

    # which fluxlike keys will respond to math between objects
    _keys_that_respond_to_math = ["flux", "model"]

    # which keys get uncertainty weighting during binning
    _keys_that_get_uncertainty_weighting = ["flux", "model", "uncertainty"]

    @property
    def residuals(self):
        """
        Calculate the residuals on the fly,
        to make sure they're always up to date.

        The residuals are calculated simply
        as the `.flux` - `.model`, so they are
        in whatever units those arrays have.

        Returns
        -------
        residuals : array, Quantity
            The 2D array of residuals (nwave, ntime).
        """
        return self.flux - self.model

    @property
    def chi_squared(self):
        """
        Calculate $\chi^2$.

        This calculates the sum of the squares of
        the uncertainty-normalized residuals,
        sum(((flux - model)/uncertainty)**2)

        Data points marked as not OK are ignored.

        Returns
        -------
        chi_squared : float
            The chi-squared value.
        """
        r = (self.flux - self.model) / self.uncertainty
        return np.sum(r[self.ok] ** 2)

    @property
    def residuals_plus_one(self):
        """
        A tiny wrapper to get the residuals + 1.

        Returns
        -------
        residuals_plus_one : array, Quantity
            The 2D array of residuals + 1 (nwave, ntime).
        """
        return self.flux - self.model + 1

    @property
    def ones(self):
        """
        Generate an array of ones that looks like the flux.
        (A tiny wrapper needed for `plot_with_model`)

        Returns
        -------
        ones : array, Quantity
            The 2D array ones (nwave, ntime).
        """
        return np.ones_like(self.flux)

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
            cheerfully_suggest(message)

    from .visualizations import (
        plot_with_model,
        plot_with_model_and_residuals,
        imshow_with_models,
        plot_one_wavelength_with_models,
        animate_with_models,
    )


# REMOVE THE RAINBOW WITH MODEL AND JUST ADD A VALIDATION STEP TO ALL MODEL-DEPENDENT THINGS?
