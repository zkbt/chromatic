from ...imports import *


def attach_model(self, model, **kw):
    """
    Attach a fluxlike model to this Rainbow,
    storing it as `self.fluxlike['model']`

    After running this once, it's OK (and faster) to simply
    update the `.model` attribute of the result.

    Parameters
    ----------
    model : np.array, u.Quantity
        An array of model values, with the same shape as 'flux'
    kw : dict
        All other keywords will be interpreted as items
        that can be added to a Rainbow. You might use this
        to attach intermediate model steps or quantities
        (for example, 'planet_model' or 'systematics_model').
    """

    # create a history entry for this action (before other variables are defined)
    h = self._create_history_entry("attach_model", locals())

    # make sure the shape is reasonable
    assert np.shape(model) == np.shape(self.flux)

    # add the model to the fluxlike array
    inputs = self._create_copy()._get_core_dictionaries()
    inputs["fluxlike"]["model"] = model

    # import here (rather than globally) to avoid recursion?
    from ..withmodel import RainbowWithModel

    # create new object
    new = RainbowWithModel(**inputs)

    # add other inputs to the model
    for k, v in kw.items():
        new.__setattr__(k, v)

    # append the history entry to the new Rainbow
    new._record_history_entry(h)

    # return the RainboWithModel
    return new
