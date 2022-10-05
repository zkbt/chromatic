from ...imports import *


def attach_model(self, model, **kw):
    """
    Attach a fluxlike model, thus making a new `RainbowWithModel.`

    The `model` will be stored in `.fluxlike['model']`.
    After running this to make a `RainbowWithModel` it's OK
    (and faster) to simply update `.fluxlike['model']` or `.model`.

    Parameters
    ----------
    model : array, Quantity
        An array of model values, with the same shape as 'flux'
    **kw : dict
        All other keywords will be interpreted as items
        that can be added to a Rainbow. You might use this
        to attach intermediate model steps or quantities
        (for example, `'planet_model'` or `'systematics_model'`).

    Returns
    -------
    rainbow : RainbowWithModel
        A new `RainbowWithModel` object, with the model attached.
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
