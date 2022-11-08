"""
Add ability for Rainbow to keep track of their own history.
"""
from ...imports import *

__all__ = [
    "_setup_history",
    "_record_history_entry",
    "_remove_last_history_entry",
    "_create_history_entry",
    "history",
]


def _setup_history(self):
    """
    Record the first history entry in the log.

    Parameters
    ----------
    h : dict
        A dictionary containing keys `name` as the name
        of the action and `inputs` as its keyword inputs.
    """
    self.metadata["history"] = []


def _record_history_entry(self, h):
    """
    Record a history entry in the log.

    Parameters
    ----------
    h : dict
        A dictionary containing keys `name` as the name
        of the action and `inputs` as its keyword inputs.
    """
    self.metadata["history"].append(h)


def _remove_last_history_entry(self):
    """
    Remove the most recent history entry from the log
    """
    try:
        self.metadata["history"].pop()
    except KeyError:
        self.metadata["history"] = []


def represent_as_copypasteable(x):
    """
    Represent a particular quantity as a copy-pastable string
    within the chromatic import environment.

    Parameters
    ----------
    x : any
        The thing to represent

    Returns
    -------
    r : str
        The thing represented as a string
    """
    if isinstance(x, u.Quantity):
        value = represent_as_copypasteable(x.value)
        unit = x.unit
        return f"u.Quantity({value})*u.Unit('{unit}')"
    elif isinstance(x, np.ndarray):
        return f"np.{repr(x)}"
    elif isinstance(x, slice):
        return (
            f"{x.start}:{x.stop}:{x.step}".replace("None", "")
            .replace("::", ":")
            .replace("::", ":")
        )
    else:
        return repr(x)


def _create_history_entry(self, name, inputs={}):
    """
    Create a history entry.

    Parameters
    ----------
    name : str
        The name of the action being applied.
    inputs : dict
        Dictionary of all inputs.

    Returns
    -------
    h : dict
        A dictionary containing keys `name` as the name
        of the action and `inputs` as its keyword inputs.
    """

    # remove "self" from the list of inputs
    inputs.pop("self", None)

    # remove None inputs
    to_remove = [k for k in inputs if inputs[k] is None]
    for k in to_remove:
        inputs.pop(k, None)

    # create a dictionary to store function call components
    if name == "__getitem__":
        w = represent_as_copypasteable(inputs["i_wavelength"])
        t = represent_as_copypasteable(inputs["i_time"])
        call = f"[{w},{t}]"
    elif name in "+-*/":
        try:
            other = f'\n{textwrap.indent(inputs["other"].history(), " ")}'
        except AttributeError:
            other = represent_as_copypasteable(inputs["other"])
        call = f"{name}{other}"
    else:
        list_of_arguments = [
            f"{k}={represent_as_copypasteable(v)}" for k, v in inputs.items()
        ]
        arguments_as_string = "\n   " + ",\n   ".join(list_of_arguments)

        if "Rainbow" in name:
            call = f"{name}({arguments_as_string})"
        else:
            call = f".{name}({arguments_as_string})"

    return call


def history(self):
    """
    Return a summary of the history of actions that have gone into this `Rainbow`.

    Returns
    -------
    history : str
        A string that does its best to try to summarize
        all the actions that have been applied to this
        `Rainbow` object from the moment it was created.
        In some (but not all) cases, it may be possible
        to copy, paste, and rerun this code to recreate
        the `Rainbow`.
    """

    calls = self.metadata["history"]
    return "(\n" + "\n".join(calls) + "\n)"
