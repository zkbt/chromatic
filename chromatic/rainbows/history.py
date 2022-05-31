"""
Add ability for Rainbow to keep track of their own history.
"""
from ..imports import *

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

    # create history entry dictionary
    h = dict(name=name, inputs=inputs)
    return h


def history(self, format="string"):
    """
    Return a summary of the history of actions that have gone into this object.

    Parameters
    ----------
    format : str
        How should the history be returned? Options include...
            'string' = a mostly copy-pastable string
            'table' = an ordered table of each individual entries

    Returns
    -------
    history : ?
        The history of actions, in the requested format

    """

    # create a dictionary to store function call components
    d = dict(names=[], arguments=[])
    for h in self.metadata["history"]:
        d["names"].append(h["name"])
        d["arguments"].append(
            ", ".join(
                [f"{k}={represent_as_copypasteable(v)}" for k, v in h["inputs"].items()]
            )
        )

    # create a table of inputs
    table = Table(d)
    if format == "table":
        return table
    elif format == "string":
        calls = [f"{row['names']}({row['arguments']})" for row in table]
        return ".".join(calls)
    else:
        raise ValueError(f"`format='{format}'` not available")
