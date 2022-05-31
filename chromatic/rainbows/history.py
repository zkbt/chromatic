"""
Add ability for Rainbow to keep track of their own history.
"""
from ..imports import *

__all__ = [
    "_record_history_entry",
    "_remove_last_history_entry",
    "_create_history_entry",
    "history",
]


def _record_history_entry(self, h):
    """
    Record a history entry in the log.

    Parameters
    ----------
    h : dict
        A dictionary containing keys `name` as the name
        of the action and `inputs` as its keyword inputs.
    """
    try:
        self.metadata["history"].append(h)
    except KeyError:
        self.metadata["history"] = [h]


def _remove_last_history_entry(self):
    """
    Remove the most recent history entry from the log
    """
    try:
        self.metadata["history"].pop()
    except KeyError:
        self.metadata["history"] = []


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
            ",".join([f"{k}={repr(v)}" for k, v in h["inputs"].items()])
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