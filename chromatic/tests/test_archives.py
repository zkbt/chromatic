from ..archives import *
from .setup_tests import *


def test_download_from_mast():
    d = download_from_mast()
    assert d[0]["Status"] == "COMPLETE"
