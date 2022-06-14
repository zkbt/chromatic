import os, pytest
import matplotlib.pyplot as plt

test_directory = "examples/"

try:
    os.mkdir(test_directory)
except FileExistsError:
    pass
