import os, pytest
import matplotlib.pyplot as plt

test_directory = "examples/"
plt.ioff()

try:
    os.mkdir(test_directory)
except FileExistsError:
    pass
