import os

test_directory = "examples/"

try:
    os.mkdir(test_directory)
except FileExistsError:
    pass
