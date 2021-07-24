# Installation

Welcome! For installing this code we assume you have a Python environment set up, into which you can install packages via `pip`. If this isn't the case, we recommend installing the Anaconda `python` distribution, and using `conda` to manage the version(s) of `python` you have installed on your computer. One (of many) tutorials about how to get started with Python is available [here](https://github.com/ers-transit/hackathon-2021-day0).

## Basic Installation

If you want to install this code just to use it, you can simply run
```
pip install chromatic-lightcurves
```
and it should install everything, along with all the dependencies needed to the code. If you previously installed this package and need to grab a newer version, you can run
```
pip install --upgrade chromatic-lightcurves
```

## Developer Installation

If you want to install this code while being able to edit and develop it, you can fork and/or clone its [github repository](https://github.com/zkbt/chromatic.git) onto your own computer and then install it directly as an editable package its local directory by running
```
git clone https://github.com/zkbt/chromatic.git
cd chromatic
pip install -e '.[develop]'
```
This will point your environment's `chromatic` package to your local folder, meaning that any changes you make in the repository will be reflected what Python sees when it tries to `import chromatic`. Including the `[develop]` after the `.` will install both the dependencies for the package itself and the extra dependencies required for development (= testing and documentation).

You can quickly test whether your installation worked, and what version you have, by running the Python code
```python
import chromatic
chromatic.version()
```

Happy `chromatic`-ing!
