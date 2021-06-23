# chromatic
Tools for visualizing spectroscopic light curves, with flux as a function of wavelength and time.

This is *super* in development right now!

## Installation
You should be able to install this by running
```
pip install git+https://github.com/zkbt/chromatic.git
```
from a UNIX prompt.

If you want to be able to modify the code yourself, please also feel free to fork/clone this repository onto your own computer and install directly from that editable package. For example, this might look like:
```
git clone https://github.com/zkbt/chromatic.git
cd chromatic
pip install -e '.[develop]'
```
This will link the installed version of the `chromatic` package to your local repository. Changes you make to the code in the repository should be reflected in the version Python sees when it tries to `import chromatic`. Including the `[develop]` will install the dependencies for the package itself, as well as the extra dependencies required for development (for testing or writing documentation). 

## Usage

The following snippet of code shows the basic structure and functionality of this code (so far):
```python
from chromatic import *

# create a simulated spectroscopic light curve
s = SimulatedRainbow(R=50,
                     dt=10*u.minute,
                     signal_to_noise=100)

# access some basic attributes
print('The wavelengths are...')
print(s.wavelength)
print('The times are...')
print(s.time)
print('The flux (as a function wavelength and time) is...')
print(s.flux)

# bin in both time and wavelength
b = s.bin(dt=0.5*u.hour, dw=0.5*u.micron)

# make a plot showing unbinned + binned flux
fi, ax = plt.subplots(2, 1, sharex=True)
imshowkw = dict( vmin=0.98, vmax=1.02)
s.imshow(ax=ax[0], **imshowkw)
plt.title('Unbinned')
b.imshow(ax=ax[1], **imshowkw)
plt.title('Binned')
plt.tight_layout()
plt.show()
```

## Contributors

This package is being developed during James Webb Space Telescope Early Release Science [Pre-Launch Data Hackathon](https://ers-transit.github.io/pre-launch-hackathon.html). Contributors who agree to follow the [Code of Conduct](https://ers-transit.github.io/code-of-conduct.html#ers-transit) are welcome to join.

- [Zach Berta-Thompson](https://github.com/zkbt)
- ...
