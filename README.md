# chromatic
Tools for visualizing spectroscopic light curves, with flux as a function of wavelength and time.

This is *super* in development right now!

## Installation
If you want to install this code just to use it, you can simply run
```
pip install git+https://github.com/zkbt/chromatic.git
```

If you want to install this code while being able to edit and develop it, you can fork and/or clone this repository onto your own computer and then install it directly as an editable package by running
```
git clone https://github.com/zkbt/chromatic.git
cd chromatic
pip install -e '.[develop]'
```
This will point your environment's `chromatic` package to point to this local folder, meaning that any changes you make in the repository will be reflected what Python sees when it tries to `import chromatic`. Including the `[develop]` will install both the dependencies for the package itself and the extra dependencies required for development (= testing and documentation).

## Contributing

If you want to contribute to this project, especially as part of the ERS Pre-Launch Hackathon, please join the #hack-chromatic channel on the ERS slack.

To be careful we don't mess with each other's stuff, let's *please* avoid doing any code development on the `main` branch or pushing any code directly to `main`. If you want to work on something to contribute, please create a new branch, develop in that branch, and the submit a pull request to merge that branch back into `main`. Christina Hedges' [Astronomy Workflow](https://christinahedges.github.io/astronomy_workflow/notebooks/1.0-basics/git-basics.html) gives a great intro if this feels complicated to you!

*And for context, Zach is a little new to trying to manage a big collaborative code project myself, so if there are things we could be doing better, please let him know!*

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
