# chromatic
Tools for visualizing spectroscopic light curves, with flux as a function of wavelength and time. Read the 🌈[documentation](https://zkbt.github.io/chromatic/)🌈 to see how it works!

It's being developed in support of the JWST Transiting Exoplanet Community Early Release Science Program ([ers-transit](https://ers-transit.github.io/)) and easier multiwavelength observations of transiting exoplanets in general, from telescopes in space or on the ground. This package is still actively being developed. Please submit Issues for bugs you notice, features that aren't clearly explained in the documentation, or functionality you'd like to see us implement.

## Installation
If you want to install this code just to use it, you can simply run

```
pip install chromatic-lightcurves
```

and it should install everything, along with all the dependencies needed to run the code. If you previously installed this package and need to grab a newer version, you can run

```
pip install --upgrade chromatic-lightcurves
```
For Developer Installation instructions, see the [documentation](https://zkbt.github.io/chromatic/installation/).

## Usage

For an ultra-quick start try
```python
from chromatic import *
r = SimulatedRainbow().inject_transit().inject_noise()
r.normalize().bin(dw=0.5*u.micron, dt=15*u.minute).imshow()
```
and then see the 🌈[documentation](https://zkbt.github.io/chromatic/)🌈  for more.


## Contributing

We welcome contributions from anyone who agrees to follow the `ers-transit` [Code of Conduct](https://ers-transit.github.io/code-of-conduct.html#ers-transit). If you're on the `ers-transit` slack, please join the #hack-chromatic channel there and say hello; otherwise, please contact Zach directly or just dive right in!

A great initial way to contribute would be to [submit an Issue](https://github.com/zkbt/chromatic/issues) about a bug, question, or suggestion you might have. If you want to contribute code, the [Developer Guide](https://zkbt.github.io/chromatic/designing/) is probably the best place to start. We know it can feel a little scary to try to contribute to a shared code package, so we try our best to be friendly and helpful to new contributors trying to learn how!

*And for context, Zach is a little new to trying to manage a big collaborative code project, so if there are things we could be doing better, please let him know!*

The approximate goal is currently (as of May 2022) to have the code ready and the documentation complete enough to submit `chromatic-lightcurves` to the [Journal of Open Source Software](https://joss.theoj.org/) early enough to be a documented tool in support of the real ERS data (so, like, early late June 2022). If you contribute before then, you'll be included on the paper!
