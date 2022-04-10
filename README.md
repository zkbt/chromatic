# chromatic
Tools for visualizing spectroscopic light curves, with flux as a function of wavelength and time. Read the 🌈[documentation](https://zkbt.github.io/chromatic/)🌈 to see how it works! 

It's being developed in support of Webb Transitng Exoplanet Community Early Release Science Program ([ers-transit](https://ers-transit.github.io/)) and easier multiwavelength observations of transiting exoplanets, in general, from telescopes in space or on the ground. This package is still activately being developed. Please submit Issues for bugs you notice, features that aren't clearly explained in the documentation, or functionality you'd like to see us implement.

## Installation
### Basic Installation

If you want to install this code just to use it, you can simply run

```
pip install chromatic-lightcurves
```

and it should install everything, along with all the dependencies needed to run the code. If you previously installed this package and need to grab a newer version, you can run

```
pip install --upgrade chromatic-lightcurves
```

### Developer Installation

If you want to install this code while being able to edit and develop it, you can fork and/or clone its github repository onto your own computer and then install it directly as an editable package its local directory by running

```
git clone https://github.com/zkbt/chromatic.git
cd chromatic
pip install -e '.[develop]'
```

This will point your environment's chromatic package to your local folder, meaning that any changes you make in the repository will be reflected what Python sees when it tries to import chromatic. Including the [develop] after the . will install both the dependencies for the package itself and the extra dependencies required for development (= testing and documentation).

## Contributing

If you want to contribute to this project, especially as part of the ERS Pre-Launch Hackathon, please join the #hack-chromatic channel on the ERS slack.

To be careful we don't mess with each other's stuff, let's *please* avoid doing any code development on the `main` branch or pushing any code directly to `main`. If you want to work on something to contribute, please create a new branch, develop in that branch, and the submit a pull request to merge that branch back into `main`. Christina Hedges' [Astronomy Workflow](https://christinahedges.github.io/astronomy_workflow/notebooks/1.0-basics/git-basics.html) gives a great intro if this feels complicated to you!

*And for context, Zach is a little new to trying to manage a big collaborative code project, so if there are things we could be doing better, please let him know!*

The approximate goal is currently (as of March 2022) to have the code ready and the documentation complete enough to submit `chromatic-lightcurves` to the [Journal of Open Source Software] early enough to be a documented tool in support of the real ERS data (so, like, early summer 2022). If you contribute before then, you'll be included on the paper!

## Usage

See the 🌈[documentation](https://zkbt.github.io/chromatic/)🌈  for example usage.

## Contributors

This package started being developed during James Webb Space Telescope Early Release Science [Pre-Launch Data Hackathon](https://ers-transit.github.io/pre-launch-hackathon.html). Contributors who agree to follow the [Code of Conduct](https://ers-transit.github.io/code-of-conduct.html#ers-transit) are welcome to join. If you're on the ers-transit slack, please join the #hack-chromatic channel there and say hello; otherwise, please contact Zach directly or just dive right in!
