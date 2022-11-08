# Reference

::: chromatic.rainbows.read_rainbow

::: chromatic.rainbows.rainbow.Rainbow
    selection:
      members:
        - Rainbow
        - __init__
        - __getattr__
        - __setattr__
        - __getitem__
::: chromatic.rainbows.withmodel.RainbowWithModel
::: chromatic.rainbows.simulated.SimulatedRainbow

## 🌈 Helpers
::: chromatic.rainbows.helpers.get
::: chromatic.rainbows.helpers.help
::: chromatic.rainbows.helpers.history
::: chromatic.rainbows.helpers.save

## 🌈 Actions
::: chromatic.rainbows.actions.align_wavelengths
::: chromatic.rainbows.actions.attach_model
::: chromatic.rainbows.actions.binning
    options:
      show_root_heading: False
::: chromatic.rainbows.actions.compare
::: chromatic.rainbows.actions.flag_outliers
::: chromatic.rainbows.actions.fold
::: chromatic.rainbows.actions.inflate_uncertainty
::: chromatic.rainbows.actions.inject_noise
::: chromatic.rainbows.actions.inject_outliers
::: chromatic.rainbows.actions.inject_spectrum
::: chromatic.rainbows.actions.inject_systematics
::: chromatic.rainbows.actions.inject_transit
::: chromatic.rainbows.actions.normalization
    options:
      show_root_heading: False
::: chromatic.rainbows.actions.operations
    options:
      show_root_heading: False
    selection:
      members:
        - __add__
        - __sub__
        - __mul__
        - __truediv__
        - __eq__
::: chromatic.rainbows.actions.remove_trends
::: chromatic.rainbows.actions.shift
::: chromatic.rainbows.actions.trim


## 🌈 Get/Timelike
::: chromatic.rainbows.get.timelike.average_lightcurve
    options:
      show_root_heading: False
::: chromatic.rainbows.get.timelike.median_lightcurve
    options:
      show_root_heading: False
::: chromatic.rainbows.get.timelike.subset
    options:
      show_root_heading: False
::: chromatic.rainbows.get.timelike.time
    options:
      show_root_heading: False

## 🌈 Get/Wavelike
::: chromatic.rainbows.get.wavelike.average_spectrum
    options:
      show_root_heading: False
::: chromatic.rainbows.get.wavelike.expected_uncertainty
    options:
      show_root_heading: False
::: chromatic.rainbows.get.wavelike.measured_scatter_in_bins
    options:
      show_root_heading: False
::: chromatic.rainbows.get.wavelike.measured_scatter
    options:
      show_root_heading: False
::: chromatic.rainbows.get.wavelike.median_spectrum
    options:
      show_root_heading: False
::: chromatic.rainbows.get.wavelike.spectral_resolution
    options:
      show_root_heading: False
::: chromatic.rainbows.get.wavelike.subset
    options:
      show_root_heading: False

## 🌈 Visualizations
::: chromatic.rainbows.visualizations.animate
    options:
      show_root_heading: False
::: chromatic.rainbows.visualizations.colors
    options:
      show_root_heading: False
::: chromatic.rainbows.visualizations.imshow
::: chromatic.rainbows.visualizations.interactive
    options:
      show_root_heading: False
::: chromatic.rainbows.visualizations.pcolormesh
::: chromatic.rainbows.visualizations.plot_lightcurves
::: chromatic.rainbows.visualizations.plot_spectra
::: chromatic.rainbows.visualizations.plot

## 🔨 Tools
::: chromatic.spectra.planck.get_planck_photons
::: chromatic.spectra.phoenix.get_phoenix_photons
::: chromatic.resampling
    options:
      show_root_heading: False
    selection:
      members:
        - bintoR
        - bintogrid
        - resample_while_conserving_flux
::: chromatic.imports
    options:
      show_root_heading: False
