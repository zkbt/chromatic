from ...imports import *


def extract_sn_from_image(images, extraction_box=3):
    """
    Use a 2D image to estimate per-wavelength S/N.

    (This is tuned for MIRI/LRS. It probably won't work on anything else.)

    Returns
    -------
    sn_extracted : array
        The 1D extracted S/N.
    """

    # extract S/N and saturation images
    sn = images["snr"]
    saturation = images["saturation"]

    # define extraction aperture
    y_pix = np.arange(sn.shape[0]) - sn.shape[0] / 2 + 0.5
    extent = [0, sn.shape[1] - 1, -sn.shape[0] / 2 + 0.5, sn.shape[0] / 2 + 0.5]
    extract = np.abs(y_pix) < extraction_box
    extract_alpha = 0.5
    print(f"{np.sum(extract)} pixels being used in the extraction")

    # display the images
    fi, ax = plt.subplots(
        2,
        2,
        figsize=(7.5, 3),
        gridspec_kw=dict(width_ratios=[8, 1]),
        sharey="row",
        sharex="col",
    )

    plt.sca(ax[0, 0])
    plt.imshow(np.log10(sn), extent=extent, cmap="gray")

    sat_color = "red"
    sat_cmap = one2another(bottom=sat_color, top=sat_color, alpha_bottom=0, alpha_top=1)
    plt.imshow(saturation, extent=extent, cmap=sat_cmap, vmax=2)
    for sign in [-1, 1]:
        plt.axhline(sign * extraction_box, alpha=extract_alpha, color="gray")

    plt.sca(ax[0, 1])
    spatial_profile = np.nanmean(sn, axis=1)
    plt.plot(spatial_profile, y_pix, color="black")
    for sign in [-1, 1]:
        plt.axhline(sign * extraction_box, alpha=extract_alpha, color="gray")

    sn_extracted = np.sqrt(np.sum(sn**2, axis=0))
    plt.sca(ax[1, 0])
    plt.plot(sn_extracted, color="black")
    plt.ylabel("S/N per pixel\nper integration")

    for i, label in zip([1, 2], ["partial", "full"]):
        is_saturated = np.nonzero(np.max(saturation, axis=0) >= i)[0]
        if len(is_saturated) > 0:
            left, right = np.min(is_saturated), np.max(is_saturated)
            plt.axvspan(
                left,
                right,
                color=sat_color,
                alpha=0.25 * i,
            )
        print(f"{len(is_saturated)} pixels experience at least {label} saturation")
    return sn_extracted


def set_image_slice(instrument):
    if instrument.lower() == "p750l":
        image_slice = slice(27, 399)
    else:
        image_slice = slice(None)
    return image_slice
