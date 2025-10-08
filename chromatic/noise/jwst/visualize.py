from ...imports import *


def plot_pandexo(t):
    """
    Make a quick summary plot of PandExo tabular results.
    """
    plt.figure(figsize=(8, 4))
    plt.title(
        f"{t.meta['pandexo_input']['Instrument']}+{t.meta['pandexo_input']['Mode']}+{t.meta['pandexo_input']['Disperser']} | PandExo"
    )
    for k in t.colnames:
        if "snr" in k:
            plt.plot(
                t["wavelength"],
                t[k],
                label=f"{k}\n(median S/N = {np.median(t[k]):.1f})",
            )
    plt.yscale("log")
    plt.legend(frameon=False, bbox_to_anchor=(1, 1))
    plt.xlabel("Wavelength (microns)")
    plt.ylabel("S/N per extracted pixel")

    N_groups = t.meta["number_of_groups_per_integration"]
    summary = f"""
    {t.meta['number_of_groups_per_integration']} groups/integration
    ({N_groups}-1)/({N_groups}+1) = {t.meta['observing_efficiency']:.2%} duty cycle
    {t.meta['time_per_group']}s/group
    {t.meta['time_per_integration']}s/integration
    {t.meta['number_of_integrations_per_transit']} integrations/transit
    {len(t)} wavelengths"""
    plt.text(1, 0, summary, transform=plt.gca().transAxes, va="bottom")


def plot_etc(t):
    """
    Make a quick summary plot of ETC tabular results.
    """
    plt.figure(figsize=(8, 4))
    i = t.meta["configuration"]["instrument"]
    plt.title(f"{i['instrument']}+{i['aperture']}+{i['disperser']} | ETC")
    for k in t.colnames:
        if "snr" in k:
            plt.plot(
                t["wavelength"],
                t[k],
                label=f"{k}\n(median S/N = {np.median(t[k]):.1f})",
            )
    plt.yscale("log")
    plt.legend(frameon=False, bbox_to_anchor=(1, 1))
    plt.xlabel("Wavelength (microns)")
    plt.ylabel("S/N per extracted pixel")

    N_groups = t.meta["configuration"]["detector"]["ngroup"]
    summary = f"""
    {N_groups} groups/integration
    ({N_groups}-1)/({N_groups}+1) = {(N_groups-1)/(N_groups+1):.2%} duty cycle
    {'?'}s/group
    {'?'}s/integration
    {len(t)} wavelengths"""
    # {t.meta['configuration']['detector']['nint']} integrations/exposure
    plt.text(1, 0, summary, transform=plt.gca().transAxes, va="bottom")


def plot_etc_and_pandexo(t_etc, t_pandexo):
    """
    Compare S/N estimates between the official ETC and Pandexo.

    Parameters
    ----------
    t_etc : Table
        The 1D tabular output from `read_etc`
    t_pandexo : Table
        The 1D tabular output from `read_pandexo`
    """
    plt.figure(figsize=(8, 4))
    for t, l in zip([t_etc, t_pandexo], ["ETC", "PandExo"]):
        for k in t.colnames:
            if "snr_per_integration" in k:
                plt.plot(t["wavelength"], t[k], label=f"{k} from {l}")

    plt.xlabel("Wavelength (microns)")
    plt.ylabel("S/N per extracted pixel per integration")
    plt.legend(frameon=False, bbox_to_anchor=(1, 1))
    # plt.yscale('log')
