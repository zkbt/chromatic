import astropy.units as u

MJy_sr = u.def_unit("MJy/sr", u.MJy / u.sr)
MJy_sr_sq = u.def_unit("(MJy/sr)^2", (u.MJy / u.sr) ** 2)
u.add_enabled_units([MJy_sr, MJy_sr_sq])


data_number_per_second = u.def_unit("DN/s")
data_number_per_second_sq = u.def_unit("(DN/s)^2", data_number_per_second**2)
for kludge_DN_unit in [data_number_per_second, data_number_per_second_sq]:
    try:
        u.add_enabled_units([kludge_DN_unit])
    except ValueError:
        pass

electrons_per_group = u.def_unit("Electrons/group")
u.add_enabled_units([electrons_per_group])

# FIXME: some of these feel super kuldgy; how do we make this
# interact more smoothly with units already defined in astropy
# (for example, DN should work, but is always flaky for me...)

u.add_enabled_aliases(
    {"microns": u.micron, "ELECTRONS": u.electron, "electrons": u.electron}
)
