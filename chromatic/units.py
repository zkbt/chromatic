import astropy.units as u

MJy_sr = u.def_unit('MJy/sr', u.MJy/u.sr)
MJy_sr_sq = u.def_unit('(MJy/sr)^2', (u.MJy/u.sr)**2)
u.add_enabled_units([MJy_sr,MJy_sr_sq])


