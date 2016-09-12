import wise_flux_correction
import sed_paramfile as sp

mags4Dust0 = sp.mags4DustWFC
wfc = wise_flux_correction.WISECorrect(bands=mags4Dust0)