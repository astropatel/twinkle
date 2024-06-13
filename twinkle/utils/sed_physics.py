

import os, re, sys, operator
import glob, string
import numpy as np

from utils import directories
from utils import mosaic_tools as mt
import scipy.interpolate as intp
import scipy.integrate as sintp

try:
    from astropy.io import ascii
    from astropy.io import fits
    from astropy import constants as con
    import astropy.units as u
except ImportError:
    print('Ummmm... Astropy doesnt seem to be installed. Well, that sucks for you.')

__author__ = 'Rahul I. Patel <ri.patel272@gmail.com>, Joe Trollo'



opj = os.path.join
# SET UP CONSTANTS
_CS = con.c.to('cm/s').value
_WIEN = con.b_wien.to('K cm').value
_H = con.h.to('erg s').value
_KB = con.k_B.to('erg/K').value

# SET UP UNIT CONVERSION.
_CM2ANG = 100000000.0
_ANG2CM = 1e-8

class PhysModels:

    def wienTEMP(self, lambda_, units='micron'):
        """
         To calculate maximum temperature of a blackbody at a given
         wavelength using Wien's law.

        Parameters:
        -----------
        lambda_: (float or np.array) reference wavelength(s) to calculate
                                    blackbody temperature
        units: (string) describes reference wavelength unit. All values in lambda_
            have to be the same. Allowed units are 'angstrom','microns','cm'
            'meters'. Units will be converted to cm for ease of calculation.

        Returns:
        --------
        Temp: (float or np.array) Temperature in Kelvin.
        """
        # CHECK UNITS
        x = lambda_
        if units == 'angstrom':
            x *= 1e-8
        elif units == 'microns':
            x *= 1e-4
        elif units == 'meters':
            x *= 1e-2
        elif units == 'cm':
            x = x
        else:
            raise ValueError('Unit was not recognized')
        Temp = _WIEN / x
        return Temp

    def blackbody(self, lambda_, p0, su2ea1=1, bands=None,
                  units='angstrom', bulk=False, **kwargs):
        """
        To calculate the blackbody irradiance in units of cgs
        Flux (erg s-1 cm-2 A-1). This module assumes projected
        emission from a source, so solid angle = pi. All physical
        constants are in units to facilitate this. Wavelength must
        be in units of angstrom, microns, cm, or meters. They will
        be converted to cm to carry out the calculation.
        units: string indicating units of input wavelength
        will be converted to cm.

                   args: 'angstrom','micron','meters','cm'
            p0: Array, or list of system parameters: [Td, Rd]
                Td: dust temperature in Kelvin
                Rd: dust radius in AU
            su2ea1: Partial flux normalization. Incorporates distance
                    to the star, and assumes 1 AU radius: (Rd[AU]/dist[pc])
                    in natural units. ((AU2cm)/(distance*pc2cm))**2

            returns flux of dust
        """

        # CHECK UNITS
        x = lambda_.copy().astype('float32')
        if units == 'angstrom':
            x *= 1e-8
        elif units == 'microns':
            x *= 1e-4
        elif units == 'meters':
            x *= 1e-2
        elif units == 'cm':
            x = x
        else:
            raise ValueError('Unit was not recognized')

        temp0 = p0[0]
        # CHECK PARAMETERS AND NsORMALIZATION
        try:
            su2ea = (p0[1] ** 2) * su2ea1
        except IndexError:
            su2ea = su2ea1
        try:
            bulk = kwargs['bulk']
        except KeyError:
            bulk = False
        # IF CONVOLVING TO FILTER TRANSMISSION, X SHOULD CONFORM TO LENGTH
        # TRANSMISSION CURVE
        if bulk:
            const = su2ea * ((2 * _H * _CS ** 2))
            fluxbb = const / ((x ** 4 * ((np.exp(_H * _CS / (_KB * np.outer(x, temp0))) - 1).transpose())).transpose())

        elif not bulk:
            const = su2ea * ((2 * _H * _CS ** 2) / x ** 4)
            fluxbb = const / (np.exp(_H * _CS / (x * _KB * temp0)) - 1)

            # IF FILTER IS GIVEN, CALCULATES FLUX FROM THROUGH GIVEN
            # BANDPASS. OTHERWISE RETURNS FLUX ARRAY FROM ALL WAVELENGHTS
            # GIVEN IN lambda_
        flux_arr = np.array([])
        if bands is not None:
            for band in bands:
                # print(band,'bbcalc')
                pband = eval('self.' + band + 'pband')
                xmin, xmax = x.min() * _CM2ANG, x.max() * _CM2ANG  # lambda_.min(), lambda_.max()
                pb_xmin, pb_xmax = pband.wavelength.min(), pband.wavelength.max()
                if xmin > pb_xmin:
                    raise ValueError('Sample size too small. Need more on blue end for %s band.' % band)
                elif xmax < pb_xmax:
                    raise ValueError('Sample size too small. Need more on red end %s band.' % band)
                else:
                    pass
                flux = np.array(self.rsr_flux(pband, x * _CM2ANG, fluxbb))
                flux_arr = np.append(flux_arr, flux / pband.pivotWavelength())
        else:
            flux_arr = fluxbb / (x * _CM2ANG)
        return flux_arr


