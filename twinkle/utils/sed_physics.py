import os
import numpy as np

# from utils import directories
# from utils import mosaic_tools as mt
#
# from astropy.io import ascii
# from astropy.io import fits
from astropy import constants as con

# import astropy.units as u

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
units_factor = {'angstrom': 1e-8, 'microns': 1e-4,
                'meters': 1e-2, 'cm': 1}


class PhysModels:

    def wienTEMP(self, lambda_, units='micron'):
        r"""
         To calculate maximum temperature of a blackbody at a given
         wavelength using Wien's law.

        Args:
            \lambda_: (float or numpy array) reference wavelength(s) to
                calculate blackbody temperature
            units: (string) describes reference wavelength unit. All values in
                \lambda_ have to be the same. Allowed units are 'angstrom',
                'microns','cm', 'meters'. Units will be converted to cm
                for ease of calculation.

        Returns:
            Temp: (float or np.array) Temperature in Kelvin.
        """
        # CHECK UNITS
        x = lambda_
        try:
            x *= units_factor[units]
        except ValueError:
            raise ValueError('Unit was not recognized')
        temp = _WIEN / x

        return temp

    def blackbody(self, lambda_, p0, su2ea1=1, bands=None,
                  units='angstrom', bulk=False, **kwargs):
        r"""
        Calculates the blackbody irradiance in cgs units (erg s⁻¹ cm⁻² Å⁻¹).

        This function computes the blackbody radiation spectrum for a dust source
        with specified temperature and radius, assuming projected emission with
        a solid angle of π. The input wavelength is converted to centimeters
        for calculations, and the physical constants are used accordingly.

        Args:
            lambda_ (array-like):
                Wavelength array at which the blackbody flux is calculated.
                Units can be specified using the `units` parameter.
            p0 (list or array-like):
                System parameters [Td, Rd]:
                - Td (float): Dust temperature in Kelvin.
                - Rd (float): Dust radius in astronomical units (AU).
            su2ea1 (float, optional):
                Partial flux normalization factor. It accounts for the distance
                to the star and assumes a 1 AU radius:
                [(Rd[AU] / dist[pc]) in natural units: ((AU to cm) / (distance in pc x pc to cm))²]
                Default is 1.
            bands (array-like, optional):
                Spectral bands over which to integrate the blackbody spectrum.
                Default is None.
            units (str, optional):
                Units of the input wavelength. Supported values are:
                - 'angstrom' (default)
                - 'micron'
                - 'meters'
                - 'cm'
            bulk (bool, optional):
                If True, performs bulk processing for efficiency. Default is False.
            kwargs (dict):
                Additional keyword arguments passed to internal calculations.

        Returns:
            numpy.ndarray:
                The blackbody flux at each wavelength in units of erg s⁻¹ cm⁻² Å⁻¹.

        Notes:
            - The function assumes a projected emission with a solid angle of π.
            - The input wavelength is internally converted to centimeters to carry out calculations.
            - This module is designed for dust emission calculations, particularly useful in astrophysical contexts.
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
