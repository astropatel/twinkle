"""
==================================================================================
 sed.py by Rahul I. Patel (ri.patel272@gmail.com)
==================================================================================
  Code with classes and definitions for useful SED plotting/calculations. Right
  now it only supports WISE, 2MASS, Johnson/Bessel, Spitzer/MIPS60,70/160 RSR's,
  although individual bandpasses can be called on using create_passband. Contains
  modules to calculate simple blackbody and photospheric grid fluxes.
  Parts of class Bandpass were taken from astLib with permission
  from Matt Hilton (http://astlib.sourceforge.net/). Please cite as necessary.

  Issues and changes that need to be made:
  1. Remove mpfit and use internal scipy fitting routine.
  2. vega2AB needs to be updated or deprecated.
  3. Add logger -- DONT' NEED FOR THIS ONE
  4. remove readcol.py -- DONE, Tested
  5. reamove astro_tools
  7. Change init to only load pband files from mags2use.
  8. Calc_temp needs better docstring
"""

import os, re, operator
import glob, string
import numpy as np

from utils import directories
from utils import mosaic_tools as mt
import scipy.interpolate as intp
import scipy.integrate as sintp

try:
    from astropy import constants as con
except ImportError:
    print 'Does not seem Astropy is installed, or at least the constants package is messed up. We kinda need this. Get to it yo.'
try:
    from astropy.io import ascii
    from astropy.io import fits
except ImportError:
    print 'Ummmm... Astropy doesnt seem to be installed. Well, that sucks for you.'

__author__ = 'Rahul I. Patel <ri.patel272@gmial.com>, Joe Trollo'

DIR = directories
opj = os.path.join
intpdir = DIR.SupportFiles()
fRSR = DIR.RSR()
AT = mt.ArrayTools()
FT = mt.FittingTools()

# SET UP CONSTANTS
_CS = con.c.to('cm/s').value
_WIEN = con.b_wien.to('K cm').value
_H = con.h.to('erg s').value
_KB = con.k_B.to('erg/K').value

# SET UP UNIT CONVERSION.
_CM2ANG = 100000000.0
_ANG2CM = 1e-8

#  DICTIONARY TO HOUSE ALL THE PHOTOSPHERIC MODELS EVENTUALLY
MegaGrid = {}

#  DICTIONARY THAT CONTAINS ALL THE EMPIRICAL STELLAR COLOR DATA
EmpDat = {}

#  DICTIONARY THAT CONTAINS ALL DATA FROM A DATA INPUT FILE
StarsDat = {}


class SEDTools:
    """SEDTools has functions to aid in modeling an SED of an astronomical source
        -- creates bandpass objects, converts stellar magnitudes to flux in cgs units,
           convolves flux with said bandpasses and calculates (for now) simple blackbody
           and stellar photosphere from Grids. Amongst a bunch of other things.
    """

    def __init__(self):

        # intpdir = os.path.join(os.getcwd(), 'Interpolation_Files')
        # fRSR = os.path.join(intpdir, 'RSR')
        self.create_passbands()

    def create_passbands(self, RSRFile=None, flat=False, waverange=(None, None), cntr=None):
        """
         Returns passband objects (no input) for all bands in WISE, 2MASS, and
         Johnson UBVRI. All passband info need to be in '.../Interpolation_Files/RSR/
         directory.

         You can serparately enter the RSRFile name along with wavelength units
         to create your own passband object. Above rules apply.
         Wavelength is automatically converted to Angstroms.

         Parameters:
         -----------
         RSRFile : str
            File name that contains the RSR information
         flat : bool
            When set True, will let you create a flat spectrum
         waverange : tuple
            2 float elements. These need to be initiated to the
                    min and max wave of the flat spectrum
         cntr: center or isophotal wavelength of the flat bandpass.
         """

        units = {'WISE': 'microns', '2MASS': 'microns', 'Johnson': 'angstroms',
                 'MIPS': 'microns', 'Tycho': 'angstroms', 'IRAS': 'angstroms',
                 'HPACS': 'angstroms', 'Akari':'angstroms','LNIRC2':'microns'}

        aspband = Bandpass
        if RSRFile is None and not flat:

            #  MICRONS
            u = units['WISE']
            self.W1pband = aspband(opj(fRSR, 'W1_WISE.dat'), inputUnits=u)
            self.W2pband = aspband(opj(fRSR, 'W2_WISE.dat'), inputUnits=u)
            self.W3pband = aspband(opj(fRSR, 'W3_WISE.dat'), inputUnits=u)
            self.W4pband = aspband(opj(fRSR, 'W4_WISE.dat'), inputUnits=u)
            self.W4truncpband = aspband(opj(fRSR, 'W4_WISE_truncated.dat'), inputUnits=u)
            self.W4stretchpband = aspband(opj(fRSR, 'W4_WISE_stretched.dat'), inputUnits=u)

            #  MICRONS
            u = units['MIPS']
            self.MIPS24pband = aspband(opj(fRSR, '24_MIPS.dat'), inputUnits=u)
            self.MIPS70pband = aspband(opj(fRSR, '70_MIPS.dat'), inputUnits=u)
            self.MIPS160pband = aspband(opj(fRSR, '160_MIPS.dat'), inputUnits=u)

            #  MICRONS
            u = units['2MASS']
            self.J2Mpband = aspband(opj(fRSR, 'J_2MASS.dat'), inputUnits=u)
            self.H2Mpband = aspband(opj(fRSR, 'H_2MASS.dat'), inputUnits=u)
            self.Ks2Mpband = aspband(opj(fRSR, 'Ks_2MASS.dat'), inputUnits=u)

            #  ANGSTROMS
            u = units['Johnson']
            self.UBpband = aspband(opj(fRSR, 'U_Bessel.dat'), inputUnits=u)
            self.BBpband = aspband(opj(fRSR, 'B_Bessel.dat'), inputUnits=u)
            self.VBpband = aspband(opj(fRSR, 'V_Bessel.dat'), inputUnits=u)
            self.RBpband = aspband(opj(fRSR, 'R_Bessel.dat'), inputUnits=u)
            self.IBpband = aspband(opj(fRSR, 'I_Bessel.dat'), inputUnits=u)

            self.UJpband = aspband(opj(fRSR, 'U_Johnson.dat'), inputUnits=u)
            self.BJpband = aspband(opj(fRSR, 'B_Johnson.dat'), inputUnits=u)
            self.VJpband = aspband(opj(fRSR, 'V_Johnson.dat'), inputUnits=u)
            self.RJpband = aspband(opj(fRSR, 'R_Johnson.dat'), inputUnits=u)
            self.IJpband = aspband(opj(fRSR, 'I_Johnson.dat'), inputUnits=u)

            u = units['Tycho']
            self.VTpband = aspband(opj(fRSR, 'V_Tycho.dat'), inputUnits=u)
            self.BTpband = aspband(opj(fRSR, 'B_Tycho.dat'), inputUnits=u)
            self.Hppband = aspband(opj(fRSR, 'Hp.dat'), inputUnits=u)

            u = units['IRAS']
            self.IRAS60pband = aspband(opj(fRSR, '60_IRAS.dat'), inputUnits=u)
            self.IRAS100pband = aspband(opj(fRSR, '100_IRAS.dat'), inputUnits=u)
            self.IRAS25pband = aspband(opj(fRSR,'25_IRAS.dat'), inputUnits=u)
            self.IRAS12pband = aspband(opj(fRSR,'12_IRAS.dat'), inputUnits=u)

            u = units['HPACS']
            self.HPACS70pband = aspband(opj(fRSR, '70_HPACS.dat'), inputUnits=u)
            self.HPACS100pband = aspband(opj(fRSR, '100_HPACS.dat'), inputUnits=u)
            self.HPACS160pband = aspband(opj(fRSR, '160_HPACS.dat'), inputUnits=u)

            u = units['Akari']
            self.Akari90pband = aspband(opj(fRSR,'90_Akari.dat'), inputUnits=u)

            u = units['LNIRC2']
            self.LpNIRC2pband = aspband(opj(fRSR,'Lp_NIRC2.dat'),inputUnits=u)
            self.HNIRC2pband = aspband(opj(fRSR, 'H_NIRC2.dat'), inputUnits=u)
            self.KpNIRC2pband = aspband(opj(fRSR, 'Kp_NIRC2.dat'), inputUnits=u)

        elif RSRFile is not None and not flat:
            self.pband = aspband(opj(fRSR, RSRFile), inputUnits=u)
        elif flat:
            self.pband = Flatbandpass(waverange, cntr)

    def cgs2Jy(self, flux,dflux,band=None,wave=None,nu=None):
        """To convert specific flux from erg/s/cm2/Ang to Jy. It assumes
        input wavelength is in Angstroms

        Parameters:
        -----------
        nu (flt/array): frequency in Hz to convert to Jansky.
        wave (flt/array): wavelength in Angstroms
        flux (flt/array): flux in erg/s/cm2/Ang

        Return:
        --------
        flux density in Jansky(ies).
        """
        _CS = con.c.to('angstrom/s').value

        if wave is None:
            wave = band.isoWavelength()


        Flux = flux * wave

        if nu is None:
            try:
                nu = band.isoFrequency()

            except UnboundLocalError:

                nu = _CS / wave # ERROR IS INCURRED BY CHOICE OF SPEED OF LIGHT -- USE CAREFULLY
            #  IF POSSIBLE TRY TO USE GIVEN FREQUENCIES
        else: nu = _CS / wave


        # y = flux*(wave**2/(_CS))*1e23
        y = (Flux / nu) * 1e23
        # DOES NOT INCLUDE UNCERTAINITY IN FREQUENCY OR WAVELENGTH
        dy = (1e23 / nu) * dflux * wave # + (flux * dwave)**2)

        return y, dy

    def mag2Jy(self, band,mag,magerr):

        fcgs,dfcgs = self.mag2fluxZPLam(eval('self.{}pband'.format(band)),mag,magerr)

        fjy, efjy = self.cgs2Jy(fcgs,dfcgs,eval('self.{}pband'.format(band)))

        return fjy,efjy

    def Jy2cgs(self, nuFnu, DnuFnu=None):
        """ To convert specific flux in Jansky to erg/s/cm2.
        It assumes the wavelength is in Angstroms.

        Parameters:
        -----------
        nuFnu: tuple(flt/array,flt/array): frequency in Hz at which flux
               is found and flux in Janskies
        DnuFnu: tuple(flt/array,flt/array): error in frequency measurement
                and error in flux measurement

        Return:
        --------
        A tuple with the flux and conversion error (Flux, dFlux) in erg/s/cm2.
        If no errors are given, dFlux = None
        """
        nu, flux = np.asarray(nuFnu[0]), np.asarray(nuFnu[1])  # 2 ELEMENT TUPLE
        y = (flux * nu) / 1e23

        if DnuFnu is None:
            dFlux = None
        else:
            dnu, dfnu = np.asarray(DnuFnu[0]), np.asarray(DnuFnu[1])
            dFlux = np.sqrt((flux * dnu) ** 2 + (nu * dfnu) ** 2) / 1e23

        return y, dFlux

    def get_eff_wavelengths(self, mag_list):
        """
        Returns dictionary of central wavelength's listed in mag_list
        Contingent on passband objects being created in create_passbands.
        Returns isophotal wavelength in RSR headers or calculates
        effective/pivot wavelength if isowavelength is
        unavailable (Tokunanga & Vacca 2005).

        Parameters:
        -----------
        mag_list (list): list of strings with keys corresponding to band from which
                         to extract wavelength. If mag_list = ['mv1','mv2',...,'mvN'],
                         then wave[i] = self.mv1pband.isoWavelength(). Check to see which
                         RSR bands are available.

        Return:
        --------
        waveDict (dict): dictionary of keys:val such that keys=mvi & val are the iso/pivot
                         wavelengths (Angstroms).
         """

        waveDict = {}
        mag2use = mag_list
        for mv in mag2use:
            temp = mv
            waveDict[temp] = eval('self.' + mv + 'pband.isoWavelength()')
        return waveDict

    def vega2AB(self, band, m_vega):
        """
        XXX - This module needs to be updated.
        Module to convert photometric magnitudes into AB mag system.
        AB magnitude is defined such that the monochromatic flux is
        measured in erg/s/cm^2/Hz,
        I.e., AB = -2.5 log(f_cgs) - 48.6 ,
        such that for any bandpass or filter, the zero point mag
        corresponds to a flux density of 3631 Jy.

        Parameters:
        -----------

        system: string. Instrument/basis used for photometry (eg. WISE, 2MASS, Johnson)
                currently supports WISE, 2MASS, Johnson
        band: string. Which pass band. (options are: J,Ks,H,W1,W2,W3,W4,V,R,B,I,Rc,Ic)
        m_vega: flt. or array Vega magnitude
         """
        # TAKEN FROM ODENWALD ET AL. 2003
        # Analysis of the Diffuse Near-Infrared Emission from 2MASS Deep Integration Data:
        # Foregrounds versus the Cosmic Infrared Background

        if band == 'J2M':
            mab = m_vega + 0.894
        elif band == 'H2M':
            mab = m_vega + 1.37
        elif band == 'Ks2M':
            mab = m_vega + 1.84
        # TAKEN FROM WISE FAQ
        # http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/

        elif band == 'W1':
            mab = m_vega + 2.699
        elif band == 'W2':
            mab = m_vega + 3.339
        elif band == 'W3':
            mab = m_vega + 5.174
        elif band == 'W4':
            mab = m_vega + 6.620
        # TAKEN FROM FREI & GUNN 1995
        # http://www.astro.utoronto.ca/~patton/astro/mags.html# conversions
        # AT THE MOMENT IT DOES NOT INCLUDE STROMGREN FILTERS-- TOO LAZY
        elif band == 'VJ':
            mab = m_vega - 0.044
        elif band == 'BJ':
            mab = m_vega - 0.163
        elif band == 'RJ':
            mab = m_vega + 0.055
        elif band == 'IJ':
            mab = m_vega + 0.309
        elif band == 'RJ':
            mab = m_vega + 0.117
        elif band == 'Ic':
            mab = m_vega + 0.342
        else:
            print 'Band ', band, ' doesn" exist. Please check bandpass.'
            mab = None

        if mab == None:
            print 'No conversion occured. Check variables'

        return mab

    def batch_vega2AB(self, mag_list, vegaMag, vegaMagErr):
        """
         Converts all the vega magnitudes in vegaMag
         and associated errors in vegaMagErr to AB magnitudes.

        Parameters:
        -----------
         mag_list (list): List of strings with passband name
         vegaMag,+Err (Dictionary): Vega magnitudes in mag_list

         Return:
         --------
         Dictionaries of converted AB mag and associated errors with keys used from input
         mag_list suffixed by '_AB' or '_ABerr'
         """
        mag2use = mag_list
        vegaMagDict, vegaMagErrDict = vegaMag, vegaMagErr

        ABmag_Dict = {}
        ABmagerr_Dict = vegaMagErrDict

        for mv in mag2use:
            # CONVERT MAG TO AB mag
            temp = mv + '_AB'  # if mv = 'W1', temp = 'W1_AB'
            try:
                ABmag_Dict[temp] = self.vega2AB(mv, vegaMagDict[mv])
            except:
                ABmag_Dict[temp] = 0

            # CONVERT MAG_ERR TO AB mag err --> no change in error
            temp = mv + '_ABerr'
            ABmagerr_Dict[temp] = vegaMagErrDict[mv]

        return [ABmag_Dict, ABmagerr_Dict]

    def mag2fluxAB(self, band, abmags, abmagserr):
        """
        Does the same thing as batch_mag2flux except only takes one band
        and can be applied to multiple measurements in the same band to produce a flux

        """

        fluxJy = (10 ** 23.0) * 10 ** (-(abmags + 48.6) / 2.5)  # converts to AB Mag
        alam = 3e-13  # converts to erg s-1 cm-2 angstrom-1 with lam in microns (for some stupid reason)
        errLmic = band.isoWavelength() * (1e-10 / 1e-6)
        # band.pivotWavelength()
        fluxWLUnits = alam * fluxJy / errLmic ** 2

        fluxJyErr = (10 ** 23.0) * 10 ** (-(abmags - abmagserr + 48.6) / 2.5)
        fluxWLUnitsErr = alam * fluxJyErr / errLmic ** 2
        fluxWLUnitsErr = fluxWLUnitsErr - fluxWLUnits

        return [fluxWLUnits, fluxWLUnitsErr]

    def batch_mag2fluxAB(self, band_list, abmags, abmagserr):
        """Tools to convert from AB mags to fluxes.
          returns a Dictionary of fluxes and flux errors with keys
          from input list.

          Parameters:
          -----------
          band_list (list): strings with magnitude names (ex: ['mv1','mv2',...,'mvn'])
          abmags/abmagserr (dict): dictionary of magnitude/error in AB mag format
                                with keys 'mvi_AB' or 'mvi_ABerr' for errors
          Return:
          -------
          dictionaries of converted AB mag to fluxes and separate dic for
                   flux errors (erg s-1 cm-2 A-1); keys: 'mvi_flux'
       """
        band2use = band_list
        fluxDict, fluxDict_err = {}, {}
        # CALUCLATE FLUX AND ERR IN FLUX AND CREATE DICTIONARIES OF EACH
        for mv in band2use:
            temp = mv + '_flux'
            magAB, magABerr = abmags[mv + '_AB'], abmagserr[mv + '_ABerr']
            fluxList = self.mag2fluxAB(eval('self.' + mv + 'pband'), magAB, magABerr)
            fluxDict[temp], fluxDict_err[temp] = fluxList[0], fluxList[1]
        # THIS CONVERTS AB MAG TO FLUX IN erg/s/cm^2/Angstrom

        return fluxDict, fluxDict_err

    def fluxZPLam2mag(self, band, flux, fluxErr):
        """Convert input flux which should be in units of ergs/cm^2/s/Angstrom to
        vega magnitudes using the zero point fluxes for each bandpass listed in
        the RSR file and hence the passband object

        Parameters:
        -----------
        band: passband object
        flux (float): Flux density in ergs/cm^2/s/Angstrom
        fluxErr: (float) error assosciated wiht the input flux
        
        Output:
        --------
        list of [mag, magerr] in vega magnitudes"""

        Fiso, dFiso = band.fluxVegaZeroPointLam()
        mag_lam = -2.5 * np.log10(flux / Fiso)
        mag_lam_err = (1. / Fiso) * (2.5 / np.log(10)) * fluxErr / 10 ** (mag_lam / -2.5)
        return mag_lam, mag_lam_err

    def mag2fluxZPLam(self, band, Vmag, VmagErr, sysErr=False):
        """Converts the input vega magnitude to specific flux in ergs/cm^2/s/Angstrom using the
        zero point fluxes obtained from the passband file from literature.

        Parameters:
        -----------
        band: passband object
        Vmag: (float) Vega magnitude
        VmagErr: (float) error associated with Vega Magnitude)
        sysErr: (bool) Includes the given bandpasses photospheric calibration systematic
                       uncertainty found in the RSR text file. Set to True if this error
                       needs to be added in quadrature to the total integrated flux
                       uncertainty. Default is False

        Output:
        --------
        list of [flux, fluxerr] in ergs cm^-2 s^-1 A^-1"""

        Fiso, dFiso = band.fluxVegaZeroPointLam()
        Flux_lam = Fiso * (10 ** (Vmag / -2.5))

        df_l1 = ((np.log(10) / -2.5) * Fiso * VmagErr)  # LOG IS NATURAL LOG
        if sysErr:
            df_l2 = dFiso
        else:
            df_l2 = 0

        dFlux_lam = 10 ** (Vmag / -2.5) * np.sqrt(df_l1 ** 2 + df_l2 ** 2)

        return Flux_lam, dFlux_lam

    def mag2fluxZPNu(self, band, Vmag, VmagErr, sysErr=False):
        """Converts the input vega magnitude to specific flux in ergs/cm^2/s/Hz using the
        zero point fluxes obtained from the passband file from literature. For this module,
        the zero point fluxes in frequency need to be there otherwise this won't work.

        Parameters:
        -----------
        band: passband object
        Vmag: (float) Vega magnitude
        VmagErr: (float) error associated with Vega Magnitude)
        sysErr: (bool) Includes the given bandpasses photospheric calibration systematic
                       uncertainty found in the RSR text file. Set to True if this error
                       needs to be added in quadrature to the total integrated flux
                       uncertainty. Default is False


        Output:
        --------
        list of [flux,fluxerr] in ergs cm^-2 s^-1 Hz^-1
        """

        Fiso, dFiso = band.fluxVegaZeroPointFreq()
        Flux_nu = Fiso * (10 ** (Vmag / -2.5))

        df_l1 = ((np.log(10) / -2.5) * Fiso * VmagErr)
        if sysErr:
            df_l2 = dFiso
        else:
            df_l2 = 0

        dFlux_nu = 10 ** (Vmag / -2.5) * np.sqrt(df_l1 ** 2 + df_l2 ** 2)

        return Flux_nu, dFlux_nu

    def batch_mag2fluxZPLam(self, band_list, Vmags, Vmagserr):
        """Converts mag to flux using mag2fluxZP using zero point fluxs.
        for multiple bands all at once. Uses mag2fluxZPLam.

        Input:
        -------
        band_list: (list) strings with magnitude names (ex:['mv1','mv2',...,'mvn'])
        Vmags/Vmagserr: dictionary of magnitude/error in Vega Mag format with keys 'mvi'
                         for both mag and errors.

        Return:
        --------
        dictionary of converted Vega magnitude to fluxes and associated errors in units of
        ergs/s/cm^2/Angstrom."""

        band2use = band_list
        fluxDict, fluxDict_err = {}, {}
        for mv in band2use:
            temp = mv + '_flux'
            magV, magVerr = Vmags[mv], Vmagserr[mv]
            fluxList = self.mag2fluxZPLam(eval('self.' + mv + 'pband'), magV, magVerr)
            fluxDict[temp], fluxDict_err[temp] = fluxList[0], fluxList[1]

        return (fluxDict, fluxDict_err)

    def rsr_eflux(self, pband, lambda_,flux,flux_up, flux_down,
                  report_type='avg'):
        """
        Calculates the integrated flux filtered through a particluar
        passband, as well as the upper and lower limits of the flux
        to give a 1 sigma uncertainty based on the fit parameters
        uncertainties

        Parameters
        ----------
        pband : passband object for a given photometric band created
              from Bandpass class
        lambda_ :(Array) Wavelength in Angstroms range of flux to
                be calculated. Range should be atleast the extent of pband
        flux: (Array) Flux of main fit
        flux_up : (Array) Flux of upper limit of fit.
        flux_down : (Array) Flux of lower limit of fit.
        report_type : (str) either 'avg' to report average of uncertainties or
                            'both'
        Returns
        -------
        flux at the particular band and uncertainties. If 'avg', returns 2 element array
        (flux, eflux). If 'both', returns 3 element array (flux, eflux_low, eflux_high)

        """

        flux = self.rsr_flux(pband,lambda_,flux)[0]
        flux_low = self.rsr_flux(pband,lambda_,flux_down)[0]
        flux_up = self.rsr_flux(pband,lambda_,flux_up)[0]

        fl, fu = flux - flux_low, flux_up - flux
        if report_type == 'avg':
            return flux,np.average([fl,fu])

        elif report_type == 'both':
            return flux, fl, fu
        else:
            raise Exception('No report_type specified. If you only want integrated flux, use rsr_flux')

    def rsr_flux(self, pband, lambda_, flux):
        """Calculates integreated flux filtered through a passband object

        Parameters:
        -----------

       pband: passband object for a given photometric band created
              from Bandpass class
       lambda_: (Array) Wavelength in Angstroms range of flux to
                be calculated. Range should be atleast the extent of pband
       flux: (Array) Flux mapped to lambda_ in erg s^-1 cm^-2?

       Return:
       -------
       Integrated flux over the specified bandpass as float or array. Integration is
       done using Simpson's rule for integration such that:
       Integrated Flux = Integral(Flux*RSR*lam * dlam)/ Integral(RSR*lam*dlam)
        """

        # FIND THE BOUNDS FOR THE WAVELENGHTS IN THIS PARTICULAR BAND
        # Bounds can't exceed the limits of the bandpass' limits
        # OTHERWISE IT WONT' INTERPOLATE CORRECTLY
        # CHECK TO SEE IF ONLY ONE FLUX OR MULTIPLE FLUXES NEED TO BE CALCULATED
        # AT THE SAME BANDPASS
        flux2 = flux.copy()
        Sn_arr = np.array([])
        band_min, band_max = pband.wavelength.min(), pband.wavelength.max()
        indzero = np.where(pband.transmission <= 0)[0]
        pband.transmission[indzero] = 1e-18
        filter_interp = intp.interp1d(np.log10(pband.wavelength), \
                                      np.log10(pband.transmission))
        if lambda_.ndim < 2:

            RSRminInd, RSRmaxInd = int(np.searchsorted(lambda_, band_min)), \
                                   int(np.searchsorted(lambda_, band_max))
            if (RSRmaxInd - RSRminInd) < 10:
                print 'Resolution less than 10.'
            # Create index array
            Ind = np.arange(RSRminInd, RSRmaxInd)
            # COLLECT ALL THE WAVELENGTH AND FLUX VALUES UNDER BANDPASS
            # FOUND IN MODEL VALUES
            lam0, flx0 = lambda_[Ind], flux2[Ind]
            # FOR EACH WAVELENGTH, CALCUALTE A NEW RSR: INTERPOLATE
            RSRNew = 10 ** filter_interp(np.log10(lam0))

            # Use Simpson's approximation to calculate total flux
            # FIND UNIQUE ELEMENTS IN WAVELENGTH BECAUSE THERE ARE REPEATED VALUES
            #  IN NEXTGEN GRIDS
            u, indu = np.unique(lam0, return_index=True)
            lam0, flx0 = lam0[indu], flx0[indu]
            lam02, RSRNew = lam0, RSRNew

            RSRNew = RSRNew[indu]
            Sn = sintp.simps(flx0 * lam0 * RSRNew, lam0) / sintp.simps(lam0 * RSRNew, lam0)

            Sn_arr = np.append(Sn_arr, Sn)
        else:
            # DETERMINES HOW MANY rows (different temperatures) THERE ARE
            for i in range(len(lambda_)):
                lambda_i, fluxi = lambda_[i], flux2[i]
                RSRminInd, RSRmaxInd = int(np.searchsorted(lambda_i, band_min)), \
                                       int(np.searchsorted(lambda_i, band_max))
                if (RSRmaxInd - RSRminInd) < 10:
                    print 'Resolution less than 10.'
                # Create index array
                Ind = np.arange(RSRminInd, RSRmaxInd)
                # COLLECT ALL THE WAVELENGTH AND FLUX VALUES UNDER BANDPASS
                # FOUND IN MODEL VALUES
                lam0i, flx0i = lambda_i[Ind], fluxi[Ind]
                # FOR EACH WAVELENGTH, CALCUALTE A NEW RSR: INTERPOLATE
                RSRNew = 10 ** filter_interp(np.log10(lam0i))
                # Use Simpson approximation to calculate total flux
                # FIND UNIQUE ELEMENTS IN WAVELENGTH BECAUSE THERE ARE REPEATED VALUES
                #  IN NEXTGEN GRIDS
                u, indu = np.unique(lam0i, return_index=True)
                lam0i, flx0i = lam0i[indu], flx0i[indu]
                lam02i, RSRNew = lam0i, RSRNew  # dset2Newi
                RSRNew = RSRNew[indu]

                Sn = sintp.simps(flx0i * lam0i * RSRNew, lam0i) / sintp.simps(lam0i * RSRNew, lam0i)

                if np.isnan(np.sum(Sn)):
                    print 'Is nan', np.isnan(np.sum(Sn))

                Sn_arr = np.append(Sn_arr, Sn)

        return Sn_arr

    def indSandwhich(self, t0, arr):
        ind = np.searchsorted(arr, t0)
        indRight = ind + 3
        indLeft = ind - 3
        if indRight > len(arr):
            indRight = int(len(arr))
        if indLeft < 0:
            indLeft = 0

        return (indLeft, indRight)

    def calc_grids(self, lambda_, p0, su2ea1=1, griddata=None, tempArr=None,
                   mag2use=None, resample=True, **kwargs):

        """Returns array of flux values at input wavelength array interpolated
          at T0.

          Parameters:
          -----------
          lambda_: (Array,float), wavelength values where flux will be
                    calculated
          T0: (float) Temperature in Kelvins to determine best grid
                      flux
          su2ea1: (float) factor that includes the distance^2 to convert from surface
                   earth flux : (1/radius)^2*(radius^2/dist^2)
          griddata: tuple of (lam_arr,flux_arr) created from get_grids
          tempArr: (np.array) Array of increasing temperature values
                    mapped to those found in griddata
          mag2use: (np.array) Array of strings with names of photometric
                    magnitude bands where fluxes are calculated. Needed
                    to call passband object.

          Returns:
          ---------
          GFluxArr; Array of fluxes calculated at lambda_ in units of
                   grid: Kurucz (ergs cm^{-2} s^{-1} A^{-1})
        """
        temp0 = p0[0]
        try:
            su2ea = (p0[1] ** 2) * su2ea1
        except:
            su2ea = su2ea1
        try:
            resample = kwargs['resample']
        except:
            resample = resample


        mags2use = mag2use
        grids, TempsArr = griddata, tempArr
        indLeft, indRight = self.indSandwhich(temp0, TempsArr)

        lam_arr_all = grids[0][indLeft:indRight + 1]
        flux_arr_all = grids[1][indLeft:indRight + 1]
        xin = TempsArr[indLeft:indRight + 1]
        flux_arr_all = flux_arr_all * su2ea
        GFluxArr = np.array([])
        # TO BE USED IF NO RSR IS NEEDED... INDIVIDUAL FLUX
        # CALUCLATION AT SINGLE WAVELENGTH
        if mag2use is None:
            fnew_all_temp = np.array([])  # WILL BE USED TO STORE CALCULATED FLUXES

            lam_arr_all, flux_arr_all = np.log10(lam_arr_all), np.log10(flux_arr_all)
            lambda_ = np.log10(lambda_)
            # GATHER NEW FLUXES FOR EACH TEMPERATURE AT GIVEN g and met
            for i in range(len(lam_arr_all)):
                ipolate = intp.interp1d(lam_arr_all[i], flux_arr_all[i])
                fnew = 10 ** ipolate(lambda_)
                if len(fnew_all_temp) == 0:
                    fnew_all_temp = np.array([fnew])
                else:
                    fnew_all_temp = np.append(fnew_all_temp, [fnew], axis=0)

            # COLUMN: INPUT WAVELENGTH, ROW: AVAILABLE TEMPERATURE
            # LOOP OVER THE INPUT WAVELENGTHS AND CALCULATE A NEW FLUX
            # IF LEN(LAMBDA_) ==1, SPECIAL CARE IS TAKEN
            for j in range(len(lambda_)):
                if fnew_all_temp.ndim < 2:
                    yin = fnew_all_temp
                else:
                    yin = fnew_all_temp[:, j]
                intpObj = intp.interp1d(xin, np.log10(yin))
                m = 10 ** intpObj(temp0)
                GFluxArr = np.append(GFluxArr, m)
        else:
            for band in mags2use:
                pband = eval('self.' + band + 'pband')

                if pband.isoWavelength() > 110000 and resample:
                    lam_arr_all, flux_arr_all = FT.resample_model(lam_arr_all, flux_arr_all,
                                                                  min(pband.wavelength),
                                                                  max(pband.wavelength),
                                                                  pband=pband)

                Mflux = self.rsr_flux(pband, lam_arr_all, flux_arr_all)
                yin = Mflux
                # Interpolation object
                intpObj = intp.interp1d(xin, np.log10(yin))
                # Interpolated flux points

                m = 10 ** intpObj(temp0)

                GFluxArr = np.append(GFluxArr, m)

        return GFluxArr

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
        Temp: (float or np.array) Temperature in Kelvin."""

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
        """To calculate the blackbody irradiance in units of cgs
            Flux (erg s-1 cm-2 A-1). This module assumes projected
            emission from a source, so solid angle = pi. All physical
            constants are in units to facilitate this. Wavelength must
            be in units of angstrom, microns, cm, or meters. They will
            be converted to cm to carry out the calculation.
            units: string indicating units of input wavelength
                   will be converted to cm
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
        # CHECK PARAMETERS AND NORMALIZATION
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

    def calcRJ_spectrum(self, xarr, yarr):

        slope = (yarr[1] - yarr[0]) / (xarr[1] - xarr[0])
        yint = yarr[0] - slope * xarr[0]

        return slope, yint

    def modifiedBB(self, lambda_, p0, su2ea1=1, bands=None,
                   units='angstrom', **kwargs):
        """calculates a modified blackbody function given a powerlaw
        in conjunction with sed.SEDTools.blackbody
        lam0 = wavelength at which emission peaks.
        """

        try:
            pwr = kwargs['beta']
        except:
            # pwr = p0[2]# PASSED AS ARGUMENT TO FIT INSTEAD
            pwr = p0[2]
        try:
            lam0 = kwargs['lam0']
        except:
            lam0 = self.W3pband.isoWavelength()

        x = lambda_.copy().astype('float32')
        x0 = float(lam0)
        if units == 'angstrom':
            x *= 1e-8
            x0 *= 1e-8
        elif units == 'microns':
            x *= 1e-4
            x0 *= 1e-4
        elif units == 'meters':
            x *= 1e-2
            x0 *= 1e-2
        elif units == 'cm':
            x = x
            x0 = x0
        else:
            raise ValueError('Unit was not recognized')

            #  temp0 = p0[0]*100.
        # p0[0] = p0[0]*100.
        # CHECK PARAMETERS AND NORMALIZATION
        try:
            su2ea = (p0[1] ** 2) * su2ea1
            # su2ea = (p0[2]**2)*su2ea1
        except:
            su2ea = su2ea1
        # p0 = np.array(p0[0],p0[
        p0 = np.array([p0[0]])
        bbflux = self.blackbody(x, p0, su2ea1=su2ea, bands=bands, units='cm')

        mod_arr = np.array([])
        if bands is not None:
            for band in bands:
                pband = eval('self.' + band + 'pband')
                wav = pband.isoWavelength() * 1e-8
                mod_arr = np.append(mod_arr, (x0 / wav))
        else:
            mod_arr = (x0 / x)

        # mod_arr = (mod_arr<1).choose(mod_arr,1) # FOR REALLY INEFFICIENT SCATTERING DUST
        mod_arr = mod_arr ** pwr
        flux = bbflux * mod_arr
        return flux

    def calcBBTemp(self, T0, lamArr, flxArr):  # lam1,lam2,Flx1,Flx2):
        """Assumes Flux in units of erg/s/cm2/A and returns temperature in kelvin.
        to be used in conjunction with Newton-Raphson algorithm
        """

        lA, fA = lamArr.copy(), flxArr.copy()
        l1, l2 = lA[0], lA[1]  # lam1,lam2
        Flux1, Flux2 = fA[0], fA[1]
        Ratio = (Flux1 / Flux2) * (l1 / l2) ** 5
        gc = _H * _CS / _KB
        func = np.exp(gc / (l2 * T0)) - Ratio * np.exp(gc / (l1 * T0)) - 1 + Ratio
        return func

    def calcModTemp(self, T0, lam0, lamArr, flxArr):
        """Calculates the temperature of a modified blackbody with a power law
        assuming the power index has been parametrized and the ratio of two fluxes
        are taken. The Flux is in erg/s/cm2/A and wavelength in cm. This is to be used
        in conjunction with a Newton-Raphson/Bisecting algorithm to find zeros for the
        temperature equation.
        """

        lam0 = lam0  # 0.00115608 # W3 [cm]
        lA, fA = lamArr.copy(), flxArr.copy()
        l1, l2 = lA[0], lA[1]
        Flux1, Flux2 = fA[0], fA[1]
        Ratio = (Flux1 / Flux2) * (l1 / l2) ** 5
        gc = _H * _CS / _KB
        func1 = np.log10(((np.exp(gc / (l2 * T0)) - 1) / (np.exp(gc / (l1 * T0)) - 1)) * Ratio ** (-1))
        func2 = (gc / (T0 * lam0)) * np.exp(gc / (lam0 * T0)) / (np.exp(gc / (lam0 * T0)) - 1)
        func = func1 / (np.log10(l2 / lam0)) + func2 - 5.
        return func

    def photosphere(self, p0, su2ea2, modelinfo,
                    wave=(2000, 1e7), gridpts=10000):
        """Calculates the photospheric emission line from grid models
            given relevant parameters and griddata and temperature array,
            and range of wavelength to calculate emission for. Usually this
            should not exceed the limits for the grid models. If limits are exceeded
            flux is extrapolated linearly in log-log space, following Rayleigh-Jeans
            law.

            Parameters:
            ----------

            p0: Array of parameters to be passed to SEDTools.calc_grids.
            su2ea2: Normalization to be passed to SEDTools.calc_grids.
            modelinfo: Multi-dimensional array to be passed to sed.SEDTools.calc_grids
            wave: tuple or array of 2 elements; min and max of wavelength for flux
                     to be calculated. Be consistent with units. (angstroms)
            gridPts: Integer value for how many grid points you want returned
                      (aka resolution)

            Return:
            -------
            [xphot, fluxphot, slope, yint]
            xphot: wavelength array
            fluxphot: photospheric flux sampled/extrapolated at xphot
            slope: slope of extrapolated Rayleigh-Jeans line
            yint: y-intercept of Rayleigh-Jeans line obviously,
                  the last two are in log space
            """

        tempStar = p0[0] * 1000.
        try:
            su2eaRJ = (p0[1] ** 2) * su2ea2
        except:
            su2eaRJ = su2ea2

        gdat = MegaGrid[modelinfo]
        tempArr = gdat[-1]
        lam_arr_all, flux_arr_all = gdat[0], gdat[1]
        griddata = (lam_arr_all, flux_arr_all)
        # DETERMINE PHOTOSPHERIC CALCULATION CUT OFF WAVELENGTH BEFORE
        #  LINEAR INTERPOLATIN BEGINS
        xlim = np.max(lam_arr_all, axis=1).min()
        lamcut = xlim

        # CREATE SAMPLING POINTS FOR GRID
        xphot = np.logspace(np.log10(wave[0]), np.log10(wave[1]), gridpts)  # angstroms

        # This is only for Kurucz: xphot[angstrom], yphot[erg s-1 cm-2 A-1]
        # CHECK TO SEE IF EXTRAPOLATION IS NECESSARY GIVEN INPUT LIMITS
        if wave[1] > xlim:
            # CUT UP X ARRAY -- grid part and extroplation part

            ind_cut1, ind_cut2 = np.where(xphot < lamcut)[0], np.where(xphot >= lamcut)[0]
            xphot1, xphot2 = xphot[ind_cut1], xphot[ind_cut2]
            # yphot1 USES GRIDS, WHILE yphot2 USES RAYLEIGH-JEANS
            yphot1 = self.calc_grids(xphot1, p0, su2ea2, griddata, tempArr)
            # USE RAYLEIGH-JEANS FOR EXTRAPOLATION
            yphot2 = (np.pi * 1.4) * su2eaRJ * _CS * _KB * tempStar / (xphot2 * _ANG2CM) ** 3
            yphot2 = yphot2 / xphot2
            # PUT IT ALL TOGETHER
            xphot = np.concatenate([xphot1, xphot2])
            yphot = np.concatenate([yphot1, yphot2])


        else:
            # slope,yint=0,0
            yphot = self.calc_grids(xphot, p0, su2ea2, griddata, tempArr)
        # #  for Kurucz, yphot is in erg s-1 cm-2, xphot is in micron
        # return [xphot, yphot, slope, yint]
        return xphot, yphot

    def scaleSED2bands(self, scbdlist, usebdlist, yphot,
                       fluxm, fluxme, synflux):

        """

        Parameters
        ----------
        scbdlist
        usebdlist
        yphot
        fluxm
        fluxme
        synflux

        Returns
        -------

        """

        nirflux, nirbands = [], []
        yphot_unsc = yphot
        norm_wise_nir = 1.
        for band in scbdlist:
            if np.any(band == np.array(usebdlist)):
                nirbands.append(band)

        nirbands = np.array(nirbands)

        if len(nirbands) != 0:
            RJ_On = True
            flux_nirbands = AT.dict2list(fluxm, nirbands, '_flux')
            eflux_nirbands = AT.dict2list(fluxme, nirbands, '_flux')

            for band in nirbands:
                nirflux.append(synflux[band])
            nirflux = np.array(nirflux)
            wts = (flux_nirbands / eflux_nirbands) **2

            norm_wise_nir = np.average(flux_nirbands / nirflux, weights=wts)
            # CHANGE NORMALIZATION FOR SYNTHETIC FLUXES
            for key in synflux:
                synflux[key] = synflux[key] * norm_wise_nir

            yphot *= norm_wise_nir
        else:
            RJ_On = False
            print 'Unable to scale SED to W1 and W2'

        return norm_wise_nir, yphot, yphot_unsc, RJ_On


    def fit_photosphere(self, xlam, yfluxdat, p0, su2ea2,
                        modinfo, magfit, func):
        """

        Parameters
        ----------
        xlam
        yfluxdat
        p0
        su2ea2
        modinfo
        magfit
        func

        Returns
        -------

        """

        gdat = MegaGrid[modinfo]
        mg4phot = magfit['photmags']
        mg4scale = magfit['scalemags']
        yflux, yfluxerr = yfluxdat

        tempArr = gdat[-1]
        lam_arr_all, flux_arr_all = gdat[0], gdat[1]
        griddata = (lam_arr_all, flux_arr_all)

        Flx2Fit = AT.dict2list(yflux, mg4phot, '_flux')
        Flx2Fiterr = AT.dict2list(yfluxerr, mg4phot, '_flux')
        lam2Fit = AT.dict2list(xlam, mg4phot)
        Flx2scale = AT.dict2list(yflux, mg4scale, '_flux')
        Flx2scaleErr = AT.dict2list(yfluxerr, mg4scale, '_flux')

        FluxSED = func(lam2Fit, p0, 1, griddata, tempArr, mg4scale)

        print 'Bands used to fit photosphere: %s'%np.str(mg4phot)
        print 'Bands used to scale photosphere: %s'%np.str(mg4scale)

        # SURFACE TO OBSERVED FLUX NORMALIZATION WEIGHTED FROM ERRORS
        FluxNorm = np.average(Flx2scale / FluxSED, weights=1. / np.array(Flx2scaleErr))
        Rad = np.sqrt(FluxNorm / su2ea2)
        p0 = np.array([p0[0], Rad])  # radius squared

        # DEFINE INPUT PARAMETERS IN ORDER TO FIT TEMPERATURE
        self.fa = {'x': lam2Fit, 'y': Flx2Fit, 'err': Flx2Fiterr,
                   'func': func, 'griddata': griddata,
                   'tempArr': tempArr, 'mag2use': mg4phot,
                   'su2ea1': su2ea2}

        parinfo = [{'value': 0., 'relstep': 0, 'limits': [0, 0], 'limited': [0, 0]}
                   for l in range(2)]
        for j in range(2): parinfo[j]['value'] = p0[j]
        parinfo[0]['relstep'] = 0.1
        parinfo[1]['relstep'] = 0.3
        parinfo[0]['limited'] = [1, 1]
        parinfo[0]['limits'] = [tempArr[0], tempArr[-1]]

        mf = mt.mpfit(FT.deviates_from_model, parinfo=parinfo,
                      functkw=self.fa, quiet=1, maxiter=200000,
                      xtol=1e-16, ftol=1e-16, gtol=1e-16)

        p0, errors = mf.params, mf.perror
        try:
            self.chi2 = mf.fnorm / mf.dof
        except ZeroDivisionError:
            self.chi2 = -1
            print 'Degrees of freedom = 0'
        print 'chi2 = %.2f' % self.chi2

        radius = p0[1]
        tempnew = p0[0] * 1000.

        print 'Fitted Stellar Radius: %.3f Rsun' %radius
        print 'Fitted Stellar Temperature: %i K' %tempnew

        return radius,tempnew,mf

    def calc_temp(self, y, yerr, x, temparr, kw_cor):
        """

        Parameters
        ----------
        y : array of fluxes
        yerr: array of uncertaintiy in fluxes
        x: string array denoting bands
        temparr: arrays of temperatures
        kw_cor: dictionary of stuff from load_wfcorrection

        Returns
        -------
        res: residuals
        FluxNormed: flux normalized to blackbody fluxes
        alpha : array of chi2 calculations.

        """
        fluxNormed = []
        kcorList = []

        #  x should be list of strings with band info
        for bd in x:
            fluxNormed.append(np.array(kw_cor['F_cor_%s' % bd]))
            kcorList.append(np.array(kw_cor['kcor_%s' % bd]))

        kcorList = np.array(kcorList).transpose()
        Predicted = np.array(fluxNormed)
        a, b = Predicted[0], Predicted[1]
        x0, y0 = y[0], y[1]
        sigx, sigy = yerr[0], yerr[1]
        alpha = ((a * x0 / sigx ** 2) + (b * y0 / sigy ** 2)) / ((a / sigx) ** 2 + (b / sigy) ** 2)

        FluxNormed = (alpha * (Predicted.transpose() * kcorList).transpose()).transpose()
        res = np.subtract(FluxNormed, y)
        return res, FluxNormed, alpha

#  DICTIONARY THAT CONTAINS ALL THE DATA IN THE INPUT JSON FILE.


class DataLogistics:
    """ Set of tools to help with the logistical aspect
    of identifying the photospheric fit to stellar photometry:
    - loading data
    """

    def __init__(self, specs=None):

        # self.W1_lim, self.W2_lim = 4.5, 2.8
        # self.W3_lim, self.W4_lim = 3.5, -0.4
        #
        # self.J_lim = -1000
        # self.H_lim, self.Ks_lim = -1000, -1000
        # self.B_lim, self.V_lim = -1000, -1000
        import pdb
        pdb.set_trace()
        workingdir = directories.WorkingDir(specs['files']['stinfo_topdir'])
        
        starfile = opj(workingdir, specs['files']['stinfo_file'])

        empfile = opj(workingdir,specs['files']['stcolor_dir'],
                               specs['files']['bv_colorfile'])
        self.loadAllStars(starfile, specs['changekeys'])
        self.loadAllModels()
        self.loadEmpiricalData(empfile)



    def loadAllStars(self, starfile, changekeys):


        """Loads all data into a dictionary from starfile.
        If changekeys is active, it will try to replace the
        saturated corrected photometry keywords to standard
        WISE keywords.
        
        Unsaturated WISE photometry can be corrected using 
        procedures found in Patel,+2014.
        
        Need to add in part about self-correcting
        the photometry instead of relying on input.
        """
        if len(StarsDat) == 0:

            dat = np.genfromtxt(starfile, names=True, dtype=None)
            dat = dat.flatten()
            colnames = list(dat.dtype.names)
            for name in colnames:
                StarsDat[name] = dat[name]


            if changekeys:
                if 'W1mC' in colnames: #np.any(colnames == 'W1mC'):
                    StarsDat['W1m'] = StarsDat.pop('W1mC')
                    StarsDat['W1me'] = StarsDat.pop('W1meC')
                if 'W1mC' in colnames: #np.any(colnames == 'W2mC'):
                    StarsDat['W2m'] = StarsDat.pop('W2mC')
                    StarsDat['W2me'] = StarsDat.pop('W2meC')

        else:
            print 'Star"s data already loaded'


    def loadEmpiricalData(self,filename):
        """

        Parameters
        ----------
        filename

        Returns
        -------

        """

        if len(EmpDat) == 0: #is not None:
            dfemp = ascii.read(filename,comment='#')
            EmpDat['dat'] = dfemp
            #dat = EmpDat['test']


    def loadAllModels(self):

        """

        Returns
        -------

        """

        allg = np.unique(StarsDat['grav'])
        allmet = np.unique(StarsDat['met'])
        allmod = np.unique(StarsDat['model'])

        print '-------------------------------'
        print '      Loading All Gridmodels   '

        for mod in allmod:
            for z in allmet:
                for g in allg:
                    if (mod, g, z) not in MegaGrid:
                        MegaGrid[(mod, g, z)] = getattr(GMod, 'get_{}Grids'.format(mod))(g, z)
                        print 'Loaded %s of g=%s, met=%s'%(mod,g,z)


        print '       Done Loading Models     '
        print '-------------------------------'


class GridModels:

    @staticmethod
    def get_NextGenGrids(grav=None, met=0, model='NextGen', ext='.txt'):

        """
        Parameters:
        ------------
        grav, met: (strings) gravity and metallity values associated with files
                   for given model --> must match format for the models as indicated
                   in their respective readme files.
        model: (string) Model name. Should be same name as the directory all
               grids are located.
        ext: (string) Grid file extensions.

        Returns:
        ------------
        (lam_arr_all, flux_arr_all, temparr): tuple
        lam_arr_all: (1X N x M)-D numpy array with N wavelength points
                     for each M temperature in the specified grav and
                     met range.

        flux_arr_all: (1xNxM)-D numpy array with N flux points from
                      grid model to each wavelength point in lam_arr_all
                      for each M temperature for the specified met and grav.
        temparr: (1xM) numpy array of sorted temperatures from grid models
                  with values divided by 1000.
        [wavelength]:  angstroms
        [flux] :       ergs cm^{-2} s^{-1} A^{-1}
        [temperature]: Kelvin/1000.
        The grids are unevenly spaced, so filler 'nan's are placed
        in both the wavelength and flux points to make things run smoothly.

         """
        # CHECK FORMAT OF GRAV AND METALLICITY

        gdir = opj(intpdir,'Models', model)

        if grav is None:  # COLLECT ALL FILES OF ANY GRAV -- MEANT TO BE USED TO KEEP GRAV AS FREE PARAMETERS
            filesGrid = glob.glob(opj(gdir, 'lteNextGen*_%.1f%s' % (met, ext)))
        else:
            grav = grav / 10.
            filesGrid = glob.glob(opj(gdir, 'lteNextGen*_%.1f_%.1f%s' % (grav, met, ext)))

        if len(filesGrid) < 1:
            raise ValueError('There were no files matching your criteria. Try again.')

        ddat = {}
        tempArr = np.array([])
        gravArr = np.array([])

        for f in filesGrid:
            lam, flux = np.loadtxt(f,unpack=True)
            gridInfo = string.split(os.path.basename(f), '_')
            # print gridInfo
            ti = float(gridInfo[1])
            dat = np.array([lam, flux])
            ddat[str(ti / 1000.)] = dat

            tempArr = np.append(tempArr, ti)

        tempArr.sort()
        tempArr = tempArr / 1000.

        # CREATE 3D MATRIX: INDEXED BY [TEMP][WAVELENGTH][FLUX]
        gridnp = np.array([])
        i = 0
        si = len(tempArr)
        while i < si:
            # START OFF BY CREATING THE FIRST STACK and JUMP TWO INDICES
            gridi = ddat[str(tempArr[i])]

            if len(gridnp) == 0:
                try:
                    grid_i1 = ddat[str(tempArr[i + 1])]
                    gridnp = np.vstack(([gridi], [grid_i1]))

                except ValueError:
                    diff = len(gridi[0]) - len(ddat[str(tempArr[i + 1])])
                    zr = np.zeros(abs(diff)).astype(int)
                    if diff < 0:
                        gridi = np.insert(gridi, zr, zr, axis=1)
                        gridnp = np.vstack(([gridi], [grid_i1]))
                    elif diff > 0:
                        grid_i1 = np.insert(grid_i1, zr, zr, axis=1)
                        gridnp = np.vstack(([gridi], [grid_i1]))
                    else:
                        print 'Grids are not matched between %d, %d' % (tempArr[i], tempArr[i + 1])

                i += 2
            else:

                try:
                    gridnp = np.vstack((gridnp, [gridi]))
                except ValueError:
                    diff = len(gridnp[0][0]) - len(gridi[0])
                    zr = np.zeros(abs(diff)).astype(int)
                    if diff < 0:
                        gridnp = np.insert(gridnp, zr, zr, axis=2)
                        gridnp = np.vstack((gridnp, [gridi]))
                    elif diff > 0:
                        gridi = np.insert(gridi, zr, zr, axis=1)
                        gridnp = np.vstack((gridnp, [gridi]))
                    else:
                        print 'Grids are not matched between %d, %d' % (tempArr[i], tempArr[i - 1])

                i += 1
        lam_arr_all = gridnp[:, 0]
        flux_arr_all = gridnp[:, 1]

        return (lam_arr_all, flux_arr_all, tempArr)

    @staticmethod
    def get_ATLAS9Grids(grav, met, model='ATLAS9', ext='.fits'):
        """
        Parameters:
        ------------
        grav, met: (strings) gravity and metallity values associated with files
                   for given model --> must match format for the models as indicated
                   in their respective readme files.
        model: (string) Model name. Should be same name as the directory all
               grids are located.
        ext: (string) Grid file extensions.

        Returns:
        ------------
        (lam_arr_all, flux_arr_all, temparr): tuple
        lam_arr_all: (1X N x M)-D numpy array with N wavelength points
                     for each M temperature in the specified grav and
                     met range.

        flux_arr_all: (1xNxM)-D numpy array with N flux points from
                      grid model to each wavelength point in lam_arr_all
                      for each M temperature for the specified met and grav.
        temparr: (1xM) numpy array of sorted temperatures from grid models
                  with values divided by 1000.
        [wavelength]:  angstroms
        [flux] :       ergs cm^{-2} s^{-1} A^{-1}
        [temperature]: Kelvin/1000.
         """
        # CHECK TO MAKE SURE GRAV AND MET ARE IN CORRECT FORMAT
        # grav: g# # , met: p# # 
        grav, met = str(grav), str(met)
        if grav.find('g') == -1:
            grav = 'g' + grav
        if met.find('p') == -1:
            met = 'p0' + met

        # CHANGE DIRECTORY TO NEEDED METALLICITY FILE
        dir = opj(intpdir,'Models', model, 'k' + met)
        # THIS SELECTS OUT ONLY THE FILES THAT MEET THE METALLICITY
        # CRITERIA
        filesGrid = glob.glob(opj(dir, 'k' + met + '*.fits'))
        if len(filesGrid) < 1:
            raise ValueError('There were no files matching your criteria. Try again.')

        ddat = {}
        tempArr = np.array([])

        for f in filesGrid:
            try:
                hdui = fits.open(f)
                ti = hdui[0].header['TEFF']
            except:
                print "Something's wrong with %s" % f
            tbdata = hdui[1].data
            colnames = np.array(hdui[1].columns.names)
            # CHECK TO SEE IF GRAVITY SELECTED IS IN FITS FILE
            isGrav = np.where(grav == colnames)[0]
            if len(isGrav) != 0:
                tempArr = np.append(tempArr, ti)
                lam, flux = tbdata['WAVELENGTH'], tbdata[grav]
                lam, flux = np.array(lam), np.array(flux)
            else:
                print "Gravity %s was not found in file %s" % (grav, f)
            dat = np.array([lam, flux])
            ddat[str(ti / 1000.)] = dat
            # ddat[str(ti)] = dat
        tempArr.sort()
        tempArr = tempArr / 1000.
        gridnp = np.array([])

        # CREATE 3D MATRIX: INDEXED BY [TEMP][WAVELEN][FLUX]
        i = 0
        # for te in tempsArr:
        si = len(tempArr)
        while i < si:
            # START OFF BY CREATING THE FIRST STACK and JUMP TWO INDICES
            if len(gridnp) == 0:
                gridnp = np.vstack(([ddat[str(tempArr[i])]],
                                    [ddat[str(tempArr[i + 1])]]))
                i += 2
            else:
                gridnp = np.vstack((gridnp, [ddat[str(tempArr[i])]]))
                i += 1
        lam_arr_all = gridnp[:, 0]
        flux_arr_all = gridnp[:, 1]

        return (lam_arr_all, flux_arr_all, tempArr)

    def treat_NextGenSpecModels(self, met=0, grav='all', ext='.spec'):
        """Use this in order to convert the raw files of the NextGen
           atmospheric models (i.e. ".spec" files) from France Allard's page.
           The raw files are listed as such:
           low-res: low-resolution (100-25,000A @ 2A) spectra
           (directly from the model iterations)
           Obtained from: ftp://ftp.hs.uni-hamburg.de/pub/outgoing/phoenix
           Spectrum files:
           Wavelengths are in vacuum, fluxes and BB(Teff) are in
           cgs (erg/s/cm^2/cm) (yes, that is correct. NOT in
           erg/s/cm^2/A). Also note that the fluxes and BB's are given
           directly, NOT as log10 as in previous versions of
           the model grid.

           line1: Teff logg [M/H] of the model
           Wavelength, flux and bbflux are stacked. This is not a column wise file.
           This module assumes all files have been unzipped. It does not discriminate
           with respect to temperature, although selections can be made to transform grids
           of a given solar metallicity and gravity.

           Input: met: string or float. "All" will choose files for all metallicities. Input of a
                       number from the models will select files of that metallicity.
                  grav: Same as met
                  ext: What the original extension of the raw files are. I suggest to keep them as .'spec'


           Returns: Individual ascii files of the grids, two column of wavelength and F_lam as in the
                    raw files. File names: 'lte_temp_grav_met.txt'
                    Flam units: erg/s/cm^2/Angstrom
                    Wavelength: Angstrom


        ********** WARNING *******************
        THERE EXIST POINTS IN THE GRIDS THAT CORRESPOND TO THE SAME DATA POINT IN BOTH LAM AND
        FLUX SPACE -- PROCEED WITH CAUTION IF NAN'S OCCUR
        **************************************
        """

        conv2ang = 1e8

        newdir = opj(intpdir, 'NextGen2')
        # os.mkdir(newdir)
        dir = '~/Desktop/PHOENIX'
        if (met == 'all') and (grav == 'all'):
            filesGrid = glob.glob(opj(dir, 'lte*.NextGen%s' % ext))
        elif (met == 'all') and (grav != 'all'):
            filesGrid = glob.glob(opj(dir, 'lte*-%.1f-*.NextGen%s' % (grav, ext)))
        elif (met != 'all') and (grav == 'all'):
            filesGrid = glob.glob(opj(dir, 'lte*-%.1f.NextGen%s' % (met, ext)))
        elif (met != 'all') and (grav != 'all'):
            filesGrid = glob.glob(opj(dir, 'lte*-%.1f-%.1f.NextGen%s' % (grav, met, ext)))
        else:
            pass

        if len(filesGrid) < 1:
            raise ValueError('There were no files matching your criteria. Try again.')

        else:
            for file in filesGrid:
                # f = open(file,'r').readlines()
                f = open(file, 'r')
                fread = f.read()
                elements = fread.strip()
                elements2 = elements.split(' ')

                fstrip = map(string.strip, elements2)
                datarr = np.array(fstrip)
                ind = np.where(datarr == '')[0]
                dat = np.delete(datarr, ind)

                temp, grav, met, nwave = float(dat[0]), float(dat[1]), \
                                         float(dat[2]), int(dat[6])

                if grav <= 5.5:
                    dat = np.delete(dat, np.arange(0, 11)).astype('float32')
                    wave = dat[np.arange(nwave)]
                    flux = dat[np.arange(nwave, nwave * 2)]
                    # indkeep = np.where(wave<=30000.)[0]
                    # wavKeep, fluxKeep = wave[indkeep],flux[indkeep]
                    wavKeep, fluxKeep = wave, flux
                    # units of erg/s/cm^2/Angstrom
                    DataOut = np.column_stack((wavKeep, ((fluxKeep) / conv2ang)))
                    fileout = opj(newdir, 'lteNextGen_%.1f_%.1f_%.1f.txt' % (temp, grav, met))
                    np.savetxt(fileout, DataOut)
                else:
                    pass


GMod = GridModels()


class Bandpass:
    def __init__(self, fileName, normalise=True, inputUnits='angstroms'):
        """This code loads a passband file with wavelength in the first
          column and the RSR transmission data for that particular
          bandpass. The first line should be in the form of
          #!NNNNN.NN[0,f]/ where NNNN.NN is the isophotal wavelength in
          angstroms for this particular bandpass. Even if no isophotal
          wavelength is listed still have the !#  part in there. The
          second line should have the column headings "wav" "trans".

          Converts the input file wavelength to angstroms. You can tell
          the code which units the file is in and it will convert to
          angstroms.
          Options are: nanometers, microns, mm, GHz inputted as a string.
        """
        self.file = fileName

        df = np.genfromtxt(self.file,skip_header=1,names=True)

        self.wavelength = df['wav']
        self.transmission = df['trans']

        if inputUnits == 'angstroms':
            pass
        elif inputUnits == 'nanometers':
            self.wavelength *= 10.0
        elif inputUnits == 'microns':
            self.wavelength *= 10000.0
        elif inputUnits == 'mm':
            self.wavelength *= 1e7
        elif inputUnits == 'GHz':
            self.wavelength = 3e8 / (self.wavelength * 1e9)
            self.wavelength *= 1e10
        else:
            raise Exception, "DAFUQ? Units no make sense"

        # Sort into ascending order of wavelength otherwise normalization will be wrong

        merged = np.array([self.wavelength, self.transmission]).transpose()
        sortedMerged = np.array(sorted(merged, key=operator.itemgetter(0)))
        self.wavelength = sortedMerged[:, 0]
        self.transmission = sortedMerged[:, 1]

        if normalise:
            self.transmission = self.transmission / np.trapz(self.transmission, self.wavelength)

        self.interpolator = intp.interp1d(self.wavelength, self.transmission, kind='linear')

        self.zmdata = self.metaPassDat(self.file)

    def metaPassDat(self,dfile):
        """
        This will extract the first line in the RSR data files that contains the zero point
        wavelengths, frequencies and fluxes, denoted by the string '#!'

        Parameters
        ----------
        dfile : string that points to the file name of the RSR bandpass.

        Returns
        -------
        list of zero point data. Check readme file in RSR file for more information.

        """

        f = open(dfile)
        first = f.readline()

        self.S = re.compile('^#\!|\[[0-9]+\]\/|\[ *\w *, *\w+ *\]\/')
        header = re.split(self.S,first)[1:-1]
        zmdata = map(string.strip,header)

        return zmdata

    def rescale(self, maxTransmission):
        """ Rescales passband so that maximum value of the transmission is equal to
            maxTransmission
            maxTransmission: float - max value to rescale transmission curve to.
        """
        self.transmission /= maxTransmission  # self.trans

    def pivotWavelength(self):
        """Calculates pivot wavelength for the passband. This is the same as equation (3) of
        Carter et al. 2009, and equation A11 in Tokunaga & Vacca 2005.
        """

        a = np.trapz(self.transmission * self.wavelength, self.wavelength)
        b = np.trapz(self.transmission / self.wavelength, self.wavelength)
        pivWavelength = np.sqrt(a / b)

        return pivWavelength

    def isoWavelength(self):
        """Extracts the isophotal wavelength from the bandpass file if it's
        there. Otherwise it uses the calculated effective wavelength of the
        bandpass. The isophotal wavelength is kept under # !NNNNN.NN[0,f]/ in
        the first line of the bandpass file.
        """
        try:
            isowavelength = float(self.zmdata[0]) # h.header[0]

        except IndexError:
            isowavelength = self.pivotWavelength()
            print 'Error in extracting isowavelength. Using pivotWavelength instead for at %.5f' %isowavelength

        return isowavelength

    def fluxVegaZeroPointLam(self):
        """Extracts the zero point fluxes and related uncertainties from the bandpass
        if it's there. Otherwise it throws and error and lets you know you are now
        attached to another object by an inclined plane, wrapped helically around an axis.
        The information is kept in the header file that looks like # !N1[0,f]/ N2[1,f]/ N3[2,f]/
        where N1 is the isophotal wavelength and the N2, N3 are the zero point fluxes
        respectively.
        """

        try:
            zp = self.zmdata[1]
            zperr = self.zmdata[2]
        except IndexError:
            print 'No zero point flux available. Check RSR file %s' % self.file

        zpflx, zpflxErr = float(zp), float(zperr)

        return zpflx, zpflxErr

    def isoFrequency(self):

        try:
            isofrequency = float(self.zmdata[3])
        except IndexError:
            print 'No other option for frequency available. Place value in header file'

        return isofrequency

    def fluxVegaZeroPointFreq(self):


        try:
            zp = self.zmdata[4]
            zperr = self.zmdata[5]
        except IndexError:
            print 'No zero point flux available. Check RSR file %s' % self.file

        zpflx, zpflxErr = float(zp), float(zperr)

        return zpflx, zpflxErr


class Flatbandpass:
    def __init__(self, waverange=(None, None), cntr=None, normalize=True, inputUnits='angstroms'):
        """Creates a passband object with a flat response curve
        
        """

        if waverange[0] is None or waverange[1] is None:
            raise ValueError('No wavelength range defined for the flat spectrum')
        else:
            self.isoWavelength = cntr
            self.pivotWavelength = cntr
            self.wavelength = np.linspace(waverange[0], waverange[1], 200)
            self.transmission = np.ones(len(self.wavelength))

            if inputUnits == 'angstroms':
                pass
            elif inputUnits == 'nanometers':
                self.wavelength *= 10.0
                self.pivotWavelength *= 10.0
            elif inputUnits == 'microns':
                self.wavelength *= 10000.0
                self.pivotWavelength *= 10000.0
            elif inputUnits == 'mm':
                self.wavelength *= 1e7
                self.pivotWavelength *= 1.7
            elif inputUnits == 'GHz':
                self.wavelength = 3e8 / (self.wavelength * 1e9)
                self.wavelength *= 1e10
                self.pivotWavelength = 3e8 / (self.pivotWavelength * 1e9)
                self.pivotWavelength *= 1e10
            else:
                raise Exception, "DAFUQ? Units no make sense"

            self.isoFrequency = _CS / self.isoWavelength
