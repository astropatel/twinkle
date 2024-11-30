"""
twinkle.py by Rahul I. Patel (ri.patel272@gmail.com)

The class Star can be instantiated for a single star, giving you a STAR object
for each star that you want to characterize.

"""

import os
import copy
import logging
import numpy as np
from . import sed
# from twinkle import sed
from astropy import constants as con

logging.basicConfig(filename='example.log', filemode='w', level=logging.DEBUG)

STools = sed.SEDTools()

__author__ = 'Rahul I. Patel <ri.patel272@gmail.com>, Joe Trollo'

#  set up constants
_CS = con.c.to('cm/s')
_H = con.h.to('erg s')
_KB = con.k_B.to('erg/K')
_RSUN = con.R_sun.to('cm')
_LSUN = con.L_sun.to('erg/s')

#  SET UP UNIT CONVERSION.
_PC2CM = 3.08568025e+18
_SOLRAD2CM = 69550000000.0
_AU2CM = 14959787070000.0
_ANG2MICRON = 0.0001
_MICRON2ANG = 1. / _ANG2MICRON
_ANG2CM = 1e-8

CONST_1 = (_SOLRAD2CM / _PC2CM) ** 2
CONST_2 = _AU2CM ** 2 / (4 * np.pi * _PC2CM ** 2)

Photometry_spCheckList = ['mags2use0', 'mags4Phot0', 'mags4scale0']

# PLOTTING LABELS

y_flux_label = r'$\lambda F_{\lambda}\ [erg\,\ s^{-1}\ cm^{-2}] $'
x_wav_label_microns = r'$\lambda\ [\mu m]$'

class Fuel:
    """
    This should be called before all others.
    Initialization process loads the input script file, uploads and saves the
    data of the input star file, empirical color file, photospheric grid models
    and is accessible through :class:`sed.DataLogistics`
    """

    def __init__(self, jfile):
        """
        Loads the user input parameter file data

        Args:
            jfile: (str) json file of inputs
        """
        #  Load data from starfile
        #  Load data from empirical color file

        if not os.path.isfile(jfile):
            raise ValueError(f'{jfile} is not a file.')
        # CHECK TO SEE IF DICTIONARY OF STELLAR DATA HAS BEEN POPULATED
        self.data_manager = sed.DataLogistics(jfile)
        logging.info(f'JSON parameterfile used:, {jfile}')


class Star:
    """
        Instantiating StarObject does a number of things, foremost of
        which is to create a StarObject for ONE particular star. The
        initialization process loads the input script file, uploads
        and saves the data in the input star file, empirical color
        file, photospheric grid models.

        In addition, it cleans up the photometry, converts photometry
        to flux, fits photometric flux to the best photospheric model,
        and possibly scales the new photosphere to a subset of the
        input photometry.

        FYI: I get bored and some of the print statements may be... odd.
    """

    def __init__(self, fuel, starname=None):
        """

        Instantiates a "twinkling star" object (funny, I know).

        Args:
            fuel (twinkle.Fuel class obj)
            starname (str):  Name of star in MainName column of stellar file
                        Providing both sid and starname will halt execution.

        """

        self.starsdat = fuel.data_manager.StarsDat
        self.MegaGrid = fuel.data_manager.MegaGrid
        self.emdat = fuel.data_manager.EmpDat # sed.EmpDat['dat']
        specs = fuel.data_manager.specs

        self.starname = starname

        #  Find star in file if name given
        if self.starname:
            try:
                self.star_dat_i = self.starsdat.loc[self.starname]
                logging.info(f'WORKING STAR:{self.starname}')
            except ValueError:
                print(f'{self.starname} was not found in the input file')

        else:
            raise ValueError('Provide either name of location index, not both.')

        #  Set saturation limits for various bands
        #  Might need to put this in JSON file.
        self.W1_lim, self.W2_lim = 4.5, 2.8
        self.W3_lim, self.W4_lim = 3.5, -0.4

        self.J_lim = -1000
        self.H_lim, self.Ks_lim = -1000, -1000
        self.B_lim, self.V_lim = -1000, -1000

        self.mags2use = []
        self.mags4Dust = []
        self.mags4Phot = []
        self.mags4scale = []

        self.vegaMagDict = {}
        self.vegaMagErrDict = {}
        self.photFlux = {}

        self.fullspectrum = None
        self.StarRadius, self.StarTemp = None, None
        self.mfit = None

        # ADD PHOTOMETRY FROM MAGNITUDE LIST IN JSON FILES
        # REMOVE SATURATED BANDS AND REPLACE NULL VALUES
        # vegaMagDict and errdict ARE FILLED UP HERE.

        self.cleanphotometry(specs)

        #  Gather up basic stellar info from file
        self.disti = 1000. / self.star_dat_i['plx']
        self.spti = self.star_dat_i['spt']
        self.met = self.star_dat_i['met']
        self.modeli = self.star_dat_i['model']

        self.su2ea = CONST_1 / self.disti ** 2
        self.su2ea_dust = CONST_2 / self.disti ** 2

        #  Obtain grav and met from file or guess based on B-V
        self.g = self.star_dat_i['grav']
        self.T0 = self.star_dat_i['temp'] * 1e-3
        bv = self.emdat['B-V']

        try:
            bvi = self.vegaMagDict['BJ'] - self.vegaMagDict['VJ']

        except KeyError:
            bvi = self.bv_unusedDict['BJ'] - self.bv_unusedDict['VJ']

        indebv = np.searchsorted(bv, bvi)
        if self.g is None:
            self.g = self.emdat['log(g)_x10'][indebv]
        if self.T0 is None:
            self.T0 = self.emdat['Tinit'][indebv] / 1000.

        # Convert Photometry to Flux
        self.fluxTup = STools.batch_mag2fluxZPLam(self.mags2use,
                                                  self.vegaMagDict,
                                                  self.vegaMagErrDict)

        self.flux, self.fluxerr = self.fluxTup
        self.wave = STools.get_eff_wavelengths(self.mags2use)

        # FIT PHOTOSPHERE

        if specs['fitphot']:
            self.fitPhotosphere(self.su2ea, self.T0,
                                self.modeli, self.g, self.met,
                                self.MegaGrid)

            params, perror = self.mfit.params, self.mfit.perror
            self.p_up, self.p_down = params + perror, params - perror

            logging.info(f'Photosphere for {self.starname} fit with '
                         f'T={self.StarTemp} K')

            wave_min, wave_max = specs['spec_sample']['wave_min'], \
                                 specs['spec_sample']['wave_max']
            gridpts = specs['spec_sample']['gridpts']
            modelinfo = (self.modeli, self.g, self.met)
            self.StarPhotosphere = STools.photosphere([self.StarTemp / 1000.,
                                                       self.StarRadius],
                                                      self.su2ea, modelinfo,
                                                      self.MegaGrid,
                                                      wave=(wave_min, wave_max),
                                                      gridpts=gridpts)

            StarPhot1 = STools.photosphere(self.p_up,
                                           self.su2ea, modelinfo,
                                           self.MegaGrid,
                                           wave=(wave_min, wave_max),
                                           gridpts=gridpts)

            StarPhot2 = STools.photosphere(self.p_down,
                                           self.su2ea, modelinfo,
                                           self.MegaGrid,
                                           wave=(wave_min, wave_max),
                                           gridpts=gridpts)

            flxes = np.array([StarPhot1[1], StarPhot2[1]])
            maxind = np.argmax(np.sum(flxes, axis=1))
            minind = np.argmin(np.sum(flxes, axis=1))

            self.StarPhotosphere_up = (StarPhot1[0], flxes[maxind])
            self.StarPhotosphere_down = (StarPhot1[0], flxes[minind])

            self.StarPhotosphere_unsc = self.StarPhotosphere

            for band in self.mags2use:
                flxt = STools.rsr_flux(getattr(STools, f'{band}pband'),
                                       *self.StarPhotosphere)[0]
                self.photFlux[f'{band}'] = flxt

        else:

            logging.info(f'Photosphere fit not requested for {self.starname}')
            if specs['scalephot']:
                new_phot = STools.scaleSED2bands(specs['phot']['scaleSEDbands'],
                                                 self.mags2use,
                                                 self.StarPhotosphere[1],
                                                 self.flux, self.fluxerr,
                                                 self.photFlux)

                phot_norm_fac, yphot, yphot_unsc, RJ_On = new_phot
                self.StarPhotosphere[1] = yphot

    def writeSED(self, filename='sed.txt',
                 comment='# lambda: Angstroms, f_lambda: erg/s/cm^2/Angstrom, Teff={:.1f}K,rad={:.3f}Rsol.\n'):
        """
        Function to write out fitted SED to file.

        Args:
            filename: (str) file of where SED will be saved.
            comment: (str) comment starting with '#' of any comment to go on
                first line of the saved SED file.

        """
        comment = comment.format(self.StarTemp, self.StarRadius)

        if not self.StarPhotosphere:
            raise ValueError('Photosphere was not created --> cant be written out.')

        with open(filename, 'w') as file:
            file.write(f'{comment}')
            file.write('lambda\t f_lambda\n')
            np.savetxt(file, np.transpose(self.StarPhotosphere), delimiter='\t')

        logging.info(f'SED for {self.starname} saved to {filename}')

    def cleanphotometry(self, specs):
        """
        The purpose of this module is to remove any photometry that is
        null and supplement any null error measurements with 5% of the
        photometric value. It also removes saturated bands for that particular
        star.

        arrays created:
        1) mags2use = all the bands that are used.
        2) mags4Phot = bands that are used to fit to model photosphere.
        3) mags4scale = bands used to scale raw model to input flux.
        4) mags4Dust = band used to calculate excess.

        Also it removes optical bands if they're late type stars.

        Args:
            specs: (dict) dictionary created using json parameterfile.

        """

        vegaMagDict_temp, vegaMagErrDict_temp = {}, {}

        mags2use0 = copy.copy(specs['phot']['mags2use0_original'])
        mags4Phot0 = copy.copy(specs['phot']['mags4Phot0_original'])
        mags4scale0 = copy.copy(specs['phot']['mags4scale0_original'])
        mags4Dust0 = copy.copy(specs['phot']['mags4Dust0'])

        #  ========================================
        #  Whether to use W3 and W2 to fit photosphere
        if specs['W3Adapt']:
            mags4Phot0 = self.W3Adopt(specs, mags4Phot0, True)
        if specs['W2Adapt']:
            mags4Phot0 = self.W2Adopt(specs, mags4Phot0, True)
        # Remove optical bands cause of late spectral type?
        # todo - this piece doesn't work
        if self.star_dat_i['NoOptical']:

            self.bv_unusedDict = {}

            for arr in Photometry_spCheckList:
                # arr = np.array(arr)
                for mv in specs['phot']['Remove_RedStars']:
                    try:
                        # todo: change string formatting from %s
                        ind = np.where(np.array(eval(arr)) == mv)[0]
                        self.bv_unusedDict[mv] = self.star_dat_i[mv + 'm']
                        exec(f'{arr} = np.delete({arr},ind)')
                    except ValueError:
                        pass
        # ====================================================
        # KEEP ONLY VALID PHOTOMETRIC MEASUREMENTS.
        for mv in mags2use0:

            try:
                temp_mf = self.star_dat_i['{}m'.format(mv)]
            except KeyError as e:
                print(f"Error: The key '{e.args[0]}' was not "
                      f"found in the mags2use0 list.")
                temp_mf = self.star_dat_i['{}'.format(mv)]
            # CHECK IF IT'S NULL. REMOVE IF SO FROM ENTIRE LIST.
            if temp_mf == 'null':
                try:
                    mags2use0.remove(mv)
                    logging.info(f'{mv} band removed from mags2use0')
                except ValueError:
                    logging.error(f'Error in removing {mv} from mags2use')
                try:
                    mags4Phot0.remove(mv)
                    logging.info(f'{mv} band removed from mags4Phot0')
                except ValueError:
                    logging.error(f'Error in removing {mv} from mags4phot')
                try:
                    mags4scale0.remove(mv)
                    logging.info(f'{mv} band removed from mags4scale0')
                except ValueError:
                    logging.error(f'Error in removing {mv} from mags4scale')
                try:
                    mags4Dust0.remove(mv)
                    logging.info(f'{mv} band removed from mags4Dust0')
                except ValueError:
                    logging.error(f'Error in removing {mv} from mags4dust')

            # OTHERWISE ADD IT.
            else:
                if '_flux' in mv:
                    fj = temp_mf
                    mv = mv.strip('_flux')
                    lam = eval('STools.{}pband.isoWavelength()'.format(mv))
                    try:
                        freq = eval('STools.{}pband.isoFrequency()'.format(mv))
                    except UnboundLocalError:
                        print(('Iso Freq. for {} will be calculated using'
                              ' iso wavelength.'.format(mv)))
                        freq = None

                    # WHEN THERE IS A NULL VALUE
                    try:
                        efj = float(self.star_dat_i[f'{mv}_fluxe'])
                    except ValueError:
                        efj = 0.05 * fj

                    # DIVIDE BY LAMBDA TO OBTAIN ERG/S/CM2/ANG
                    fcgs, efcgs = np.array(STools.Jy2cgs((fj, efj), (freq, 0), (lam, 0))) / lam

                    fmag, efmag = STools.fluxZPLam2mag(eval('STools.{}pband'.format(mv)),
                                                       fcgs, efcgs)

                    vegaMagDict_temp[mv] = fmag
                    vegaMagErrDict_temp[mv] = efmag
                    print(f'{mv} was converted to magnitude')
                    # Assumes the rest is a magnitude
                else:
                    vegaMagDict_temp[mv] = temp_mf
                    if self.star_dat_i[f'{mv}me'] == 'null':
                        vegaMagErrDict_temp[mv] = 0.05 * (vegaMagDict_temp[mv])
                    else:
                        vegaMagErrDict_temp[mv] = float(self.star_dat_i[f'{mv}me'])

        # CHECK SATURATION LIMITS AND REMOVE FROM ALL LISTS
        # ========================================
        mags2use0 = [s.replace('_flux', '') for s in mags2use0]
        mags4Dust0 = [s.replace('_flux', '') for s in mags4Dust0]
        mags4Phot0 = [s.replace('_flux', '') for s in mags4Phot0]
        mags4scale0 = [s.replace('_flux', '') for s in mags4scale0]

        if specs['satcheck']:
            self.mags2use = self.keep_unsatmags(vegaMagDict_temp, mags2use0)
            self.mags4Dust = self.keep_unsatmags(vegaMagDict_temp, mags4Dust0)
            self.mags4Phot = self.keep_unsatmags(vegaMagDict_temp, mags4Phot0)
            self.mags4scale = self.keep_unsatmags(vegaMagDict_temp, mags4scale0)
        else:
            self.mags2use = mags2use0
            self.mags4Dust = mags4Dust0
            self.mags4Phot = mags4Phot0
            self.mags4scale = mags4scale0

        self.mags2use = list(np.sort(self.mags2use))
        self.mags4Dust = list(np.sort(self.mags4Dust))
        self.mags4Phot = list(np.sort(self.mags4Phot))
        self.mags4scale = list(np.sort(self.mags4scale))

        # REMOVE NON-USED MAGNITUDES FROM DICTIONARY
        #  ========================================
        for mv in mags2use0:
            tmp = mv.strip('_flux')  # just in case _flux is used.
            self.vegaMagDict[tmp] = vegaMagDict_temp[tmp]
            self.vegaMagErrDict[tmp] = vegaMagErrDict_temp[tmp]

    def calc_excessflux(self):
        """
        Calculates the excess flux at the passbands listed in mags4Dust
        using the photospheric fit.

        flux and uncertainty for excesses are stored in :attr:`Star.fluxEx` and
        :attr:`Star.efluxEx`.
        Wavelengths are in :attr:`Star.WaveEx` in angstroms.

        """

        # CREATE NEW tauA DICTIONARY FOR EXCESS FLUX
        self.excessFlux = {}
        self.excessFlux_e = {}
        self.excessFlux_wave = {}

        # CALCULATE EXCESS FLUX AND ADD TO TAUA
        for band in self.mags4Dust:
            exflux = self.flux[band + '_flux'] - self.photFlux[band]
            print((band, ' excess flux ', exflux, 'erg/s/cm2/ang'))
            self.excessFlux[band + '_flux'] = exflux
            self.excessFlux_wave[band] = self.wave[band]
            self.excessFlux_e[band + '_flux'] = self.fluxerr[band + '_flux']

        self.fluxEx = np.array(list(zip(*sorted(self.excessFlux.items())))[1])
        magsorder, self.waveEx = np.array(list(zip(*sorted(self.excessFlux_wave.items()))))
        self.waveEx = self.waveEx.astype('float64')
        self.efluxEx = np.array(list(zip(*sorted(self.excessFlux_e.items())))[1])

    def keep_unsatmags(self, vegaDict, magsCheck):
        """
        This will remove any saturated photometry from the
        string lists (mags2use,mags4scale,mags4Phot,mags4Dust).

        Parameters
        ----------
        vegaDict: (dict) Photometry of bands for a particular star.
        magsCheck: (list/arr) subset of bands to check.

        Returns
        -------
        magsCheck1: (list/arr) subset of magsCheck that don't have
                    saturated photometry.
        """

        magsCheck1 = magsCheck
        if not isinstance(magsCheck1, list):
            magsCheck1 = list(magsCheck1)

        mags_temp = np.array(magsCheck).copy()
        for mv in mags_temp:
            try:
                if vegaDict[mv] < eval('self.%s_lim' % mv):
                    magsCheck1.remove(mv)
                    logging.info('%s band is saturated:%.4f<%.4f'
                                 % (mv, vegaDict[mv], eval('self.%s_lim' % mv)))
                else:
                    pass
            except:
                pass

        return magsCheck1

    def W3Adopt(self, specs, pmaglist=None, simple=True):
        """
        Determines whether the WISE W3 photometry
        should be included in the list of  photometric
        points to be used to pin down the photosphere.
        v1 uses the previously determined W1-W3 and W2-W3 cuts
        as well as the excess SNRs (Sigma_{E[Wi-W3]}) to
        determine whether the photometry is photospheric for
        a given star....
        Or whether the photometry is saturated or not.

        Parameters
        ----------
        specs: (dict) dictionary created using json file.
        pmaglist : (list) includes the string list of photometry
                        used to pin down photosphere.
        simple : (bool) if simple is set, it uses the saturation limit to
                 to determine whether W3 should be used.

        Returns
        -------
        plist: (list) new list of photospheric string values
        """

        if pmaglist is None:
            raise ValueError('No maglist specified for W3Adopt')
        else:
            plist = pmaglist

        if simple:
            if self.star_dat_i['W3m'] > 3.8:
                plist = np.append(plist, 'W3')
        else:
            W13_cut = specs['WISE_excess']['W13_cut']
            W23_cut = specs['WISE_excess']['W23_cut']

            w1w3snr = self.star_dat_i['W1W3SNR']
            w2w3snr = self.star_dat_i['W2W3SNR']

            try:
                w3snr1 = float(w1w3snr)
            except ValueError:
                w3snr1 = -100

            try:
                w3snr2 = float(w2w3snr)
            except ValueError:
                w3snr2 = -100

            if (w3snr1 is None) and (w3snr2 is None):
                pass

            elif w3snr1 < W13_cut and w3snr1 > -1 * W13_cut and w3snr1 != -100:
                if (w3snr2 < W23_cut and w3snr1 > -1 * W13_cut and w3snr2 != -100) or (w3snr2 == -100):
                    plist = np.append(plist, 'W3')

            elif w3snr2 < W23_cut and w3snr2 > -1 * W23_cut and w3snr2 != -100:
                if (w3snr1 < W13_cut and w3snr2 > -1 * W23_cut and w3snr1 != -100) or (w3snr1 == -100):
                    plist = np.append(plist, 'W3')
            elif w3snr1 < W13_cut and w3snr1 > -1 * W13_cut and w3snr1 != -100 and (
                    w3snr2 < W23_cut and w3snr2 > -1 * W23_cut and w3snr2 != -100):
                plist = np.append(plist, 'W3')
            else:
                pass

        return plist

    def W2Adopt(self, specs, pmaglist=None, simple=True):
        r"""Determines whether the WISE W2 photometry
        should be included in the list of  photometric
        points to be used to pin down the photosphere.
        v1 uses the previously determined W1-W2 cuts
        as well as the excess SNRs (Sigma_{E[W1-W2]}) to
        determine whether the photometry is photospheric.

        Args:
            specs (dict):  dictionary created using json file.
            pmaglist (list): includes the string list of photometry used to pin
                down photosphere.
            simple (bool): if simple is set, it uses the saturation limit to
                determine whether W3 should be used.

        Returns:
            plist (list): new list of photospheric string values
        """

        if pmaglist is None:
            raise ValueError('No maglist specified for W2Adopt')
        else:
            plist = pmaglist

        if simple:
            if self.star_dat_i['W2m'] > 3.8:
                plist = np.unique(np.append(plist, 'W2'))


        else:
            W12_cut = specs['WISE_excess']['W12_cut']
            w1w2snr = self.star_dat_i['W1W2SNR']
            try:
                w2snr = float(w1w2snr)
            except ValueError:
                w2snr = -100  # THIS VALUE IF THERE IS NO SNR VALUE LISTED

            if w2snr <= W12_cut:
                plist = np.append(plist, 'W2')

        return plist

    def resetFullSpectrum(self):
        """
        Need doc string
        Returns
        -------

        """
        self.fullspectrum = self.StarPhotosphere[1][:]

        return

    def fitPhotosphere(self, sconst, T0, model, grav, met, MegaGridModels):
        r"""
        Fit a stellar photosphere model to bands in the `mags4Phot` list/array.
        This uses `mpfit` from `mosaic_tools.py`.

        Results are saved in the following internal variables:
        - **StarRadius**: Stellar radius in Solar Radius units.
        - **StarTemp**: Temperature of the fitted star in Kelvin.

        Args:
            sconst (float): Scaling factor for raw fluxes. It equals:
                (Solar Radius [cm] /  1 pc [cm])^2 * (1/Stellar Distance [pc]^2)
            T0 (float):
                Initial guess for the temperature, in units of Kelvin/1000.
            model (str):
                Name of the model to use. Options are `"Atlas9"` or `"NextGen"`.
            grav (float):
                Surface gravity, expressed as :math:`10\log(g)\ [\mbox{cm}/s^2]`.
            met (float):
                Logarithmic metallicity. A value of `0` represents solar metallicity.
            MegaGridModels (dict): all atmospheric grid models with metallicity,
                grav, and temperature keys.
        Returns:
            None:
                This function does not return a value but updates internal variables:
                - **StarRadius**: Stellar radius in Solar Radius units.
                - **StarTemp**: Temperature of the fitted star in Kelvin.

        """

        modeltype = (model, grav, met)
        if np.array(self.mags4scale).size == 0:
            raise ValueError('mags4scale is empty.')
        if np.array(self.mags4Phot).size == 0:
            raise ValueError('mags4Phot is empty.')
        mfitlist = {'photmags': np.array(self.mags4Phot),
                    'scalemags': np.array(self.mags4scale)}

        rawfluxdat = (self.flux, self.fluxerr)

        fit_dat = STools.fit_photosphere(self.wave, rawfluxdat,
                                         [T0], sconst, modeltype,
                                         mfitlist, STools.calc_grids,
                                         MegaGridModels)

        self.StarRadius, self.StarTemp = fit_dat[0], fit_dat[1]
        self.mfit = fit_dat[2]

    def plot_photrange(self, ax, color='orange', **kwargs):
        """
        Plot photosphere limits based on uncertainties in fit parameters
        Parameters

        Parameters
        ----------
        ax : matplotlib axis object
        color: matplotlib color identifier.
        kwargs : additional matplotlib kwargs.


        """
        xlam = self.StarPhotosphere_down[0]
        ax.fill_between(xlam * _ANG2MICRON,
                        self.StarPhotosphere_down[1] * xlam,
                        self.StarPhotosphere_up[1] * xlam, color=color,
                        lw=0, **kwargs)

    def plot_photosphere(self, ax, pointsize=4, lcolor='orange', pcolor='orange',
                         marker='o', linestyle='--', lw=2, **kwargs):

        """
        Plot photospheric data.

        Parameters
        ----------
        ax : axis object (e.g. ax = plt.figure().add_subplot(111))
        pointsize (float) : marker size
        color (string) : marker color
        marker (string) : marker type
        linestyle (string) : linesytle
        lw (float): linewidth
        kwargs: matplotlib kw

        """
        xlam, yflux = self.StarPhotosphere

        # PLOT PHOTOSPHERE CONTINUUM
        ax.plot(xlam * _ANG2MICRON, yflux * xlam, color=lcolor, ls=linestyle, lw=lw, **kwargs)

        for i, band in enumerate(self.mags2use):
            if i == 0:
                ax.plot(self.wave[band] * _ANG2MICRON, self.photFlux[band] * self.wave[band],
                        marker=marker, mfc=pcolor, mec=pcolor, ms=pointsize, zorder=10,
                        **kwargs)
            else:
                ax.plot(self.wave[band] * _ANG2MICRON, self.photFlux[band] * self.wave[band],
                        marker=marker, mfc=pcolor, mec=pcolor, ms=pointsize, zorder=10)


    def plot_observedData(self, ax, ms=4, lcolor='k', pcolor='g',
                          markerp='o', markernp='*', capsize=0, linestyle='-',
                          lw=1, **kwargs):

        """
        Plot empirical SED data.

        Parameters
        ----------
        ax : axis object (e.g. ax = plt.figure().add_subplot(111))
        pointsize (float): marker size
        lcolor (string) : line color
        pcolor (string) : point color
        markerp (string) : marker type
        markernp (string) :
        capsize
        linestyle (string): linestyle
        lw (float): linewidth
        """

        xlam, ylam = self.StarPhotosphere
        if self.fullspectrum is not None:
            ylam = self.fullspectrum

        ax.plot(xlam * _ANG2MICRON, ylam * xlam, color=lcolor,
                ls=linestyle, lw=lw, **kwargs)

        for i, (band, lam) in enumerate(list(self.wave.items())):
            flx = self.flux[band + '_flux']
            flxerr = self.fluxerr[band + '_flux']

            if np.any(band == np.array(self.mags4Phot)):
                pfmt = '%s%s' % (pcolor, markerp)
            else:
                pfmt = '%s%s' % (pcolor, markernp)
            if i == 0:
                ax.errorbar(lam * _ANG2MICRON, flx * lam, yerr=flxerr * lam,
                            fmt=pfmt, ms=ms, **kwargs)
            else:
                ax.errorbar(lam * _ANG2MICRON, flx * lam, yerr=flxerr * lam,
                            fmt=pfmt, ms=ms)
