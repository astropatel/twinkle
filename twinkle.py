"""
***********************************************************************
    twinkle.py by Rahul I. Patel (ri.patel272@gmail.com)
***********************************************************************
Issues and changes that need to be made:
 2. Add attributes description to doc string
 3. Fix plot functions to use kwargs
 5. Have code use more than Bt-Vt to interpolate mamajek's file

"""

import os
import copy
import json
import logging
import numpy as np

import sed



logging.basicConfig(filename='example.log', filemode='w', level=logging.DEBUG)

try:
    from astropy import constants as con
except ImportError:
    print 'Does not seem Astropy is installed, or at least the constants package is messed up. We kinda need this. Get to it yo.'

try:
    import matplotlib.pyplot as plt
except ImportError:
    print 'Matplotlib doesnt seem to be installed detected. Fine by me, but now you cant use the awesome' \
          'plotting function we have. Sucks for you.'

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

    def __init__(self, jfile, sid=None, starname=None):
        """
        Instantiates a "twinkling star" object (funny, I know).

        Parameters
        ----------
        jfile: (str) json file of inputs
        sid: (int) index in reference to star in stellar file
        starname: (str) Name of star in MainName column of stellar file
                    Providing both sid and starname will halt execution.
        """

        self.sid = sid
        self.starname = starname

        #  Load data from starfile
        #  Load data from empirical color file
        if not os.path.isfile(jfile):
            raise ValueError('%s aint a file yo'%jfile)
        try:
            script = open(jfile).read()
            specs = json.loads(script)
        except ValueError as Err:
            logging.info('JSON is trippin cause of, ', Err)
            raise ValueError(Err)

        # CHECK TO SEE IF DICTIONARY OF STELLAR DATA HAS BEEN POPULATED
        if not sed.StarsDat:
            # THIS RUNS THE ENTIRE PROCESS OF LOADING AND SAVING THE EMPIRICAL DATA
            sed.DataLogistics(specs)

        self.starsdat = sed.StarsDat
        self.emdat = sed.EmpDat['dat']

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

        #  Find star in file if name given otherwise
        #  find it based on input index.
        if self.sid is None and self.starname is not None:
            ind = np.where(self.starsdat['MainName'] == self.starname)[0]
            if len(ind) == 0:
                raise ValueError("%s was not found in the input file"
                                 % self.starname)
            else:
                self.sid = ind[0]
                self.starname = self.starsdat['MainName'][self.sid]

        elif self.sid is None and self.starname is None:
            raise ValueError('No object index or ID was provided')

        elif self.sid is not None and self.starname is not None:
            raise ValueError('Provide either name of location index, not both.')
        logging.info('WORKING STAR:%s'%self.starname)

        # ADD PHOTOMETRY FROM MAGNITUDE LIST IN JSON FILES
        # REMOVE SATURATED BANDS AND REPLACE NULL VALUES
        # vegaMagDict and errdict ARE FILLED UP HERE.
        self.cleanphotometry(specs)

        #  ========================================
        #  Gather up basic stellar info from file
        #  ========================================
        self.disti = 1000. / self.starsdat['plx'][self.sid]
        self.spti = self.starsdat['spt'][self.sid]
        self.met = 0
        self.modeli = self.starsdat['model'][self.sid]

        self.su2ea = CONST_1 / self.disti ** 2
        self.su2ea_dust = CONST_2 / self.disti ** 2

        #  Obtain grav and met from file or guess based on B-V
        try:
            self.g = self.starsdat['grav'][self.sid]
            self.T0 = self.starsdat['temp'][self.sid] / 1000.

        except KeyError:
            bv = self.emdat['B-V']
            bvi = self.vegaMagDict['BJ'] - self.vegaMagDict['VJ']
            indebv = np.searchsorted(bv, bvi)

            self.T0 = self.emdat['Tinit'][indebv] / 1000.
            self.g = self.emdat['log(g)_x10'][indebv]

        # ================================================================
        #                      Convert Photometry to Flux
        # ================================================================

        self.fluxTup = STools.batch_mag2fluxZPLam(self.mags2use,
                                                  self.vegaMagDict,
                                                  self.vegaMagErrDict)

        self.flux, self.fluxerr = self.fluxTup  # self.fluxTup[0], self.fluxTup[1]
        self.wave = STools.get_eff_wavelengths(self.mags2use)

        # ================================================================
        #                      FIT PHOTOSPHERE
        # ================================================================

        if specs['fitphot']:
            self.fitPhotosphere(self.su2ea, self.T0,
                                self.modeli, self.g, self.met)
            params, perror = self.mfit.params, self.mfit.perror
            self.p_up, self.p_down = params + perror, params - perror


            logging.info('Photosphere for %s fit with T=%s K'
                         %(self.starname, self.StarTemp))

            wave_min, wave_max = specs['spec_sample']['wave_min'], \
                                 specs['spec_sample']['wave_max']
            gridpts = specs['spec_sample']['gridpts']
            modelinfo = (self.modeli, self.g, self.met)
            self.StarPhotosphere = STools.photosphere([self.StarTemp / 1000., self.StarRadius],
                                                      self.su2ea, modelinfo,
                                                      wave=(wave_min, wave_max),
                                                      gridpts=gridpts)

            StarPhot1 = STools.photosphere(self.p_up,
                                           self.su2ea, modelinfo,
                                           wave=(wave_min, wave_max),
                                           gridpts=gridpts)
            StarPhot2= STools.photosphere(self.p_down,
                                          self.su2ea, modelinfo,
                                          wave=(wave_min, wave_max),
                                          gridpts=gridpts)

            flxes = np.array([StarPhot1[1],StarPhot2[1]])
            maxind = np.argmax(np.sum(flxes, axis=1))
            minind = np.argmin(np.sum(flxes, axis=1))

            self.StarPhotosphere_up = (StarPhot1[0],flxes[maxind])
            self.StarPhotosphere_down = (StarPhot1[0],flxes[minind])


            self.StarPhotosphere_unsc = self.StarPhotosphere

            for band in self.mags2use:
                flxt = STools.rsr_flux(getattr(STools, '{}pband'.format(band)),
                                       *self.StarPhotosphere)[0]
                self.photFlux['%s' % band] = flxt

        else:

            logging.info('Photosphere fit not requested for %s.'%self.starname)
            if specs['scalephot']:
                new_phot = STools.scaleSED2bands(specs['phot']['scaleSEDbands'],
                                                 self.mags2use, self.StarPhotosphere[1],
                                                 self.flux, self.fluxerr, self.photFlux)

                phot_norm_fac, yphot, yphot_unsc, RJ_On = new_phot
                self.StarPhotosphere[1] = yphot


    def writeSED(self, filename='sed.txt',
                 comment='# lambda: Angstroms, f_lambda: erg/s/cm^2/Angstrom,\
                          Teff={:.1f}K,rad={:.3f}Rsol.\n'):
        """
        Function to write out fitted SED to file.

        Parameters
        ----------
        filename: (str) file of where SED will be saved.
        comment: (str) comment starting with '#' of any comment to go
                 on first line of the saved SED file.

        """
        comment = comment.format(self.StarTemp, self.StarRadius)

        if not self.StarPhotosphere:
            raise ValueError('Photosphere was not created --> cant be written out.')

        with open(filename,'w') as file:
            file.write('%s' % comment)
            file.write('lambda\t f_lambda\n')
            np.savetxt(file, np.transpose(self.StarPhotosphere), delimiter='\t')


        logging.info('SED for %s saved to %s'%(self.starname,filename))


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

        Parameters
        ----------
        specs: (dict) dictionary created using json file.

        """

        vegaMagDict_temp, vegaMagErrDict_temp = {}, {}

        mags2use0 = copy.copy(specs['phot']['mags2use0_original'])
        mags4Phot0 =  copy.copy(specs['phot']['mags4Phot0_original'])
        mags4scale0 = copy.copy(specs['phot']['mags4scale0_original'])
        mags4Dust0 = copy.copy(specs['phot']['mags4Dust0'])

        #  ========================================
        #  Whether to use W3 and W2 to fit photosphere
        if specs['W3Adapt']:
            mags4Phot0 = self.W3Adopt(specs, mags4Phot0, True)
        if specs['W2Adapt']:
            mags4Phot0 = self.W2Adopt(specs, mags4Phot0, True)

        # Remove optical bands cause of late spectral type?
        if self.starsdat['NoOptical'][self.sid] :

            for arr in Photometry_spCheckList :
                # arr = np.array(arr)
                for mv in specs['phot']['Remove_RedStars']:
                    try:
                        ind = np.where(np.array(eval(arr)) == mv)[0]
                        exec ('%s = np.delete(%s,ind)' % (arr, arr))
                    except ValueError:
                        pass
        # ========================================
        # KEEP ONLY VALID PHOTOMETRIC MEASUREMENTS.
        for mv in mags2use0:

            if self.starsdat['{}m'.format(mv)][self.sid] == 'null':
                try:
                    mags2use0.remove(mv)
                    logging.info('{} band removed from mags2use0'.format(mv))
                except ValueError:
                    loggging.error('Error in removing {} from mags2use'.format(mv))
                try:
                    mags4Phot0.remove(mv)
                    logging.info('{} band removed from mags4Phot0'.format(mv))
                except ValueError:
                    loggging.error('Error in removing {} from mags4phot'.format(mv))
                try:
                    mags4scale0.remove(mv)
                    logging.info('{} band removed from mags4scale0'.format(mv))
                except ValueError:
                    loggging.error('Error in removing {} from mags4scale'.format(mv))
                try:
                    mags4Dust0.remove(mv)
                    logging.info('{} band removed from mags4Dust0'.format(mv))
                except ValueError:
                    loggging.error('Error in removing {} from mags4dust'.format(mv))

            else:
                # CHECK IF IT IS MAG OR FLUX
                if '_flux' in mv:
                    freq = eval('stoo.{}pband.isoFrequency()'.format(mv))
                    fj = float(self.starsdat['%s_flux'%mv][self.sid])
                    # WHEN THERE IS A NULL VALUE
                    try:
                        efj = float(self.starsdat['{}_fluxe'.format(mv)][self.sid])
                    except ValueError:
                        efj = 0.05 * fj
                    fcgs, efcgs= np.array(STools.Jy2cgs((freq, fj), (0, efj))) \
                                 / eval('STools.%spband.isoWavelength()' % mv)

                    vegaMagDict_temp[mv] = fcgs
                    vegaMagErrDict_temp[mv] = efcgs

                # Assumes the rest is a magnitude
                else:
                    vegaMagDict_temp[mv] = float(self.starsdat['{}m'.format(mv)][self.sid])
                    if self.starsdat['%sme' % mv][self.sid] == 'null':
                        vegaMagErrDict_temp[mv] = 0.05 * (vegaMagDict_temp[mv])
                    else:
                        vegaMagErrDict_temp[mv] = float(self.starsdat['{}me'.format(mv)][self.sid])

        # CHECK SATURATION LIMITS AND REMOVE FROM ALL LISTS
        # ========================================
        if specs['satcheck']:
            self.mags2use = self.keep_unsatmags(vegaMagDict_temp, mags2use0)
            self.mags4Dust = self.keep_unsatmags(vegaMagDict_temp, mags4Dust0)
            self.mags4Phot = self.keep_unsatmags(vegaMagDict_temp, mags4Phot0)
            self.mags4scale = self.keep_unsatmags(vegaMagDict_temp, mags4scale0)
        else:
            self.mags2use   = mags2use0
            self.mags4Dust  = mags4Dust0
            self.mags4Phot  = mags4Phot0
            self.mags4scale = mags4scale0


        # REMOVE NON-USED MAGNITUDES FROM DICTIONARY
        #  ========================================
        for mv in mags2use0:
            tmp = mv
            self.vegaMagDict[tmp] = vegaMagDict_temp[tmp]
            self.vegaMagErrDict[tmp] = vegaMagErrDict_temp[tmp]



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
        if not isinstance(magsCheck1,list):
            magsCheck1 = list(magsCheck1)

        mags_temp = np.array(magsCheck).copy()
        for mv in mags_temp:
            try:
                if vegaDict[mv] < eval('self.%s_lim' % mv):
                    magsCheck1.remove(mv)
                    logging.info('%s band is saturated:%.4f<%.4f'
                                 %(mv,vegaDict[mv],eval('self.%s_lim' % mv)))
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
            if self.starsdat['W3m'][self.sid] > 3.8:
                plist = np.append(plist, 'W3')
        else:
            W13_cut = specs['WISE_excess']['W13_cut']
            W23_cut = specs['WISE_excess']['W23_cut']

            w1w3snr = self.starsdat['W1W3SNR'][self.sid]
            w2w3snr = self.starsdat['W2W3SNR'][self.sid]

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
        """Determines whether the WISE W2 photometry
        should be included in the list of  photometric
        points to be used to pin down the photosphere.
        v1 uses the previously determined W1-W2 cuts
        as well as the excess SNRs (Sigma_{E[W1-W2]}) to
        determine whether the photometry is photospheric.

        Parameters
        ----------
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
            raise ValueError('No maglist specified for W2Adopt')
        else:
            plist = pmaglist

        if simple:
            if self.starsdat['W2m'][self.sid] > 3.8:
                plist = np.unique(np.append(plist, 'W2'))


        else:
            W12_cut = specs['WISE_excess']['W12_cut']
            w1w2snr = self.starsdat['W1W2SNR'][self.sid]
            try:
                w2snr = float(w1w2snr)
            except ValueError:
                w2snr = -100  # THIS VALUE IF THERE IS NO SNR VALUE LISTED

            if w2snr <= W12_cut:
                plist = np.append(plist, 'W2')

        return plist

    def fitPhotosphere(self, sconst, T0, model, grav, met):
        """
        Fit a stellar photosphere model to bands in the list/arr
        mags4Phot. This uses mpfit in mosaic_tools.py.

        Results are saved in following internal variables:
        StarRadius: Stellar radius in Solar Rad. units.
        StarTemp: Temperature of fitted star in Kelvin.

        Parameters
        ----------
        sconst: (float) scaling factor for raw fluxes that equals:
                [Solar Radius (cm) / 1 parsec (cm)] ** 2 / [star distance (pc)] ** 2
        T0:     (float) Guess temperature in units of Kelvin/1000
        model:  (str) name of model to use (either Atlas9 or NextGen)
        grav:   (float) Surface gravity in units of 10*log(g [cm/s^2])
        met:    (float) log(metallicity) (0 = solar).

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
                                         mfitlist, STools.calc_grids)

        self.StarRadius, self.StarTemp = fit_dat[0], fit_dat[1]
        self.mfit = fit_dat[2]

    def plot_photrange(self,ax,color='orange',**kwargs):
        """
        Plot photosphere limits based on uncertainties in fit parameters
        Parameters

        Parameters
        ----------
        ax
        color
        kwargs

        Returns
        -------

        """
        xlam = self.StarPhotosphere_down[0]
        ax.fill_between(xlam * _ANG2MICRON,
                        self.StarPhotosphere_down[1]*xlam,
                        self.StarPhotosphere_up[1] *xlam,color=color,alpha=0.4)

    def plot_photosphere(self,ax,pointsize=4,lcolor='orange',pcolor='orange',
                         marker='o',linestyle='--',lw=2,**kwargs):

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
        ax.plot(xlam * _ANG2MICRON, yflux * xlam, color=lcolor, ls=linestyle, lw=lw,**kwargs)

        for band in self.mags2use:
            ax.plot(self.wave[band] * _ANG2MICRON, self.photFlux[band] * self.wave[band],
                    marker=marker,mfc=pcolor,mec=pcolor,ms=pointsize,zorder=10)



    def plot_observedData(self,ax,pointsize=4,lcolor='k', pcolor='g',
                          markerp='o',markernp='*',capsize=0,linestyle='-',lw=1):

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

        Returns
        -------

        """

        xlam, ylam = self.StarPhotosphere
        if self.fullspectrum is not None:
            ylam = self.fullspectrum

        ax.plot(xlam * _ANG2MICRON, ylam * xlam,color=lcolor, ls=linestyle,lw=lw)

        for band,lam in self.wave.iteritems():
            flx = self.flux[band + '_flux']
            flxerr = self.fluxerr[band + '_flux'] * lam

            if np.any(band==np.array(self.mags4Phot)):
                pfmt = '%s%s'%(pcolor,markerp)
            else:
                pfmt= '%s%s'%(pcolor,markernp)

            ax.errorbar(lam * _ANG2MICRON, flx * lam, yerr=flxerr,fmt=pfmt,
                        capsize=capsize,ms=pointsize,zorder=10)
