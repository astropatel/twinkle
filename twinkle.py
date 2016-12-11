"""
***********************************************************************
    twinkle.py by Rahul I. Patel (ri.patel272@gmail.com)
***********************************************************************
Issues and changes that need to be made:
 1. Include Logger.
 2. Add attributes description to doc string
 3. Fix plot functions to use kwargs
 4. Remove readcol.,py
 5. Have code use more than Bt-Vt to interpolate mamajek's file


"""



import copy
import json
import os
import sys
import sed
import numpy as np

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

__author__ = 'Rahul I. Patel <ri.patel272@gmial.com>, Joe Trollo'

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

CONST_1 = (_SOLRAD2CM/ _PC2CM) ** 2
CONST_2 = _AU2CM ** 2 / (4 * np.pi * _PC2CM ** 2)

Photometry_spCheckList = ['mags2use0', 'mags4Phot0', 'mags4scale0']


class StarObject:
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
        assert os.path.isfile(jfile), '%s aint a file yo.' % jfile
        try:
            script = open(jfile).read()
            specs = json.loads(script)
        except ValueError as Err:
            print 'JSON is trippin cause of, ', Err
            raise ValueError(Err)

        if sed.StarsDat is not None:
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
                sys.exit("No STAR named %s was not found"
                         " in the input file" % self.starname)
            else:
                self.sid = ind[0]
                self.starname = self.starsdat['MainName'][self.sid]

        elif self.sid is None and self.starname is None:
            sys.exit('No object index or ID was provided')

        elif self.sid is not None and self.starname is not None:
            sys.exit('Provide either name of location index, not both.')

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

            print 'Photosphere for %s fit with T=%s K' % (self.starname, self.StarTemp)

            wave_min, wave_max = specs['spec_sample']['wave_min'], \
                                 specs['spec_sample']['wave_max']
            gridpts = specs['spec_sample']['gridpts']
            modelinfo = (self.modeli, self.g, self.met)
            self.StarPhotosphere = STools.photosphere([self.StarTemp / 1000., self.StarRadius],
                                                      self.su2ea, modelinfo,
                                                      wave=(wave_min, wave_max),
                                                      gridpts=gridpts)
            self.StarPhotosphere_unsc = self.StarPhotosphere
            for band in self.mags2use:
                flxt = STools.rsr_flux(getattr(STools, '{}pband'.format(band)),
                                       *self.StarPhotosphere)[0]
                self.photFlux['%s' % band] = flxt

        else:

            print 'Photosphere fit not requested.'
            if specs['scalephot']:
                new_phot = STools.scaleSED2bands(specs['phot']['scaleSEDbands'],
                                                 self.mags2use, self.StarPhotosphere[1],
                                                 self.flux, self.fluxerr, self.photFlux)

                phot_norm_fac, yphot, yphot_unsc, RJ_On = new_phot
                self.StarPhotosphere[1] = yphot


    def writeSED(self, filename='sedtest.txt',
                 comment='# lambda: Angstroms, f_lambda: erg/s/cm^2/Angstrom\n'):
        """
        Function to write out fitted SED to file.

        Parameters
        ----------
        filename: (str) file of where SED will be saved.
        comment: (str) comment starting with '#' of any comment to go
                 on first line of the saved SED file.

        """
        assert self.StarPhotosphere, 'Photosphere not created and cant be written out.'
        with open(filename,'w') as file:
            file.write()
            file.write('%s' % comment)
            file.write('lambda\t f_lambda\n')
            np.savetxt(file, np.transpose(self.StarPhotosphere), delimiter='\t')


        print 'SED saved to %s' % filename


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
        if self.starsdat['NoOptical'][self.sid] == 'Yes':
            for arr in Photometry_spCheckList:
                arr = np.array(arr)
                for mv in specs['phot']['Remove_RedStars']:
                    try:
                        ind = np.where(eval(arr) == mv)[0]
                        exec ('%s = np.delete(%s,ind)' % (arr, arr))
                    except ValueError:
                        pass
        # ========================================
        # KEEP ONLY VALID PHOTOMETRIC MEASUREMENTS.
        for mv in mags2use0:
            tmp = mv

            if self.starsdat['%sm' % mv][self.sid] == 'null':
                try:
                    mags2use0.remove(mv)
                except ValueError:
                    print 'Error in removing %s from mags2use' % mv
                try:
                    mags4Phot0.remove(mv)
                except ValueError:
                    print 'Error in removing %s from mags4phot' % mv
                try:
                    mags4scale0.remove(mv)
                except ValueError:
                    print 'Error in removing %s from mags4scale' % mv
                try:
                    mags4Dust0.remove(mv)
                except ValueError:
                    print 'Error in removing %s from mags4dust' % mv

            else:
                vegaMagDict_temp[tmp] = float(self.starsdat['%sm' % mv][self.sid])
                if self.starsdat['%sme' % mv][self.sid] == 'null':
                    vegaMagErrDict_temp[tmp] = 0.05 * (vegaMagDict_temp[tmp])
                else:
                    vegaMagErrDict_temp[tmp] = float(self.starsdat['%sme' % mv][self.sid])

        # CHECK SATURATION LIMITS AND REMOVE FROM ALL LISTS
        # ========================================
        self.mags2use = self.keep_unsatmags(vegaMagDict_temp, mags2use0)
        self.mags4Dust = self.keep_unsatmags(vegaMagDict_temp, mags4Dust0)
        self.mags4Phot = self.keep_unsatmags(vegaMagDict_temp, mags4Phot0)
        self.mags4scale = self.keep_unsatmags(vegaMagDict_temp, mags4scale0)

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
        mags_temp = np.array(magsCheck).copy()
        for mv in mags_temp:
            try:
                if vegaDict[mv] < eval('self.%s_lim' % mv):
                    magsCheck1.remove(mv)
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
            sys.exit('No maglist specified for W3Adopt')
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
            sys.exit('No maglist specified for W2Adopt')
        else:
            plist = pmaglist

        if simple:
            if self.starsdat['W2m'][self.sid] > 3.8:
                plist = np.append(plist, 'W2')

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
        assert np.array(self.mags4scale).size, 'mags4scale is empty.'
        assert np.array(self.mags4Phot).size, 'mags4Phot is empty.'

        mfitlist = {'photmags': np.array(self.mags4Phot),
                    'scalemags': np.array(self.mags4scale)}

        rawfluxdat = (self.flux, self.fluxerr)

        fit_dat = STools.fit_photosphere(self.wave, rawfluxdat,
                                         [T0], sconst, modeltype,
                                         mfitlist, STools.calc_grids)

        self.StarRadius, self.StarTemp = fit_dat[0], fit_dat[1]
        self.mfit = fit_dat[2]


    def plot_photosphere(self,ax,pointsize=4,lcolor='orange',pcolor='orange',
                         marker='o',linestyle='--',lw=2):

        """

        Parameters
        ----------
        ax
        pointsize
        color
        marker
        linestyle
        lw

        """
        xlam, yflux = self.StarPhotosphere

        # PLOT PHOTOSPHERE CONTINUUM
        ax.plot(xlam * _ANG2MICRON, yflux * xlam, color=lcolor, ls=linestyle, lw=lw)

        for band in self.mags2use:
            ax.plot(self.wave[band] * _ANG2MICRON, self.photFlux[band] * self.wave[band],
                    marker=marker,mfc=pcolor,ms=pointsize)



    def plot_observedData(self,ax,pointsize=4,lcolor='k', pcolor='g',
                          markerp='o',markernp='*',capsize=0,linestyle='-',lw=1):

        """
        
        Parameters
        ----------
        ax
        pointsize
        lcolor
        pcolor
        markerp
        markernp
        capsize
        linestyle
        lw

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
                        capsize=capsize,ms=pointsize)


