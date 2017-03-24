# STILL REQUIRES ADDENDUM TO INCLUDE EMPIRICAL2SYNTHETIC OFFSET
# ADJUSTMENT

#  =========================================================
import os, sed, random
import argparse, copy, json
import matplotlib.pyplot as plt, numpy as np
import matplotlib.ticker as mtick, math as ma
from scipy.optimize import bisect
import mosaic_tools as mt
from readcol import *
from matplotlib.ticker import MaxNLocator, MultipleLocator
import load_wfcorrection as lwf
from mpl_toolkits.axes_grid1 import Grid
import sed_paramfile as sp


try:
    from astropy import constants as con
except ImportError:
    print 'Does not seem Astropy is installed, or at least the constants package is messed up. We kinda need this. Get to it yo.'


__author__ = 'Rahul I. Patel <ri.patel272@gmail.com>'
#  =========================================================
#     Definition to obtain Tbb assuming flux correction
#  =========================================================




def find_RJPoints(wavedict, photkeys):
    cutoff = 20000  # 20000.
    keys = wavedict.keys()
    values = wavedict.values()
    waveVals = list(at.dict2list(wavedict, photkeys))
    key4RJ = []

    while len(key4RJ) != 2:
        if len(waveVals) == 0:
            break
        ind = np.where(waveVals == np.array(waveVals).max())[0][0]
        vali = waveVals.pop(ind)
        if vali > cutoff:
            ind2 = np.where(vali == values)[0][0]
            key_i = keys[ind2]
            key4RJ.append(key_i)

    return key4RJ

    
def plot_photosphere(ax, mg4d, xphot, yphot, wave,
                     PhotDust_Flux, a2m, pointsize, lw):

    ax.plot(xphot * a2m, yphot * xphot, 'b--', lw=lw)
    mags4Dust2 = mg4d

    for mv in mags4Dust2:
        ax.plot(wave[mv] * a2m, PhotDust_Flux[mv] * wave[mv], 'bo', ms=pointsize)

    return


def plot_blackbody(ax, wave, ExcessFlux, dust_lambda, lightW4_line, a2m, ptsz, lw):
    # print Exbool_Dict
    ax.plot(dust_lambda * a2m, lightW4_line * dust_lambda, 'm-.', linewidth=lw)
    for key, val in ExcessFlux.iteritems():
        if Exbool_Dict[key]:
            ax.errorbar(wave[key] * a2m, val * wave[key],
                        yerr=0.7 * val * wave[key], uplims=True, ecolor='red',
                        capsize=5, elinewidth=5, capthick=2)
        else:
            ax.errorbar(wave[key] * a2m, val * wave[key], yerr=ExcessFluxerr[key] * wave[key], \
                        label='Observed-Photosphere', fmt='mD', capsize=10, ms=ptsz, mfc='white', mec='magenta',
                        mew=1.5)
    return


def plot_observedData(ax, mg4p, xphot, fullspectrum, wave,
                      flux, fluxerr, a2m, ptsze, cpsze, lw):

    ax.plot(xphot * a2m, fullspectrum, 'k-', lw=lw, label='Full Spectrum')
    for band, lam in wave.iteritems():
        flx = flux[band + '_flux']
        flxerr = fluxerr[band + '_flux'] * lam

        if np.any(band == np.array(mg4p)):
            pt_fmt = 'go'

        else:
            pt_fmt = 'g*'
        ax.errorbar(lam * a2m, flx * lam, yerr=flxerr, fmt=pt_fmt, capsize=cpsze, ms=ptsze)

    return


def plot_annotations(ax, star, spti, tempStar, tempDust,
                     beta, p0, W3Excess_bool, Exfunc, fontsize, minorftsize):
    try:
        tempStar = int(round(tempStar, 0))
    except TypeError:
        pass
    ax.annotate(r'%s' % (star), xy=(0.13, .33),
                xycoords='axes fraction',
                fontsize=fontsize, family='Times New Roman')

    ax.annotate(r'%s' % (spti), xy=(0.13, 0.26),
                xycoords='axes fraction',
                fontsize=minorftsize, family='Times New Roman')

    ax.annotate(r'$T_{*} = $%sK' % (tempStar), xy=(0.13, 0.19),
                xycoords='axes fraction',
                fontsize=fontsize, family='Times New Roman', color='blue')
    # ax.annotate(r'$\chi^2$ = %.1f'%(chi2),xy = (0.48, .63),
    #            xycoords = 'axes fraction',
    #            fontsize=ftsize, family='Times New Roman',color='blue')

    ax.annotate(r'$T_{BB} = $%dK' % (tempDust), xy=(0.13, 0.12),
                xycoords='axes fraction',
                fontsize=fontsize, family='Times New Roman', color='magenta')

    if W3Excess_bool or Exfunc == 'modifiedBB':
        beta = p0[2]
        ax.annotate(r'$\beta=$%.1f' % beta, xy=(0.13, 0.05),
                    xycoords='axes fraction',
                    fontsize=fontsize, family='Times New Roman', color='red')
    else:
        pass

    return


def plot_annotationsW2(ax, star, spti, tempStar, tempDust, beta, p0,
                       W3Excess_bool, Exfunc, fontsize, minorftsize):
    try:
        tempStar = int(round(tempStar, 0))
    except TypeError:
        pass
    ax.annotate(r'%s' % (star), xy=(0.65, .90),
                xycoords='axes fraction', ha='right', va='top',
                fontsize=fontsize, family='Times New Roman')

    ax.annotate(r'%s' % (spti), xy=(0.90, 0.90),
                xycoords='axes fraction', ha='right', va='top',
                fontsize=minorftsize, family='Times New Roman')

    ax.annotate(r'$T_{*} = $%sK' % (tempStar), xy=(0.90, 0.80),
                xycoords='axes fraction', ha='right', va='top',
                fontsize=fontsize, family='Times New Roman', color='blue')
    # ax.annotate(r'$\chi^2$ = %.1f'%(chi2),xy = (0.48, .63),
    #            xycoords = 'axes fraction',
    #            fontsize=ftsize, family='Times New Roman',color='blue')

    ax.annotate(r'$T_{BB} = $%dK' % (tempDust), xy=(0.90, 0.70),
                xycoords='axes fraction', ha='right', va='top',
                fontsize=fontsize, family='Times New Roman', color='magenta')

    if W3Excess_bool or Exfunc == 'modifiedBB':
        beta = p0[2]
        ax.annotate(r'$\beta=$%.1f' % beta, xy=(0.90, 0.60),
                    xycoords='axes fraction', ha='right', va='top',
                    fontsize=fontsize, family='Times New Roman', color='red')
    else:
        pass

    # try:
    #     ax.annotate(r'$R_{d} \simeq $%d AU, $\tau$ = %.1f'%(radius_dust,ma.log10(tau)), xy=(0.25,0.12),
    #               xycoords='axes fraction',
    #               fontsize=ftsize, family = 'Times New Roman',color='red')
    # except ValueError:
    #     ax.annotate('No Excess', xy=(0.25,0.12),
    #               xycoords='axes fraction', fontsize=ftsize, family = 'Times New Roman',color='red')
    return


parser = argparse.ArgumentParser(description='Flags to run sed plotting')
#  Let multiple star names to be added to flag "-s" to be run in the sed plot code
parser.add_argument('-s', nargs='*', help='Identifiers of stars to plot')
parser.add_argument('--noshow', action="store_true", help='Used to suppress onscreen plot')
parser.add_argument('--savf', action="store_true", help='Used to save figures')
parser.add_argument('--nw', action="store_true", help='Used to suppress file writing of model summary')
p_args = parser.parse_args()

#  ====================================
#        SETUP
#  ====================================
np.seterr(all='ignore')

ft = mt.FittingTools()
at = mt.ArrayTools()
PT = mt.PlottingTools()
STools = __import__('sed').SEDTools()

W13_cut = sp.W13_cut
W23_cut = sp.W23_cut

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

# CONST_1 & CONST_2 used for fitting in loop

CONST_1 = (_SOLRAD2CM / _PC2CM) ** 2
CONST_2 = _AU2CM ** 2 / (4 * np.pi * _PC2CM ** 2)
#  ====================================
#             CONSTANTS
#  ====================================

wave_min = sp.wave_min * _MICRON2ANG
wave_max = sp.wave_max * _MICRON2ANG
# PLOT LIMITS IN erg/s/cm^2
ylim_up, ylim_low = sp.ylim_up, sp.ylim_low
xmax = sp.xmax   # PLOT LIMITS IN MICRONS
gridPad = 0.0  # GRID PADDING BETWEEN CELLS
plotsize = sp.plotsize


#  =======================================================================================
#                                    FILES
#  =======================================================================================
timeNOW = strftime('%Y%m%d_%H%M%S', gmtime())
file_write = os.path.join(os.getcwd(), 'DebrisDisks', 'XMatch', 'FullSky',
                          'disk_stars_info_%s.txt' % timeNOW)
ran = str(random.randint(1, 100))[0:4]

if sp.write2file and not p_args.nw:
    f = open(file_write, 'w')

    header = 'Object \t W3flux_cgs \t W3eflux_cgs \t W3PhotFlux_cgs \t W4flux_cgs \t W4eflux_cgs \t W4PhotFlux_cgs ' + \
             '\t ExW3Corrected_cgs \t ExW4Corrected_cgs \t W3ExcessFraction \t W4ExcessFraction \t RelativeW3Flux ' + \
             '\t RelativeW4Flux \t e_RelativeW4Flux \t e_RelativeW3Flux \t W3flux_mJy \t W3eflux_mJy ' + \
             '\t W3PhotFlux_mJy \t W4flux_mJy \t W4eflux_mJy \t W4Photflux_mJy \t W3Upperlim? \t W4Upperlim? ' + \
             '\t KcorW3 \t KcorW4 \t Log(L) \t Tstar \t Rstar \t chi2_* \t Rdust \t Tdust_calc \t fd\n'

    f.write('%s' % header)
    print 'Writing to: ', os.path.basename(file_write)

else:
    pass
#  -------------------------------------------------------------------------------------------
dfbvst = readcol(sp.fileBVStandards, asdict=True, verbose=False)
sptST, teffST, logLST = dfbvst['SpT'], dfbvst['Teff'], dfbvst['logL']
# =======================================================================================
#                                  DATA
# =======================================================================================
script = open('/Users/rpatel/Dropbox/Research/sed_paramfile.json').read()
specs_from_file = json.loads(script)

sed.DataLogistics(specs_from_file, changekeys=True)
stdat = sed.StarsDat

star_arr = stdat['MainName']
plx = stdat['plx']
dist = 1000. / plx

W3OnlyExcessFlagArr = np.array(['NNNYNN', 'NNNNYN', 'NNNYYN', 'NNNYNY', 'NNNYYY'])

#  =======================================================================================
#   SETUP BANDS: THIS TELLS THE PROGRAM WHICH PHOTOMETRIC BANDS
#   TO USE IN ORDER TO DO THE FITTING, ETC.
#  
#   mags2use: all the ones that will be used ever
#   mags2Phot: ones to use to fit the photosphere
#   mags4scale: ones to use in order to scale the surface to earth photospheric flux
#   mags4Dust: to use to fit the dust blackbody
#  =======================================================================================

# print 'WFC Start'
# wfc = wise_flux_correction.WISECorrect(bands=mags4Dust0)
# print 'WFC DONE'
#  =======================================================================================
#   SET UP WHICH STARS TO USE AND INITIAL METALLICITY AND GRAVITY
#   THE LATTER ARE DUMMY VARIABLES USED IN STARTUP. CODE ADAPTS TO INDIVIDUAL STELLAR
#   PARAMETERS. PREFERENCE IS TO SORT FILE WRT "MODELNAME,METALLICITY,GRAVITY" TO
#   INCREASE EFFICIENCY OF MULTIPLE STARS
#  =======================================================================================
if p_args.s is None:
    try:  # LOAD FROM FILE
        df_selectstars = readcol(sp.file_select, asdict=True, verbose=False)
        select_stars = df_selectstars['NAME']
    except:
        select_stars = np.array([])
else:  # command line input
    select_stars = p_args.s

f_ind = np.array([])
#  ADD STARS TO "SELECT_STARS" AND ONLY THOSE STARS WILL BE FITTED
if len(select_stars) != 0:
    # STARS WILL BE SELECTED BASED ON INDEX IN ARRAY IF ONLY CERTAIN STARS WILL BE USED, THEIR INDICES WILL BE
    # CATALOGUED.
    for i in range(len(select_stars)):
        f_indi = np.where(select_stars[i] == star_arr)[0]
        f_ind = np.append(f_ind, f_indi)
else:

    f_ind = np.arange(len(star_arr))

nobbFit = np.array([])

print '=============================================================\n'
print '        Fitting has begun. Enjoy the experience.\n'

plot_i = 0  # INITIALIZE PLOT CELL

if sp.plot_grid:
    fig = plt.figure(figsize=plotsize)
    ax2 = fig.add_subplot(111)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])

    grid = Grid(fig, rect=111, nrows_ncols=(sp.pROW, sp.pCOL), axes_pad=gridPad, label_mode='L')
    gridnum = 1  # TO ADD TO THE FILE NUMBER

for i, useind in enumerate(f_ind):

    #  =====================================================================================
    #                 Set Up individual stellar parameters
    #  =====================================================================================
    i = int(useind)
    star = star_arr[i]
    print "\nFitting for %s" % star
    spti = stdat['spt'][i]
    g, met = stdat['grav'][i], stdat['met'][i]
    modeli = stdat['model'][i]
    T0 = stdat['temp'][i] / 1000.
    disti = dist[i]

    SOBJ = sed.StarObject(stdat, i)
    # (1/Dist)^2 factor to be multiplied with Radius^2. Unitless.
    # Takes into account solar units, so radius only needs to be in solar units
    # This is to facilitate the fitting procedure.

    su2ea2 = CONST_1 / disti ** 2
    su2ea_dust = CONST_2 / disti ** 2
    b1RJ, b2RJ = 'N/A/', 'N/A'

    p0 = np.array([T0])  # FOR PHOTOSPHERE
    #  =====================================================================================
    #                 Reset Photometry Choices to ALL
    #  =====================================================================================

    mags2use0 = copy.copy(sp.mags2use0_original)
    mags4Phot0 = copy.copy(sp.mags4Phot0_original)
    mags4scale0 = copy.copy(sp.mags4scale0_original)

    #  =====================================================================================
    #                 Check whether to use W2 and/or W3 photometry
    #  =====================================================================================

    if sp.W3Adapt:
        mags4Phot0 = SOBJ.W3Adopt(mags4Phot0, True)

    if sp.W2Adapt:
        mags4Phot0 = SOBJ.W3Adopt(mags4Phot0, True)

    Photometry_spCheckList = ['mags2use0', 'mags4Phot0', 'mags4scale0']

    #  Spectral type check: Don't use listed mags in spRemove_RedStars
    #  modulize this afterwards.
    #  SPTcheck = spti[0]
    if stdat['NoOptical'][i] == 'Yes':
        for arr in Photometry_spCheckList:
            for mv in sp.Remove_RedStars:
                try:
                    ind = np.where(eval(arr) == mv)[0]
                    exec ('%s = np.delete(%s,ind)' % (arr, arr))
                except ValueError:
                    pass

    # =====================================================================================
    #                      CLEAN UP PHOTOMETRY (SATURATION AND NULLS)
    #  =====================================================================================

    #  CLEAN UP NULLS AND REMOVE SATURATED STARS
    #  CREATES DICTIONARY AND FILTERED LISTS
    SOBJ.cleanphotometry()

    #  =====================================================================================
    #             CHANGE PHOTOMETRY BASED ON ZERO POINT OFFSETS FROM COLOR TRENDS
    #  =====================================================================================

    # CHECK SATURATION LIMITS TO USE FOR DUST FITTING
    if star == 'HIP117972':
        SOBJ.vegaMagDict['W2'] += 0.02
        SOBJ.vegaMagDict['W3'] -= 0.05
        SOBJ.vegaMagDict['W4'] -= 0.01

    # =====================================================================================
    #                      Convert Photometry to Flux
    #  =====================================================================================

    # ALL OF THESE ARE DICTIONARIES -- FLUX IS IN erg/s/cm^2/Angstrom, wave in Angstrom
    fluxTup = STools.batch_mag2fluxZPLam(SOBJ.mags2use, SOBJ.vegaMagDict,
                                         SOBJ.vegaMagErrDict)
    flux, fluxerr = fluxTup[0], fluxTup[1]
    wave = STools.get_eff_wavelengths(SOBJ.mags2use)

    #  =====================================================================================
    #                   Begin Stellar Photosphere Fit
    #  =====================================================================================

    modeltype = (modeli,g,met)
    mfitlist = {'photmags':mags4Phot0,
                'scalemags':SOBJ.mags4scale}
    rawfluxdat = (flux,fluxerr)

    fit_dat = STools.fit_photosphere(wave, rawfluxdat, p0, su2ea2,
                           modeltype, mfitlist, STools.calc_grids)

    radius, tempnew = fit_dat[0],fit_dat[1]
    p0 = [round(tempnew) / 1.e3, radius]
    #p0 = [8840./ 1.e3, radius]
    #  =====================================================================================
    #                    Calculate Photosphere Line
    #  =====================================================================================

    # FULL BBODY - ANGSTROMS
    # set conv = 1/_ANG2MICRON
    # convert input wavelength into angstrom
    # XPHOT: Angstrom, YPHOT: erg s^-1 cm^-2 A^-1, same with slope and yint
    xphot, yphot = STools.photosphere(p0, su2ea2, modeli,
                                      wave=(wave_min, wave_max),gridpts=sp.gridpts)




    # **********************************************************************************
    #                   SCALE SED TO WISE FLUXES
    # **********************************************************************************
    # CREATE DICTIONARY OF SYNTHETIC FLUXES FROM NEW PHOTOSPHERE
    synpFlux = {}  # synthetic photospheric flux
    for band in SOBJ.mags2use:

        flxt = STools.rsr_flux(eval('STools.%spband' % band), xphot, yphot)[0]
        synpFlux['%s' % band] = flxt

    dat = STools.scaleSED2bands(sp.scaleSEDbands, SOBJ.mags2use, yphot,
                                flux, fluxerr, synpFlux)
    norm_wise_nir, yphot, yphot_unsc, RJ_On = dat

    # **********************************************************************************
    #                   WRITES SED TO FILE
    # **********************************************************************************
    if sp.write_SED:
        print 'Saving SED information for %s\n' % star
        SEDWrite = os.path.join(os.getcwd(), 'DebrisDisks', star + '_SED_%i.txt' % tempnew)
        fSED = open(SEDWrite, 'w')
        fSED.write('# lambda: Angstroms, f_lambda: erg/s/cm^2/Angstrom\n')
        fSED.write('lambda\t f_lambda\n')
        np.savetxt(fSED, np.transpose((xphot, yphot)), delimiter='\t')
        fSED.close()

    #  =====================================================================================
    #                    Calculate BLACKBODY PARAMETERS
    #  =====================================================================================

    # ===============SET UP BOOLEANS AND EXTRA FLUXES THAT MAY BE ADDED IN =================

    dust_lambda = xphot  # np.logspace(ma.log10(wave_min),ma.log10(wave_max),sp.gridpts)

    calcBB_bool = False
    fitBB_bool = False
    scaleBB_bool = False
    Nothing_bool = False
    dontdraw_bool = False
    W3Excess_bool = False
    i_other = int(useind)
    # CHECK IF THERE ARE LONGER WAVELENGTH FLUXES TO FIT THE SED WITH
    if stdat['OtherFlux'][i_other] != 'None' and sp.Longwave_Bool:
        OtherBands = stdat['OtherFlux'][i_other]
        band_split = OtherBands.split(',')  # ARRAY WITH NAMES FOR BANDS TO BE ADDED

        for band in band_split:
            SOBJ.mags2use.append(band)
            SOBJ.mags4Dust.append(band)
            Oflux = df[band + '_flux'][i_other].split('pm')
            flux[band + '_flux'] = float(Oflux[0])
            fluxerr[band + '_flux'] = float(Oflux[1])
            wave[band] = df[band + '_lam'][i_other] * _MICRON2ANG

    else:
        pass
    # ===============CONTINUE WITH EVERYTHING ELSE =================

    # CALCULATE PHOTOSPHERIC FLUX AT 12 AND 24 MICRONS THROUGH WISE FILTER:
    # flux: [erg s-1 cm-2 A-1], wavelength: [Angstrom]

    # ====================================================================================
    PhotDust_Flux = {}  # Photospheric flux for wavelengths in mags4Dust
    BBDust_Flux = {}  # Flux of blackbody convolved with bands in mags4Dust
    ExcessFluxerr = {}  # 1sigma errors of Excess flux for bands in mags4Dust
    ExcessFlux = {}  # Excess flux for bands in mags4Dust
    Lam_Excess = {}  # wavelength dictionary of bands in mags4Dust
    Exbool_Dict = {}  # TRUE if excess measurement is 3sig upperlimit, otherwise it's not

    beta = 1.0  # FILLER PLACEMENT
    # ========CHECK HERE TO SEE IF W3/W4 FLUXES NEED TO BE PUSHED TO 3SIG LIMIT ==========
    mags4DustTemp = np.array(SOBJ.mags4Dust).copy()
    for mv in mags4DustTemp:
        svMaBool = True

        PhotDust_tmp = STools.rsr_flux(eval('STools.' + mv + 'pband'), xphot, yphot)[0]
        ExFlux_tmp = flux[mv + '_flux'] - PhotDust_tmp

        if (ExFlux_tmp < 0):  # IF EXCESS ID NEGATIVE
            ExFlux_tmp = (flux[mv + '_flux'] + 3 * fluxerr[mv + '_flux']) - PhotDust_tmp
            Exbool_Dict[mv] = False  # remove after fitting HIP21547

            if ExFlux_tmp < 0:
                SOBJ.mags4Dust.remove(mv)
                print mv, 'removed from mags4Dust'
                svMaBool = False
            else:
                Exbool_Dict[mv] = True  # IDENTIFY THAT IT'S A 3SIGMA UPPER LIMIT
                svMaBool = True
                # mags4Dust.remove(mv)

        else:  # IF EXCESS WAS ALWAYS POSITIVE
            Exbool_Dict[mv] = False

        if svMaBool:  # STORE DATA FOR THIS BAND IF POSITIVE/3SIG UPPER LIMIT POSITIVE
            PhotDust_Flux[mv] = PhotDust_tmp
            ExcessFlux[mv] = ExFlux_tmp
            ExcessFluxerr[mv] = fluxerr[mv + '_flux']
            Lam_Excess[mv] = wave[mv]
    N_Excess = np.array(ExcessFlux.values())
    mg4d = np.array(SOBJ.mags4Dust)

    #  PhotDust_Flux['W1'] = STools.rsr_flux(eval('STools.W1pband'), xphot, yphot)[0]
    #  PhotDust_Flux['W2'] = STools.rsr_flux(eval('STools.W2pband'), xphot, yphot)[0]
    #  ========DEPENDING ON ABOVE, DETERMINE IF EXCESSES NEED TO BE FIT OR SCALED TO BB=====

    if np.any(nobbFit == star):
        W3Excess_bool = True

    elif len(np.where(N_Excess > 0)[0]) == 0 or len(SOBJ.mags4Dust) == 0:
        print 'There are no WISE mags to fit BB or no positive excess values for %s' % star
        Nothing_bool = True

    elif len(SOBJ.mags4Dust) == 1 or \
            (np.any(mg4d == 'W3') == False and np.any(mg4d == 'W2') and np.any(mg4d == 'W4')):
        scaleBB_bool = True  # SCALE A SINGLE TEMP BLACKBODY TO SINGLE DATA POINT
        try:
            #  mags4Dust.remove('W2')
            print 'W2 removed from mags4Dust'
        except ValueError:
            pass
        sp.Exfunc = 'blackbody'

    else:
        tmin, tmax = 10, stdat['temp'][i]
        # SORT BY WAVELENGTH INCREASING
        sortedBands = np.array(sorted(Lam_Excess.items(), key=lambda x: x[1]))
        bandSorted = sortedBands[:, 0]
        lamSorted = sortedBands[:, 1].astype('float32') * _ANG2CM
        flxSorted = []
        testExFlux = []

        if len(SOBJ.mags4Dust) > 3:
            fitBB_bool = True
            # Exfunc = 'blackbody'


        elif len(SOBJ.mags4Dust) <= 3 and len(SOBJ.mags4Dust) > 0:
            for bandS in bandSorted:
                flxSorted.append(ExcessFlux[bandS])
                bandSTest = STools.cgs2Jy(wave[bandS], eval('STools.%spband.isoFrequency()' % bandS), \
                                          ExcessFlux[bandS])
                testExFlux.append(bandSTest)
            flxSorted, testExFlux = np.array(flxSorted), np.array(testExFlux)
            #  THE FOLLOWING ASSUMES EITHER 2 OR 3 EXCESS FLUX POINTS TO CONSTRAIN DUST
            tef = testExFlux

            if len(SOBJ.mags4Dust) == 3 and sp.Exfunc == 'blackbody':
                try:
                    if (tef[0] < tef[1]) or (tef[0] > tef[1] > tef[2]):
                        fitBB_bool = True
                        # Exfunc = 'blackbody'

                    else:
                        SOBJ.mags4Dust.remove('W2'), PhotDust_Flux.pop('W2')  # , BBDust_Flux.pop('W2')
                        ExcessFluxerr.pop('W2'), ExcessFlux.pop('W2')
                        Lam_Excess.pop('W2'), Exbool_Dict.pop('W2')
                        bandSorted = list(bandSorted)
                        bandSorted.remove('W2')
                        bandSorted = np.array(bandSorted)
                        lamSorted, flxSorted, testExFlux = lamSorted[1:], flxSorted[1:], testExFlux[1:]
                except IndexError:
                    pass


            elif len(SOBJ.mags4Dust) == 3 and sp.Exfunc == 'modifiedBB':
                fitBB_bool = True

            if len(SOBJ.mags4Dust) == 2:
                bandlow, bandhi = bandSorted[0], bandSorted[1]
                lam1, lam2 = lamSorted  # Lam_Excess[bandlow]*_ANG2CM, Lam_Excess[bandhi]*_ANG2CM
                # flx1,flx2 = ExcessFlux[bandlow], ExcessFlux[bandhi]
                # arr1,arr2 = np.array([lam1,lam2]),np.array([flx1,flx2])
                arr1, arr2 = lamSorted, flxSorted  # np.array([lam1,lam2]),np.array([flx1,flx2])
                fa, fb = STools.calcBBTemp(tmin, arr1, arr2), STools.calcBBTemp(tmax, arr1, arr2)
                # W3Ex,W4Ex = STools.cgs2Jy(wave[bandlow],STools.W3pband.isoFrequency(),ExcessFlux[bandlow]),\
                #             STools.cgs2Jy(wave[bandhi],STools.W4pband.isoFrequency(),ExcessFlux[bandhi])

                # if (fa/fb)<0 and (W3Ex/W4Ex)>0 and len(mags4Dust)==2:
                if (fa / fb) < 0 and (testExFlux[0] / testExFlux[1]) > 0:  # W3Ex/W4Ex
                    calcBB_bool = True
                    sp.Exfunc = 'blackbody'

                elif (fa / fb) > 0:
                    fitBB_bool = True
                    beta = 1.0
                    #  for mv in magPSEDust:# HERE USE 3SIGMA UPPERLIMIT
                    #  Exbool_Dict[mv] = True
                    #  PhotDust_Flux[mv] = STools.rsr_flux(eval('STools.'+mv+'pband'), xphot, yphot)[0]
                    #  ExcessFlux[mv] = (flux[mv+'_flux'] + 3*fluxerr[mv+'_flux']) - PhotDust_Flux[mv]
                    #  Lam_Excess[mv] = wave[mv]
                    #  ExcessFluxerr[mv] = fluxerr[mv+'_flux']
                    #  fitBB_bool = True
                    #  mags4Dust = mags4Dust + magPSEDust
                    #  beta = 0

            else:
                pass

        #      ========================START FIT OR SCALING ========================================

    sortedBands = np.array(sorted(Lam_Excess.items(), key=lambda x: x[1]))
    bandSorted = sortedBands[:, 0]
    lamSorted = sortedBands[:, 1].astype('float32') * _ANG2CM
    # bandlow, bandhi = bandSorted[0],bandSorted[1]
    flxSorted = []
    testExFlux = []

    PhotDust_Fluxarr = np.array(at.dict2list(PhotDust_Flux, SOBJ.mags4Dust))
    Excess_Flxerrarr = np.array(at.dict2list(fluxerr, SOBJ.mags4Dust, '_flux'))
    Excess_Flxarr = np.array(at.dict2list(ExcessFlux, SOBJ.mags4Dust))
    print sp.Exfunc

    if fitBB_bool:  # IF THE BLACKBODY IS GOING TO BE FIT

        try:
            lam0 = wave[lam0band]
        except:
            lam0 = wave['W2']  # REFERENCE WAVELENGTH USED FOR MODIFIED BLACKBODY
        print 'lam0=W2'

        lamMin = min(Lam_Excess.values())
        tempdust = STools.wienTEMP(lamMin, units='angstrom')
        for mv in SOBJ.mags4Dust:
            BBDust_Flux[mv] = eval(
                'STools.' + sp.Exfunc + '(dust_lambda, np.array([tempdust]),1,np.array([mv]),beta=beta,lam0=lam0)')
        # mv = 'W2'

        BBDust_Fluxarr = np.array(at.dict2list(BBDust_Flux, SOBJ.mags4Dust))
        indNoNeg = np.where(Excess_Flxarr > 0)[0]
        # if len(indNoNeg) != 0:
        FluxNorm_dust = np.average(Excess_Flxarr[indNoNeg] / BBDust_Fluxarr[indNoNeg],
                                   weights=1. / Excess_Flxerrarr[indNoNeg])
        Rad_dust = ma.sqrt(FluxNorm_dust / su2ea_dust)

        if sp.Exfunc == 'blackbody':
            p0_dust = np.array([tempdust, Rad_dust])
            nparams = 2
            fa_Dust = {'x': dust_lambda, 'y': Excess_Flxarr, 'err': Excess_Flxerrarr,
                       'func': eval('STools.' + sp.Exfunc), 'su2ea1': su2ea_dust, 'bands': SOBJ.mags4Dust}
            parinfo_dust = [{'value': 0., 'relstep': 0, 'limits': [0, 0], 'limited': [0, 0], 'fixed': 0} for m in
                            range(nparams)]
            for k in range(nparams): parinfo_dust[k]['value'] = p0_dust[k]
            parinfo_dust[0]['relstep'] = 0.3
            parinfo_dust[1]['relstep'] = 0.2

        elif sp.Exfunc == 'modifiedBB':

            p0_dust = np.array([tempdust, Rad_dust, beta])
            nparams = 3
            fa_Dust = {'x': dust_lambda, 'y': Excess_Flxarr, 'err': Excess_Flxerrarr,
                       'func': eval('STools.' + sp.Exfunc), 'su2ea1': su2ea_dust,
                       'bands': SOBJ.mags4Dust, 'lam0': lam0}
            parinfo_dust = [{'value': 0., 'relstep': 0, 'limits': [0, 0], 'limited': [0, 0], 'fixed': 0} for m in
                            range(nparams)]
            for k in range(nparams): parinfo_dust[k]['value'] = p0_dust[k]
            parinfo_dust[2]['limits'] = [0, 10]
            parinfo_dust[0]['relstep'] = 0.1
            parinfo_dust[1]['relstep'] = 0.1

        else:
            print 'Specify function to fit to data for Excess emission'
            sys.exit()

        m_dust = mt.mpfit(ft.deviates_from_model, parinfo=parinfo_dust, functkw=fa_Dust,
                          quiet=1)  # ,maxiter=20000)# ,xtol = 1e-13,ftol=1e-13,gtol = 1e-13)
        p0_dust = m_dust.params

        if sp.Exfunc == 'blackbody':
            p0_dust = np.append(p0_dust, 0)  # JUST TO HAVE 3 PARAMETERS.. .0 FOR BETA DOESN'T MATTER == 1
        tempnew_dust = p0_dust[0]

    elif calcBB_bool:  # CALCULATE A BLACKBODY TO 2 FLUX POINTS -- MAINLY W3 AND W4

        alphaCon = 2 * _H * _CS ** 2
        gammaCon = _H * _CS / _KB

        lamArr = np.array([Lam_Excess[bandlow], Lam_Excess[bandhi]])
        flxArr = np.array([ExcessFlux[bandlow], ExcessFlux[bandhi]])

        flxArrErr = np.array([ExcessFluxerr[bandlow], ExcessFluxerr[bandhi]])

        res, FluxNormed, alpha = STools.calc_temp(flxArr, flxArrErr, SOBJ.mags4Dust,
                                           lwf.wfc.tempArr, lwf.wfc.kw_cor)
        chi2Dust = np.sum((res / flxArrErr) ** 2, axis=1)
        tempnew_dust = lwf.wfc.tempArr[np.where(chi2Dust.min() == chi2Dust)[0][0]]  # FIND NEW TEMPERATURES

        for bd in ExcessFlux.keys():
            kcor = 10 ** lwf.wfc.kw_cor['IP_fc%s' % bd](ma.log10(tempnew_dust))
            freal = ExcessFlux[bd] / kcor
            ExcessFlux[bd] = freal
        su2ea_dust = alpha[np.where(chi2Dust.min() == chi2Dust)][0]

        p0_dust = np.array([tempnew_dust])

    elif W3Excess_bool:
        alphaCon = 2 * _H * _CS ** 2
        gammaCon = _H * _CS / _KB
        l1band, l2band = 'W3', 'W4'
        # l1band, l2band ='W1','W2'
        lam1, lam2 = wave[l1band] * _ANG2CM, wave[l2band] * _ANG2CM
        lamArr = np.array([wave[l1band], wave[l2band]]) * _ANG2CM
        flxArr = np.array([ExcessFlux[l1band], ExcessFlux[l2band]])
        args = (lam1, lamArr, flxArr)
        tempnew_dust = bisect(STools.calcModTemp, 1., 10000., args=args)
        su2ea_dust = (ExcessFlux[l1band] * lam1 ** 5 / (_ANG2CM * alphaCon)) * \
                     (ma.exp(gammaCon / (lam1 * tempnew_dust)) - 1)
        beta = (gammaCon / (tempnew_dust * lam1)) * \
               (ma.exp(gammaCon / (lam1 * tempnew_dust)) / (ma.exp(gammaCon / (lam1 * tempnew_dust)) - 1.)) \
               - 5.
        p0_dust = np.array([tempnew_dust, 1])

        print 'beta = %.3f' % beta

    elif scaleBB_bool:  # IF ONLY SCALING TO A BB -- single flux point
        lamMin = max(Lam_Excess.values())
        # lamMin = wave['W3']
        tempnew_dust = STools.wienTEMP(lamMin, units='angstrom')
        # tempnew_dust = 272.
        for mv in SOBJ.mags4Dust:
            kcor = 10 ** lwf.wfc.kw_cor['IP_fc%s' % mv](ma.log10(tempnew_dust))
            BBDust_Flux[mv] = STools.blackbody(dust_lambda, np.array([tempnew_dust]),
                                               1, np.array([mv])) / kcor
            ExcessFlux[mv] = ExcessFlux[mv] / kcor

        BBDust_Fluxarr = np.array(at.dict2list(BBDust_Flux, SOBJ.mags4Dust))
        Excess_Flxarr = np.array(at.dict2list(ExcessFlux, SOBJ.mags4Dust))
        FluxNorm_dust = np.average(Excess_Flxarr / BBDust_Fluxarr,
                                   weights=1. / Excess_Flxerrarr)
        p0_dust = np.array([tempnew_dust, 1])
        # su2ea_dust = FluxNorm_dust
        su2ea_dust = Excess_Flxarr[0] / BBDust_Fluxarr[0]

    elif dontdraw_bool:
        FluxNorm_dust = -1
        Rad_dust = 1
        tempnew_dust = 1
        p0_dust = np.array([tempnew_dust, 0])
        print 'Not Drawing BB for %s' % star

    elif Nothing_bool:
        FluxNorm_dust = -1
        Rad_dust = 1
        tempnew_dust = -1
        p0_dust = np.array([tempnew_dust, 0, 0])
        print 'No Excess --> %s' % star

    else:
        print 'Something went wrong with %s' % star

    # =====================================================================================
    #                 CALCULATE & SET UP LINES TO PLOT
    #  =====================================================================================
    norm_star = (radius) ** 2 * su2ea2
    Lbol_star = (4 * ma.pi * (radius * _RSUN) ** 2) * np.trapz((yphot / _ANG2CM) / norm_star,
                                                                   xphot * _ANG2CM)
    if not dontdraw_bool:
        if W3Excess_bool:
            lightW4_line = STools.modifiedBB(dust_lambda, p0_dust, su2ea_dust, lam0=lam0, beta=beta)

        else:
            try:
                lightW4_line = eval('STools.' + sp.Exfunc + '(dust_lambda,p0_dust,su2ea_dust,lam0=lam0)')
            except NameError:
                lightW4_line = eval('STools.' + sp.Exfunc + '(dust_lambda,p0_dust,su2ea_dust)')

                #             if sp.Exfunc=='modifiedBB':
                pass
                # ftmp = STools.rsr_flux(STools.W2pband,dust_lambda,lightW4_line)
                # scale = ExcessFlux['W2']/ftmp
                #  MULITPLY ASSUMED RADIUS BY SCALING FACTOR TO SCALE FLUX FROM CALCULATED SPECTRUM
                #  TO W1. TYPICALLY THERE WILL BE ONLY ONE DOF, SO FIT WILL NTO BE ACCURATE.
                # p0_dust[1] = ma.sqrt(scale)*p0_dust[1]
                # lightW4_line = eval('STools.'+sp.Exfunc+'(dust_lambda,p0_dust,su2ea_dust,lam0=lam0)')

        fullspectrum = (yphot + lightW4_line) * xphot

        # tau = np.trapz(lightW4_line,dust_lambda)/np.trapz(yphot,xphot)
        #  try:
        #  dust_fullint = (p0_dust[1]**2*su2ea_dust*2*_c*con._kb*tempnew_dust) / (3*xphot.max()*_micron2cm)**3
        #  except:
        #  dust_fullint = (su2ea_dust*2*_c*con._kb*tempnew_dust) / (3*xphot.max()*_micron2cm)**3
        #  star_fullint = (p0[1]**2*su2ea2*2*_c*con._kb*tempnew)/(3*xphot.max()*_micron2cm)**3
        tau = np.trapz(lightW4_line, dust_lambda) / np.trapz(yphot, xphot)
        radius_dust = (278.3 / tempnew_dust) ** 2 * ma.sqrt(Lbol_star / _LSUN)

    else:
        indmax = np.where(yphot == yphot.max())[0][0]
        Fstar_max, wavestar_max = yphot[indmax], xphot[indmax]
        tempnew_dust = STools.wienTEMP(wave['W3'], units='angstrom')
        radius_dust = (278.3 / tempnew_dust) ** 2 * ma.sqrt(Lbol_star / _LSUN)
        PhotW3 = STools.rsr_flux(eval('STools.W3pband'), xphot, yphot)[0]
        ExcessFluxW3 = flux['W3_flux'] - PhotW3
        tau = (ExcessFluxW3 * wave['W3']) / (Fstar_max * wavestar_max)
        fullspectrum = (yphot + lightW4_line) * xphot

    # B60_lineraw = STools.blackbody(dust_lambda,np.array([Tdust60i,1]))
    # flux60_raw = STools.rsr_flux(STools.IRAS60pband,dust_lambda,B60_lineraw)
    # tmp = flux60i
    # scale = flux60i/flux60_raw
    # B60_lineNew = scale*B60_lineraw
    # f2460 = STools.rsr_flux(STools.W4pband,dust_lambda,B60_lineNew)
    # f2460kcor = 10**wfc.kw_cor['IP_fcW4'](ma.log10(Tdust60i))
    # f2460 = f2460/f2460kcor

    #  =================================================================
    #              BOLOMETRIC LUMINOSITIES
    #  =================================================================

    W3There = np.any(np.array(SOBJ.mags2use) == 'W3')
    W4There = np.any(np.array(SOBJ.mags2use) == 'W4')
    w4phot_cgs = STools.rsr_flux(STools.W4pband, xphot, yphot)[0]
    w3phot_cgs = STools.rsr_flux(STools.W3pband, xphot, yphot)[0]
    w4phot_mjy = STools.cgs2Jy(nu=STools.W4pband.isoFrequency(),
                               wave=STools.W4pband.isoWavelength(), flux=w4phot_cgs)
    w3phot_mjy = STools.cgs2Jy(nu=STools.W3pband.isoFrequency(),
                               wave=STools.W3pband.isoWavelength(), flux=w3phot_cgs)

    kcor_w3 = 10 ** lwf.wfc.kw_cor['IP_fcW3'](ma.log10(tempnew_dust))
    kcor_w4 = 10 ** lwf.wfc.kw_cor['IP_fcW4'](ma.log10(tempnew_dust))

    if W4There:
        w4f_cgs, w4fe_cgs = flux['W4_flux'], fluxerr['W4_flux']
        w4f_mjy, w4fe_mjy = STools.cgs2Jy(nu=STools.W4pband.isoFrequency(),
                                          wave=STools.W4pband.isoWavelength(),
                                          flux=w4f_cgs),\
                            STools.cgs2Jy(nu=STools.W4pband.isoFrequency(),
                                          wave=STools.W4pband.isoWavelength(),
                                          flux=w4fe_cgs)
        exw4corr_cgs = (w4f_cgs - w4phot_cgs) / kcor_w4
        exw4corr_mjy = (w4f_mjy - w4phot_mjy) / kcor_w4
        W4excessfraction = (w4f_cgs - w4phot_cgs) / w4f_cgs
        RelativeW4flux = w4f_cgs / w4phot_cgs
        e_RelaW4Flux = w4fe_cgs / w4phot_cgs
        try:
            w4upbool = str(Exbool_Dict['W4'])
        except KeyError:
            w4upbool = 'N/A'
    elif not W4There:
        w4f_cgs, w4fe_cgs = 0, 0
        w4f_mjy, w4fe_mjy = 0, 0
        w4upbool = 'N/A'
        exw4corr_cgs, exw4corr_mjy = 0, 0
        W4excessfraction, RelativeW4flux, e_RelaW4Flux = 0, 0, 0

    else:
        pass

    if W3There:
        w3f_cgs, w3fe_cgs = flux['W3_flux'], fluxerr['W3_flux']
        w3f_mjy, w3fe_mjy = STools.cgs2Jy(nu=STools.W3pband.isoFrequency(),
                                          wave=STools.W3pband.isoWavelength(),
                                          flux=w3f_cgs), \
                            STools.cgs2Jy(nu=STools.W3pband.isoFrequency(),
                                          wave=STools.W3pband.isoWavelength(),
                                          flux=w3fe_cgs)
        exw3corr_cgs = (w3f_cgs - w3phot_cgs) / kcor_w3
        exw3corr_mjy = (w3f_mjy - w3phot_mjy) / kcor_w3
        W3excessfraction = (w3f_cgs - w3phot_cgs) / w3f_cgs
        RelativeW3flux = w3f_cgs / w3phot_cgs
        e_RelaW3Flux = w3fe_cgs / w3phot_cgs

        try:
            w3upbool = str(Exbool_Dict['W3'])
        except KeyError:
            w3upbool = 'N/A'
    elif not W3There:
        w3f_cgs, w3fe_cgs = 0, 0
        w3f_mjy, w3fe_mjy = 0, 0
        w3upbool = 'N/A'
        exw3corr_cgs, exw3corr_mjy = 0, 0
        W3excessfraction, RelativeW3flux, e_RelaW3Flux = 0, 0, 0
    else:
        pass

    if sp.write2file and not p_args.nw:
        dat_str = '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t' + \
                  '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t' + \
                  '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s\n'

        f.write(dat_str % (star, w3f_cgs, w3fe_cgs, w3phot_cgs, w4f_cgs, w4fe_cgs, w4phot_cgs, \
                           exw3corr_cgs, exw4corr_cgs, W3excessfraction, W4excessfraction, \
                           RelativeW3flux, RelativeW4flux, e_RelaW4Flux, e_RelaW3Flux, \
                           w3f_mjy, w3fe_mjy, w3phot_mjy, w4f_mjy, w4fe_mjy, w4phot_mjy, \
                           w3upbool, w4upbool, kcor_w3, kcor_w4, \
                           ma.log10(Lbol_star / _LSUN), tempnew, radius, chi2, radius_dust, tempnew_dust, tau))

        try:
            print 'beta= ', p0_dust[2]
        except:
            pass
    else:
        pass

    #  =====================================================================================
    #                  PLOTTING
    #  =====================================================================================
    # PUT INTO ARRAY



    if sp.plot_any:
        fluxy = at.dict2list(flux, SOBJ.mags2use, '_flux')  # flux in erg s-1, cm-2 A-1
        fluxerry = at.dict2list(fluxerr, SOBJ.mags2use, '_flux')
        wavex = at.dict2list(wave, SOBJ.mags2use)  # WAVEX IS IN ANGSTROMS AT THIS POINT

        x = wavex
        y = fluxy * wavex
        x = wavex * _ANG2MICRON
        yerr = fluxerry * wavex

        if not sp.plot_single and sp.plot_grid:
            ptsize = 4
            ftsize = 12
            minorftsize = 12
            cps = 6

            if plot_i == int(sp.pROW * sp.pCOL):
                fmt = '.eps'
                ax2.set_xlabel(r'$\lambda (\mu m)$', fontsize=25,
                               family='sans-serif', labelpad=30)
                ax2.set_ylabel(r'$\lambda F_{\lambda} [erg\, s^{-1} cm^{-2}] $', fontsize=25,
                               family='sans-serif', labelpad=40)
                plt.subplots_adjust(left=.12, right=.96, bottom=.10,
                                    top=.93, hspace=0, wspace=0)
                save_name = os.path.join(os.getcwd(), 'DebrisDisks',
                                         'SEDGrids' + ran + '_' + str(gridnum) + fmt)
                if p_args.savf:
                    plt.savefig(save_name)
                    plt.clf()
                    plt.close()

                fig = plt.figure(figsize=plotsize)
                ax2 = fig.add_subplot(111)
                ax2.set_xticklabels([])
                ax2.set_yticklabels([])
                grid = Grid(fig, rect=111, nrows_ncols=(sp.pROW, sp.pCOL),
                            axes_pad=gridPad, label_mode='L')
                plot_i = 0
                gridnum += 1

            ax = grid[plot_i]
            ax.set_xlim([.2, xmax])
            ax.set_ylim([ylim_low, ylim_up])
            ax.loglog()
            PT.plot_setup(ax, majortick_size=10, minortick_size=3, \
                          ticklabel_fontsize=14, majortick_width=2, \
                          minortick_width=1, axes_linewidth=1.5)
            # PT.plot_setup(ax, majortick_size=10, minortickson=False,\
            #              ticklabel_fontsize=14, majortick_width=2,\
            #              axes_linewidth=1.5)
            ax.tick_params(axis='y', which='minor', left='off', right='off')
            ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
            plot_i += 1

            #  ==============================================================================================
            #                PLOT ANYTHING DEALING WITH PHOTOSPHERE
            #  ==============================================================================================
            if RJ_On:
                ax.plot((xphot * _ANG2MICRON), (yphot * xphot), 'b--', lw=1)
            else:
                plot_photosphere(ax, SOBJ.mags4Dust, xphot, yphot, wave,
                                 PhotDust_Flux, _ANG2MICRON, ptsize, lw=1)
            # ==============================================================================================
            #                PLOT ANYTHING DEALING WITH BLACKBODY (DUST)
            #  ==============================================================================================

            if not dontdraw_bool:
                if RJ_On:
                    ind_plotbb = np.where(dust_lambda >= 0)[0]
                else:
                    ind_plotbb = np.where(dust_lambda >= 2000.0)[0]
                dust_lambda, lightW4_line = dust_lambda[ind_plotbb], lightW4_line[ind_plotbb]
                plot_blackbody(ax, wave, ExcessFlux, dust_lambda,
                               lightW4_line, _ANG2MICRON, ptsize, 2)

            else:
                pass
            #  ==============================================================================================
            #                PLOT ANYTHING TO DEAL WITH OBSERVED DATA
            #  ==============================================================================================

            xi, yi = xphot, fullspectrum
            plot_observedData(ax, SOBJ.mags4Phot, xi, yi, wave,
                              flux, fluxerr,_ANG2MICRON, ptsize, cps, 1)

            #  ==============================================================================================
            #  plot_annotations(ax, star,spti,tempnew, tempnew_dust,beta,
            #                   p0_dust, W3Excess_bool,sp.Exfunc, ftsize, minorftsize)
            plot_annotationsW2(ax, star, spti, tempnew, tempnew_dust,beta,
                               p0_dust, W3Excess_bool, sp.Exfunc, ftsize,minorftsize)
            if not p_args.noshow:
                plt.show()


        elif not sp.plot_grid and sp.plot_single:
            save_name = os.path.join(os.getcwd(), 'DebrisDisks', str(star))
            ptsize = 9
            ftsize = 20
            minorftsize = 18
            cps = 11
            fig = plt.figure()
            ax = fig.add_subplot(111)
            print 'Plotting for %s' % star
            fig.set_size_inches(9.5, 7.5)
            ax.set_xlim([.2, xmax])
            ax.set_ylim([10 ** -14.5, ylim_up])
            ax = fig.add_subplot(111)
            PT.plot_setup(ax, minortickson=False)
            ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))
            ax.yaxis.set_minor_locator(MultipleLocator(1))
            ax.yaxis.set_major_locator(MultipleLocator(10))
            ax.tick_params(axis='y', which='minor', left='off', right='off')
            plt.loglog()

            if W3Excess_bool:
                ax.set_ylabel(r'$\lambda F_{\lambda} [erg\, s^{-1} cm^{-2}] $',
                              fontsize=35, family='sans-serif')
            else:
                ax.set_ylabel(r'$\lambda F_{\lambda} [erg\, s^{-1} cm^{-2}] $',
                              fontsize=35, family='sans-serif')

            ax.set_xlabel(r'$\lambda (\mu m)$', fontsize=35, family='sans-serif')
            plt.subplots_adjust(left=.19, right=.96, bottom=.15, top=.93)
            ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))

            #  ==============================================================================================
            #                PLOT ANYTHING DEALING WITH PHOTOSPHERE
            #  ==============================================================================================
            if RJ_On:
                ax.plot((xphot * _ANG2MICRON), (yphot * xphot), 'b--', lw=1)
            else:
                plot_photosphere(ax, SOBJ.mags4Dust, xphot, yphot, wave,
                                 PhotDust_Flux, _ANG2MICRON, ptsize, lw=1)

            # ==============================================================================================
            #                PLOT ANYTHING DEALING WITH BLACKBODY (DUST)
            #  ==============================================================================================

            if not dontdraw_bool:
                if RJ_On:
                    ind_plotbb = np.where(dust_lambda >= 0)[0]
                else:
                    ind_plotbb = np.where(dust_lambda >= 2000.0)[0]
                dust_lambda, lightW4_line = dust_lambda[ind_plotbb], lightW4_line[ind_plotbb]
                plot_blackbody(ax, wave, ExcessFlux, dust_lambda,
                               lightW4_line, _ANG2MICRON, ptsize, 2)

            else:
                pass
            #  ==============================================================================================
            #                PLOT ANYTHING TO DEAL WITH OBSERVED DATA
            #  ==============================================================================================

            xi, yi = xphot, fullspectrum
            plot_observedData(ax, mags4Phot0, xi, yi, wave, flux,
                              fluxerr, _ANG2MICRON, ptsize, cps, 1)

            #  ==============================================================================================

            plot_annotationsW2(ax, star, spti, tempnew, tempnew_dust, beta, p0_dust,
                               W3Excess_bool, sp.Exfunc, ftsize, minorftsize)

            if not p_args.noshow:
                plt.show()

            if p_args.savf:
                plt.savefig(save_name)
                plt.clf()
                plt.close()

if sp.write2file and not p_args.nw:
    f.close()
    print 'Data Written to: ', os.path.basename(file_write)
else:
    pass