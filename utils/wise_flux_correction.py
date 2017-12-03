import directories
import gzip
import math as ma
import matplotlib.pyplot as plt
import numpy as np
import operator
import os
import pdb
import scipy.integrate as sintp
import scipy.interpolate as intp
import sed
import sys
from glob import glob

import mosaic_tools as mt
from astro_tools import Constants
from readcol import*

con = Constants()
topdir = os.path.join(os.getcwd(), 'Interpolation_Files')
fRSR = os.path.join(topdir, 'RSR')
AT = mt.ArrayTools()
FT = mt.FittingTools()
ang2cm = con._ang2cm
ang2micron = con._ang2micron
micron2ang = 1./ang2micron
cm2ang = 1./ang2cm

class WISECorrect:
    """This module is designed to deal with corrections to WISE fluxes when integrating
    over the WISE bandpasses. The current WISE bandpasses use a zero point flux at the iso
    wavelength for a vega like spectrum. Any non-vega like spectrum integrated over the band-
    passes requires flux correction because of over estimation."""


    def __init__(self,xwave=None,bands=None):
        self.STools = sed.SEDTools()
        #IT'S SLOW SOMEWHERE BETWEEN HERE...
        #pdb.set_trace()
        self.base = os.path.join(directories.NextGen(),'Spectra')
        self.load_Vega()
        self.xwave=xwave
        if self.xwave is None:
            self.xwave = self.VegaModel[0]
            vega_interp = intp.interp1d(np.log10(self.xwave),np.log10(self.VegaModel[1]))
            self.xwave = np.arange(self.xwave[0],self.xwave[-1],10.)
            y = 10**vega_interp(np.log10(self.xwave))
            self.VegaModel = (self.xwave, y)
            y = None
            
        else: pass
        self.tempArr = np.logspace(ma.log10(20),ma.log10(1000),100)
        self.tarr,self.flxarr,self.wavearr = self.blackbody_spectrums(self.tempArr)
        self.kw_cor = {}
        for bd in bands:
            print bd
            bd_wave = eval('self.STools.%spband.isoWavelength()'%bd)
            vflux_bd = self.FindModelFlux(bd_wave,self.VegaModel)
            fluxUnSc_bd,fluxScaled_bd = self.blackbody_fluxes(self.tarr,self.wavearr,self.flxarr,bd)# convolved flux at band bd, and flux evaluated at isoBand wavelength for scaled blackbody to the
            #...AND HERE
            kcor_bd = vflux_bd/fluxScaled_bd #calculate correction factor array for band bd
            self.kw_cor['kcor_%s'%bd] = kcor_bd
            self.kw_cor['IP_F%s'%bd] = intp.interp1d(np.log10(self.tempArr),np.log10(fluxUnSc_bd)) #interpolator to evaluate corrected flux at any temperature
            self.kw_cor['IP_fc%s'%bd] = intp.interp1d(np.log10(self.tempArr),np.log10(kcor_bd)) #interpolator to evaluate correction factor at any temperature.
            self.kw_cor['F_cor_%s'%bd] = fluxUnSc_bd


    def load_Vega(self,T=10000,G=4.0,Z=0):
        """Loads a NextGen model spectrum similar in shape and form to a Vega Spectrum
        The default parameters are set to a T = 9600K, log(g) = 4.0 and solar metallicty.
        The Vega spectrum is known to be a metallicity of -0.5 dex, but NextGen does not
        have that metallicity at that temperature. Solar metallicity is sufficient for the
        purposes here. Works in conjunction with directories.py

        Parameters
        ----------
        T,G,Z: Scalar, temperature, log(g) and metallicity

        Returns
        -------
        Array of size [M x M] with wavelength and flux in units of angstroms and
        ergs/s/cm^2/Angstroms (Flam) for a model Vega spectrum (no dust)"""


        file = gzip.open(glob(os.path.join(self.base,'Z-%.1f/lte%i-%.1f*'%(Z,T/100.0,G)))[0],'rb')

        header = file.readline()
        temp, grav, met = [float(x) for x in header.split()[:3]]
        header = file.readline()
        npoints = int(header.split()[0])

        data = file.read()

        w = np.array(data.split()[0:npoints]).astype(np.float)
        f = np.array(data.split()[npoints:2*npoints]).astype(np.float)*ang2cm

        #newVega = np.array(FT.resample_model(w,f,w.min(),w.max()))

        self.VegaModel = np.array([w,f]) #newVega

        file.close()
        return

    def Vega_bandFlux(self,model,band):
        """Calculates the Filter convolved flux of a model using a particular bandpass.
        Parameters
        ----------
        model: 2xM array which contains two M-length arrays for the wavelength and flux of
                model spectrum
        band: string of the bandpass at which to calculate the integrated flux. List of
              acceptable bandpass names can be found in sed.py

        Returns
        -------
        flux: scalar flux of model convovled over bandpass"""

        wave,flux = model
        pband = eval('self.STools.'+band+'pband')

        flux = self.STools.rsr_flux(pband,wave,flux)[0]
        return flux

    def FindModelFlux(self, wave0,model,units='angstroms'):
        """find flux of model at specified wavelength. Assumes units of model are in angstroms
        and converts input units to angstroms
        model = [wave,flux]
        """

        if units == 'microns':
            wave *= micron2ang
        elif units == 'cm':
            wave *= cm2ang
        else: pass

        wav,flux = model
        ipolate = intp.interp1d(np.log10(wav),np.log10(flux))
        f_return = 10**ipolate(np.log10(wave0))
        return f_return

    #MIGHT BE GOING WRONG HERE. -- NOT SURE -- DO I NEED TO CONVOLVE?
    def ScaleModel2Vega(self,RawModelSpectrum,pband,VegaFlux=None):
        vflux = VegaFlux
        SpecRaw,WavRaw = RawModelSpectrum #WAVELNEGHT AND FLUX FOR UNSCALED MODEL SPECTRUM
        if vflux is None:
            vflux = self.VegaModel[1]
        try:
            Vinteg = self.STools.rsr_flux(pband,self.VegaModel[0],vflux)[0]
        except:
            pdb.set_trace()
        RawModelinteg = self.STools.rsr_flux(pband,WavRaw,SpecRaw)[0]

        scaledModel = (Vinteg/RawModelinteg)*np.array(SpecRaw)

        return scaledModel

    def blackbody_spectrums(self,tempArr,xwave=None):
        waveArr,fluxArr = [],[]
        xwave = xwave
        if xwave is None:
            xwave = self.xwave

        if np.isscalar(tempArr):
            tempArr = np.array([tempArr])
        for i in xrange(len(tempArr)):
            p0 = np.array([tempArr[i]])
            flux = self.STools.blackbody(xwave,p0)#band integrated flux
            waveArr.append(xwave)
            fluxArr.append(flux)

        return tempArr,np.array(fluxArr),np.array(waveArr)


    def blackbody_fluxes(self,tempArr,waveArr,fluxArr,band,convolved=False):
        """Gives you Real and Convolved (shifted/uncorrected) fluxes for blackbodiees at
        tempearatures in tempArr using Raw Surface fluxes in fluxArr"""
        
        print 'start bbflux stuff'
        fluxScaled = []
        fluxUnscaled = []
        #fluxConvolved = []
        
        if np.isscalar(tempArr):
            tempArr = np.array([tempArr])
        else: pass
        #vega's convolved flux in the band= band

        pband = eval('self.STools.'+band+'pband')
        bandWave = pband.isoWavelength()


        for i in xrange(len(tempArr)):
            wavei = waveArr[i]
            fluxi = fluxArr[i]#blackbody spectrum at temperature tempArr[i]

            scaledFlux_i = self.ScaleModel2Vega([fluxi,wavei],pband) #Blackbody spectrum scaled to vega
            fluxScaledi = self.FindModelFlux(bandWave,[wavei,scaledFlux_i]) #blackbody scaled flux evaluated at wavei
            if convolved:
                fluxUnSc = self.STools.rsr_flux(pband,wavei,fluxi)[0]#band integrated flux
            else:
                fluxUnSc = self.FindModelFlux(bandWave,[wavei,fluxi])

            fluxScaled.append(fluxScaledi)
            fluxUnscaled.append(fluxUnSc)
        print 'end bbflux stuff'
        return np.array(fluxUnscaled), np.array(fluxScaled)

    def calc_CorrectedBlackbody(self,wave, p0,kw=None,scale=1):
        kwarg = kw
        temp,A0 = p0
        temp0 = log10(temp*100)
        su2ea = (A0**2)*scale
        newfluxes = []
        A0 *= scale
        for band in wave:
            this_kcor = 10**kwarg['IP_fc%s'%band](temp0)
            this_flux = this_kcor*10**(kwarg['IP_F%s'%band](temp0))
            this_flux *= su2ea

            newfluxes.append(this_flux)

        return newfluxes
