#!/scisoft/bin/python

import directories
import numpy as np
from glob import glob
from scipy.interpolate import griddata, interp1d
import sed
import gzip
import os
import pyfits
from astro_tools import constants

import pdb
tDefault = (-np.inf,np.inf)
gDefault = (-np.inf,np.inf)
zDefault = (-0.0,0.0)
ang2cm = constants()._ang2cm

class gridModels:

    def __init__(self,tRange=tDefault,gRange=gDefault,zRange=zDefault,family='NextGen'):
        """Grid models object. initialize with the range of parameters wanted
           for interpolation, and desired model family. This starts by loading
           models and their parameters and saving them to object attributes.
        """

        if (family == 'NextGen'):
            self.base = os.path.join(directories.NextGen(),'Spectra')
        elif (family == 'ATLAS9'):
            self.base = directories.Atlas9()
        else:
            raise ValueError
        self.family = family
        self.params, self.models = self.read(tRange,gRange,zRange,True)

    def read(self,tRange,gRange,zRange,init=False):
        """Reads models from files based on family (Atlas9 or NextGen) in its
           native file structure. This module should always read at least eight
           models, and will return the parameters of each model and the models
           themselves (2D arrays of wavelength in microns and flux in
           ergs/s/cm^2/Ang).

           Ranges are given in tuples containing the desired max and min. The
           returned set will be inclusive of the bounds set -- that is, given
           a set of integer parameters, specifying bounds of (1.5,3.5) will
           return (1,2,3,4) and not (2,3).
        """

        if init:
            self.tRange = tRange
            self.gRange = gRange
            self.zRange = zRange

        tMin, tMax = tRange
        gMin, gMax = gRange
        zMin, zMax = zRange

        tArr = np.array([])
        gArr = np.array([])
        zArr = np.array([])
        mArr = []

        if (self.family == 'NextGen'):
            print 'Reading NextGen models...'

            zList = np.sort(np.array([os.path.basename(x).split('Z')[1] for x in glob(os.path.join(self.base,'Z-[0-9].[0-9]'))]).astype(np.float))
            zIndLo, zIndHi = max(np.searchsorted(zList,zMin,side='right')-1,0), min(np.searchsorted(zList,zMax,side='left')+1,len(zList))
            if (zIndHi-zIndLo < 2):
                if zIndLo > 0: zIndLo -= 1
                else: zIndHi += 1
            for Z in zList[zIndLo:zIndHi]:
                tList = np.unique(np.array([os.path.basename(x).split('-')[0].replace('lte','') for x in glob(os.path.join(self.base,'Z%.1f'%Z,'*'))]).astype(np.float)*100.0)
                tIndLo, tIndHi = max(np.searchsorted(tList,tMin,side='right')-1,0), min(np.searchsorted(tList,tMax,side='left')+1,len(tList))
                if (tIndHi-tIndLo < 2):
                    if tIndLo > 0: tIndLo -= 1
                    else: tIndHi += 1
                for T in tList[tIndLo:tIndHi]:
                    gList = np.sort(np.array([os.path.basename(x).split('-')[1] for x in glob(os.path.join(self.base,'Z%.1f'%Z,'lte%i*'%(T/100.0)))]).astype(np.float))
                    gIndLo, gIndHi = max(np.searchsorted(gList,gMin,side='right')-1,0), min(np.searchsorted(gList,gMax,side='left')+1,len(gList))
                    if (gIndHi-gIndLo < 2):
                        if gIndLo > 0: gIndLo -= 1
                        else: gIndHi += 1
                    for G in gList[gIndLo:gIndHi]:

                        #print 'Reading Z = %.1f, T = %i, G = %.1f'%(Z,T,G)

                        file = gzip.open(glob(os.path.join(self.base,'Z%.1f/lte%i-%.1f*'%(Z,T/100.0,G)))[0],'rb')

                        header = file.readline()
                        temp, grav, met = [float(x) for x in header.split()[:3]]
                        header = file.readline()
                        npoints = int(header.split()[0])

                        data = file.read()

                        w = np.array(data.split()[0:npoints]).astype(np.float)
                        f = np.array(data.split()[npoints:2*npoints]).astype(np.float)*ang2cm
                        model = np.array([w,f])

                        zArr = np.append(zArr,met)
                        tArr = np.append(tArr,temp)
                        gArr = np.append(gArr,grav)
                        mArr.append(model)

                        file.close()

        elif (self.family == 'Atlas9'):
            print 'Reading Atlas9 models...'

            signdict = {1.0:'p',0.0:'p',-1.0:'m'}
            zList = np.sort(np.array([os.path.basename(x).replace('kp','+').replace('km','-') for x in glob(os.path.join(self.base,'k*'))]).astype(np.float)/10.0)
            zIndLo, zIndHi = max(np.searchsorted(zList,zMin,side='right')-1,0), min(np.searchsorted(zList,zMax,side='left')+1,len(zList))
            if (zIndHi-zIndLo < 2):
                if zIndLo > 0: zIndLo -= 1
                else: zIndHi += 1
            for Z in zList[zIndLo:zIndHi]:
                tList = np.unique(np.array([os.path.basename(x).replace('.fits','').split('_')[1] for x in glob(os.path.join(self.base,'k%s%.2i'%(signdict[np.sign(Z)],abs(Z)*10.0),'*'))]).astype(np.float))
                tIndLo, tIndHi = max(np.searchsorted(tList,tMin,side='right')-1,0), min(np.searchsorted(tList,tMax,side='left')+1,len(tList))
                if (tIndHi-tIndLo < 2):
                    if tIndLo > 0: tIndLo -= 1
                    else: tIndHi += 1
                for T in tList[tIndLo:tIndHi]:
                    img = pyfits.open(os.path.join(self.base,'k%s%.2i'%(signdict[np.sign(Z)],abs(Z)*10.0),'k%s%.2i_%i.fits'%(signdict[np.sign(Z)],abs(Z)*10.0,T)))
                    thisT = img[0].header['TEFF']
                    gList = np.sort(np.array([x.replace('g','') for x in img[1].data.columns.names[1:]]).astype(np.float)/10.0)
                    gIndLo, gIndHi = max(np.searchsorted(gList,gMin,side='right')-1,0), min(np.searchsorted(gList,gMax,side='left')+1,len(gList))
                    if (gIndHi-gIndLo < 2):
                        if gIndLo > 0: gIndLo -= 1
                        else: gIndHi += 1
                    for G in gList[gIndLo:gIndHi]:

                        #print 'Reading Z = %.1f, T = %i, G = %.1f'%(Z,T,G)

                        w = np.array(img[1].data['WAVELENGTH'])
                        f = np.array(img[1].data['g%.2i'%(G*10.0)])
                        model = np.array([w,f])

                        zArr = np.append(zArr,img[0].header['LOG_Z'])
                        tArr = np.append(tArr,img[0].header['TEFF'])
                        gArr = np.append(gArr,G)
                        mArr.append(model)

                    img.close()

        points = np.transpose(np.array([tArr,gArr,zArr]))

        return points, mArr

    def convolve(self,bandpasses):
        """Given a list of bandpasses, this will convolve all models and store
           the resulting fluxes in self.fluxes. For each bandpass, an array of
           fluxes -- one entry per model -- will be generated, and then a list of
           these arrays -- one array per bandpass -- is saved. These are
           automatically used for interpolation.
        """
        fluxlist = []
        fluxdict = {}
        ST = sed.SEDTools()

        for bandpass in bandpasses:
            fluxes = []
            for wav, flx in self.models:
                fluxes = np.append(fluxes,ST.rsr_flux(bandpass,wav,flx))
            fluxlist.append(fluxes)
            fluxdict[bandpass] = fluxes

        self.fluxdict = fluxdict
        self.fluxes = fluxlist

    def evaluate(self,x,scale,T,G,Z=zDefault[0]):
        """Given a temperature, gravity, and metallicity, this returns an array
           of photometry values interpolated from griddata -- one for every
           bandpass given to the convolve routine.
        """
        if T < self.tRange[0]: T = self.tRange[0]
        elif T > self.tRange[1]: T = self.tRange[1]
        if G < self.gRange[0]: G = self.gRange[0]
        elif G > self.gRange[1]: G = self.gRange[1]
        if Z < self.zRange[0]: Z = self.zRange[0]
        elif Z > self.zRange[1]: Z = self.zRange[1]
        newvalues = np.array([])
        newpoints = np.array([[T,G,Z]])
        for table in self.fluxes:
            newvalue = griddata(self.params,table,newpoints,method='linear')[0]*scale
            newvalues = np.append(newvalues,newvalue)

        return newvalues

    def dicteval(self,bandpasses,scale,T,G,Z=zDefault[0]):
        """Same as above but uses a dictionary to return only some of the
           convolved fluxes.
        """
        if T < self.tRange[0]: T = self.tRange[0]
        elif T > self.tRange[1]: T = self.tRange[1]
        if G < self.gRange[0]: G = self.gRange[0]
        elif G > self.gRange[1]: G = self.gRange[1]
        if Z < self.zRange[0]: Z = self.zRange[0]
        elif Z > self.zRange[1]: Z = self.zRange[1]
        newvalues = np.array([])
        newpoints = np.array([[T,G,Z]])
        #print scale,T,G,Z
        for bandpass in bandpasses:
            newvalue = griddata(self.params,self.fluxdict[bandpass],newpoints,method='linear')[0]*scale
            newvalues = np.append(newvalues,newvalue)

        return newvalues

    def make(self,wavelengths,scale,T,G,Z=None):
        """Given an array of wavelength points and parameters for a model,
           this interpolates a new model from the nearest eight models using
           trilinear interpolation and returns an array of fluxes (in
           ergs/s/cm^2/Ang) to match the input wavelengths.
        """

        fluxes = []
        zFixed = False
        if Z == None:
            Z = zDefault[0]
            zFixed = True
        params, models = self.read((T,T),(G,G),(Z,Z))
        assert len(params) == 8
        if zFixed == True:
            models = [x for i,x in enumerate(models) if params[i][2]==zDefault[0]]
            params = [x for x in params if x[2]==zDefault[0]]
            assert len(params) == 4
        wav10 = np.log10(wavelengths)
        for i, (wav, flx) in enumerate(models):
            wav,flx = np.log10(wav), np.log10(flx)
            interp = interp1d(wav,flx,kind='linear')
            fluxes.append(np.power(10.0,interp(wav10)))

        tArr, gArr, zArr = np.transpose(params)

        tLo, tHi = min(tArr), max(tArr)
        gLo, gHi = min(gArr), max(gArr)
        zLo, zHi = min(zArr), max(zArr)
        tDiff, gDiff, zDiff = tHi-tLo, gHi-gLo, zHi-zLo

        tLo, tHi = (T-tLo)/tDiff, (tHi-T)/tDiff
        gLo, gHi = (G-gLo)/gDiff, (gHi-G)/gDiff
        if zFixed == False:
            zLo, zHi = (Z-zLo)/zDiff, (zHi-Z)/zDiff
        else:
            zLo, zHi = 1.0, 1.0

        if zFixed == True:
            flux = fluxes[0]*zHi*tHi*gHi+ \
                   fluxes[1]*zHi*tHi*gLo+ \
                   fluxes[2]*zHi*tLo*gHi+ \
                   fluxes[3]*zHi*tLo*gLo
        else:
            flux = fluxes[0]*zHi*tHi*gHi+ \
                   fluxes[1]*zHi*tHi*gLo+ \
                   fluxes[2]*zHi*tLo*gHi+ \
                   fluxes[3]*zHi*tLo*gLo+ \
                   fluxes[4]*zLo*tHi*gHi+ \
                   fluxes[5]*zLo*tHi*gLo+ \
                   fluxes[6]*zLo*tLo*gHi+ \
                   fluxes[7]*zLo*tLo*gLo

        return flux*scale
