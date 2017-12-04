ALL RSR'S ARE STRICTLY THAT.. ALL LAMBDA*RSR has been divided out

Header portion in each RSR File: #!N0[0,f]/ N1[1,f]/ N2[2,f]/ N3[3,f]/ N4[4,f]/ N5[5,f]/

N0 = Isophotal wavelength (Angstroms)
N1 = vega zero point flux in erg/s/cm^2/Angstrom
N2 = uncertainty in vega zero point flux in erg/s/cm^2/Angstrom
N3 = Isophotal frequency (Hz)
N4 = vega zero point flux in erg/s/cm^2/Hz
N5 = uncertainty in vega zero point flux in erg/s/cm^2/Hz


AB mag conversion:

http://www.astro.utoronto.ca/~patton/astro/mags.html

http://www.aerith.net/astro/color_conversion.html

http://casa.colorado.edu/~ginsbura/filtersets.htm__

2MASS:

x: wavelength (microns)
y: relative spectral response

U,B,V,I:

x: wavelength (angrstoms)
y: Relative spectral response

Johnson zero points from http://www.stsci.edu/hst/observatory/documents/isrs/scs8.rev.pdf

Taken from
http://spiff.rit.edu/classes/phys440/lectures/filters/filters.html

WISE:
w1,w2,w3,w4
x:wavelength (microns)
y:rsr (not normalized)
# Relative Spectral Response per erg (Equal-Energy); Jarrett et al. 2011
taken from http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#WISEZMA

The following RSRs are modified RSRs from Brown+2014
The monochromatic flux densities are changed based on
the new vega AB mag = 6.66

    truncated W4:
    Identical to lab RSR except it has zero transmission lambda<21.1 microns

    stretched W4: 
    Identiacal to lab except all the wavelengths have been revised upward by 3.3% (i.e., ∆λ = 0.033λ)
    
    


IRAS:
60
x:wavelength(Angstrom)
y:rsr




SPITZER: 
MIPS: FROM http://irsa.ipac.caltech.edu/data/SPITZER/docs/mips/calibrationfiles/spectralresponse/
Lambda in microns
wavelength is the one gien in MIPS handbook


HERSCHEL:
from SVO profile service
pivot wavelength used
angstroms

AKARI:
from SVO profile service
wavelength in angstrom
Using lam_eff for wavelength


NIRC2:
http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=Keck/NIRC2.H&&mode=browse&gname=Keck&gname2=NIRC2#filter

filter data from NIRC2 website. 
central wavelength and zero point is from SVO 
